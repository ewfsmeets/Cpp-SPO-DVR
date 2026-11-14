#!/usr/bin/env bash
# Mutation check of the CRP library tests (docs/crp_library.md "Mutation check"): applies each
# mutation of mutations.py to a git worktree of HEAD, rebuilds the library tests there (with the
# test of its file reader, tests/fortran_text) and runs ctest -R "CrpLib|FortranText"; every
# mutation must make at least one test fail, except those marked [SAFETY MARGIN] (a redundant
# guard: removing it alone cannot change a result). The source tree is not touched.
# The patterns match the code as it was when a mutation was added: update them when that code
# changes (a pattern that no longer matches is reported as "not applied").
#
# Usage: tools/crp_mutation_check/run.sh [scratch_dir [index ...]]
#   scratch_dir  default: a new mktemp -d directory;  index: only these mutations (default: all)
# Environment: OMP_NUM_THREADS (build jobs and test threads, default 8);
#   FETCHCONTENT_SOURCE_DIR_GOOGLETEST and FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (sources for an
#   offline configure; default <repo>/build/_deps/googletest-src and .../nlohmann_json-src, cloned
#   from GitHub if those are missing).
# Run time: 15-16 min (28 mutations at 4 threads; three runs on 2026-11-14).
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
scratch="${1:-$(mktemp -d)}"
shift $(($# > 0 ? 1 : 0))
selected=("$@")
# Sources for an offline configure: from the environment, else the main build's downloads if present.
gtest_source="${FETCHCONTENT_SOURCE_DIR_GOOGLETEST:-${repo}/build/_deps/googletest-src}"
json_source="${FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON:-${repo}/build/_deps/nlohmann_json-src}"
if [[ -n "${FETCHCONTENT_SOURCE_DIR_GOOGLETEST:-}" && ! -d "${gtest_source}" ]]; then
  echo "error: FETCHCONTENT_SOURCE_DIR_GOOGLETEST='${gtest_source}' is not a directory; point it at the googletest sources or unset it" >&2
  exit 1
fi
if [[ -n "${FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON:-}" && ! -d "${json_source}" ]]; then
  echo "error: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON='${json_source}' is not a directory; point it at the nlohmann_json sources or unset it" >&2
  exit 1
fi
source_options=()
[[ -d "${gtest_source}" ]] && source_options+=(-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST="${gtest_source}")
[[ -d "${json_source}" ]] && source_options+=(-DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON="${json_source}")
wt="${scratch}/worktree"
mkdir -p "${scratch}"
[[ -d "${wt}" ]] || git -C "${repo}" worktree add --detach "${wt}" HEAD > /dev/null
if [[ ! -f "${wt}/build/build.ninja" ]]; then
  cmake -S "${wt}" -B "${wt}/build" -G Ninja -DCMAKE_BUILD_TYPE=Release "${source_options[@]}" \
    > "${scratch}/configure.log" 2>&1 || { echo "configure failed: ${scratch}/configure.log"; exit 1; }
fi
# The library's tests, and those of its file reader (tests/fortran_text, tests FortranText.*)
crplib_targets=$(cd "${wt}/tests/crplib" && ls test_*.cpp | sed 's/\.cpp$//; s/^/crplib_/')
reader_targets=$(cd "${wt}/tests/fortran_text" && ls test_*.cpp | sed 's/\.cpp$//; s/^/fortran_text_/')
targets="${crplib_targets} ${reader_targets}"
threads="${OMP_NUM_THREADS:-8}"  # build jobs and OpenMP threads (at most 8 on this machine)
build() { cmake --build "${wt}/build" -j "${threads}" --target ${targets} > "${scratch}/build.log" 2>&1; }
run() { OMP_NUM_THREADS="${threads}" ctest --test-dir "${wt}/build" -R "CrpLib|FortranText" > "${scratch}/ctest.log" 2>&1; }

echo "baseline:"
build || { echo "  build failed: ${scratch}/build.log"; exit 1; }
run && echo "  all pass" || { echo "  baseline fails: ${scratch}/ctest.log"; exit 1; }
survived=0
n=$(python3 -I "${here}/mutations.py" list)
indices=("${selected[@]}")
[[ ${#indices[@]} -gt 0 ]] || indices=($(seq 0 $((n - 1))))
for i in "${indices[@]}"; do
  git -C "${wt}" checkout -q -- src
  name=$(python3 -I "${here}/mutations.py" "${wt}" "${i}") || { echo "${i}: ${name} (not applied)"; continue; }
  if ! build; then echo "${i} ${name}: BUILD FAILED"; continue; fi
  if run; then
    if [[ "${name}" == *"[SAFETY MARGIN]"* ]]; then
      echo "${i} ${name}: survived (expected: a redundant margin, no test can tell)"
    else
      echo "${i} ${name}: SURVIVED"
      survived=$((survived + 1))
    fi
  else
    echo "${i} ${name}: killed by $(grep -cE '\(Failed\)' "${scratch}/ctest.log") tests," \
      "e.g. $(grep -E '\(Failed\)' "${scratch}/ctest.log" | head -1 | sed -E 's/.*(CrpLib|FortranText)\.//; s/ .*//')"
  fi
done
git -C "${wt}" checkout -q -- src
git -C "${repo}" worktree remove --force "${wt}"
exit $((survived > 0))
