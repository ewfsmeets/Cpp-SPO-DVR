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
# Output: one line per mutation. Exit 0 every mutation caught, 1 one survived or the baseline
#   build or tests failed, 2 bad arguments.
# Environment: OMP_NUM_THREADS (build jobs and test threads, default 8);
#   FETCHCONTENT_SOURCE_DIR_GOOGLETEST and FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (sources for an
#   offline configure; default <repo>/build/_deps/googletest-src and .../nlohmann_json-src, cloned
#   from GitHub if those are missing; tools/scratch_build.sh offline_sources). Needs the checkout's
#   Python environment (`uv sync`, README "Dependencies"): the scratch worktree's tests use its .venv.
# Run time: 15-16 min (28 mutations at 4 threads; three runs on 2026-11-14).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
count="$(python3 -I -B "${here}/mutations.py" list)"

usage() {
  echo "usage: run.sh [scratch_dir [index ...]]   (index: 0 to $((count - 1)))" >&2
  exit 2
}
scratch="${1:-$(mktemp -d)}"
shift $(($# > 0 ? 1 : 0))
indices=("$@")
for i in "${indices[@]}"; do
  [[ "${i}" =~ ^(0|[1-9][0-9]*)$ && "${i}" -lt "${count}" ]] || usage
done
((${#indices[@]} > 0)) || mapfile -t indices < <(seq 0 $((count - 1)))

# shellcheck source=../scratch_build.sh
source "${repo}/tools/scratch_build.sh"
offline_sources GOOGLETEST NLOHMANN_JSON
wt="${scratch}/worktree"
mkdir -p "${scratch}"
[[ -d "${wt}" ]] || git -C "${repo}" worktree add --detach "${wt}" HEAD > /dev/null
if [[ ! -f "${wt}/build/build.ninja" ]] &&
  ! cmake -S "${wt}" -B "${wt}/build" -G Ninja -DCMAKE_BUILD_TYPE=Release "${source_options[@]}" \
    -DPython3_EXECUTABLE="${repo}/.venv/bin/python" \
    > "${scratch}/configure.log" 2>&1; then
  echo "error: the configure failed; see ${scratch}/configure.log" >&2
  exit 1
fi
# The library's tests, and those of its file reader (tests/fortran_text, tests FortranText.*).
targets=()
for test_source in "${wt}"/tests/crplib/test_*.cpp; do targets+=("crplib_$(basename "${test_source}" .cpp)"); done
for test_source in "${wt}"/tests/fortran_text/test_*.cpp; do targets+=("fortran_text_$(basename "${test_source}" .cpp)"); done
threads="${OMP_NUM_THREADS:-8}"  # build jobs and OpenMP threads (at most 8 on this machine)
build() { cmake --build "${wt}/build" -j "${threads}" --target "${targets[@]}" > "${scratch}/build.log" 2>&1; }
run() { OMP_NUM_THREADS="${threads}" ctest --test-dir "${wt}/build" -R "CrpLib|FortranText" > "${scratch}/ctest.log" 2>&1; }

echo "baseline:"
build || { echo "error: the baseline build failed; see ${scratch}/build.log" >&2; exit 1; }
run || { echo "error: the baseline tests fail; see ${scratch}/ctest.log" >&2; exit 1; }
echo "  all pass"
survived=0
for i in "${indices[@]}"; do
  git -C "${wt}" checkout -q -- src
  name=$(python3 -I -B "${here}/mutations.py" "${wt}" "${i}") || { echo "${i}: ${name} (not applied)"; continue; }
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
