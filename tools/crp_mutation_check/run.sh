#!/usr/bin/env bash
# Mutation check of the CRP library tests (docs/crp_library.md "Mutation check"): applies each
# mutation of mutations.py to a git worktree of HEAD, rebuilds the library tests there and runs
# ctest -R CrpLib; every mutation must make at least one test fail. The source tree is not touched.
# The patterns match the code as it was when a mutation was added: update them when that code
# changes (a pattern that no longer matches is reported as "not applied").
#
# Usage: tools/crp_mutation_check/run.sh [scratch_dir]   (default: a new mktemp -d directory)
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
scratch="${1:-$(mktemp -d)}"
wt="${scratch}/worktree"
mkdir -p "${scratch}"
[[ -d "${wt}" ]] || git -C "${repo}" worktree add --detach "${wt}" HEAD > /dev/null
gtest_source="${repo}/build/_deps/googletest-src"  # reuse the main build's download if present
gtest_option=()
[[ -d "${gtest_source}" ]] && gtest_option=(-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST="${gtest_source}")
if [[ ! -f "${wt}/build/build.ninja" ]]; then
  cmake -S "${wt}" -B "${wt}/build" -G Ninja -DCMAKE_BUILD_TYPE=Release "${gtest_option[@]}" \
    > "${scratch}/configure.log" 2>&1 || { echo "configure failed: ${scratch}/configure.log"; exit 1; }
fi
targets=$(cd "${wt}/tests/crplib" && ls test_*.cpp | sed 's/\.cpp$//; s/^/crplib_/')
build() { cmake --build "${wt}/build" -j 8 --target ${targets} > "${scratch}/build.log" 2>&1; }
run() { OMP_NUM_THREADS=8 ctest --test-dir "${wt}/build" -R CrpLib > "${scratch}/ctest.log" 2>&1; }

echo "baseline:"
build || { echo "  build failed: ${scratch}/build.log"; exit 1; }
run && echo "  all pass" || { echo "  baseline fails: ${scratch}/ctest.log"; exit 1; }
survived=0
n=$(python3 -I "${here}/mutations.py" list)
for ((i = 0; i < n; ++i)); do
  git -C "${wt}" checkout -q -- src
  name=$(python3 -I "${here}/mutations.py" "${wt}" "${i}") || { echo "${i}: ${name} (not applied)"; continue; }
  if ! build; then echo "${i} ${name}: BUILD FAILED"; continue; fi
  if run; then
    echo "${i} ${name}: SURVIVED"
    survived=$((survived + 1))
  else
    echo "${i} ${name}: killed by $(grep -cE '\(Failed\)' "${scratch}/ctest.log") tests," \
      "e.g. $(grep -E '\(Failed\)' "${scratch}/ctest.log" | head -1 | sed 's/.*CrpLib\.//; s/ .*//')"
  fi
done
git -C "${wt}" checkout -q -- src
git -C "${repo}" worktree remove --force "${wt}"
exit $((survived > 0))
