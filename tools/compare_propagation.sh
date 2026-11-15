#!/usr/bin/env bash
# Checks that a change which may alter rounding (summation order, FFT algorithm or library)
# keeps the results within a tolerance: builds a base commit and the working tree, runs a short
# propagation (50 time steps of the hard wall case, averages after every analysis step) with
# both, and compares averages.dat column by column with tools/compare/compare_averages.py
# (largest relative difference; a value that is not finite fails). For changes that must not
# alter any bit use tools/check_bitwise_refactor.sh.
#
# Usage: tools/compare_propagation.sh [base-ref] [tolerance]   (defaults: HEAD, 1e-13)
#   tolerance: a finite number >= 0, checked before the builds.
# Output: the largest relative difference per column and one "OK" line, or "DIFFERENT: ...".
#   Exit 0 within the tolerance, 1 different or error, 2 usage.
# Environment: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (nlohmann_json sources for an offline
#   configure; default <repo>/build/_deps/nlohmann_json-src, cloned from GitHub if that is
#   missing), OMP_NUM_THREADS (threads of the runs), TMPDIR (scratch builds),
#   CMAKE_BUILD_PARALLEL_LEVEL (build jobs, read by cmake); the builds: tools/scratch_build.sh.
# Run time: about 1 min (two Release builds of spodvr and two runs; 58 s and 42 s measured
#   2026-11-15, offline, CMAKE_BUILD_PARALLEL_LEVEL=4 and OMP_NUM_THREADS=4 on 8 of the 16
#   hardware threads).
# Docs: README.md, section "Testing".
set -euo pipefail

usage() {
  echo "usage: $0 [base-ref] [tolerance]   (tolerance: a finite number >= 0, default 1e-13)" >&2
  exit 2
}
[[ $# -le 2 ]] || usage
base_ref="${1:-HEAD}"
tolerance="${2:-1e-13}"
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
compare_dir="$repo/tools/compare"
# The tolerance rule of compare_averages.py, applied before the builds rather than after them.
python3 -I -B -c 'import sys; sys.path.insert(0, sys.argv[1]); from compare_averages import parse_tolerance
parse_tolerance(sys.argv[2])' "$compare_dir" "$tolerance" 2>/dev/null || usage
if ! git -C "$repo" rev-parse --verify --quiet "$base_ref^{commit}" >/dev/null; then
  echo "error: '$base_ref' is not a commit in $repo; give a branch, tag or commit hash" >&2
  exit 1
fi
source "$repo/tools/scratch_build.sh"
scratch="$(mktemp -d)"
trap 'git -C "$repo" worktree remove --force "$scratch/base-src" >/dev/null 2>&1 || true; rm -rf "$scratch"' EXIT

echo "Building base ($base_ref) and working tree ..."
git -C "$repo" worktree add --detach "$scratch/base-src" "$base_ref" >/dev/null 2>&1
build_spodvr "$scratch/base-src" "$scratch/base-build"
build_spodvr "$repo" "$scratch/new-build"

python3 -I -B -c '
import json, sys
config = json.load(open(sys.argv[1]))
config.update({"tEnd": 100.0, "nStep": 50, "dtAve": 10.0})
json.dump(config, open(sys.argv[2], "w"), indent=2)
' "$repo/tests/reference/hardwall_short/input.json" "$scratch/propagation.json"

echo "Running the short propagation ..."
for version in base new; do
  "$scratch/$version-build/spodvr" "$scratch/propagation.json" --no-color --output-dir "$scratch/$version" >/dev/null
done

python3 -I -B "$compare_dir/compare_averages.py" "$scratch/base/averages.dat" "$scratch/new/averages.dat" "$tolerance"
