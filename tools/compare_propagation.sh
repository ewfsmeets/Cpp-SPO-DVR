#!/usr/bin/env bash
# Checks that a change which may alter rounding (summation order, FFT algorithm or library)
# keeps the results within a tolerance: builds a base commit and the working tree, runs a short
# propagation (50 time steps of the hard wall case, averages after every analysis step) with
# both, and compares averages.dat column by column (largest relative difference). For changes
# that must not alter any bit use tools/check_bitwise_refactor.sh instead.
#
# Usage: tools/compare_propagation.sh [base-ref] [tolerance]   (defaults: HEAD, 1e-13)
set -euo pipefail

base_ref="${1:-HEAD}"
tolerance="${2:-1e-13}"
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
scratch="$(mktemp -d)"
trap 'git -C "$repo" worktree remove --force "$scratch/base-src" >/dev/null 2>&1 || true; rm -rf "$scratch"' EXIT

configure_and_build() {  # <source dir> <build dir>
  local json_src="$repo/build/_deps/nlohmann_json-src"
  local reuse=()
  [[ -d "$json_src" ]] && reuse=(-DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON="$json_src")
  cmake -S "$1" -B "$2" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_BUILD_TESTS=OFF -DSPODVR_BUILD_TOOLS=OFF \
    "${reuse[@]}" >/dev/null
  cmake --build "$2" --target spodvr >/dev/null
}

echo "Building base ($base_ref) and working tree ..."
git -C "$repo" worktree add --detach "$scratch/base-src" "$base_ref" >/dev/null 2>&1
configure_and_build "$scratch/base-src" "$scratch/base-build"
configure_and_build "$repo" "$scratch/new-build"

python3 -c '
import json, sys
config = json.load(open(sys.argv[1]))
config.update({"tEnd": 100.0, "nStep": 50, "dtAve": 10.0})
json.dump(config, open(sys.argv[2], "w"), indent=2)
' "$repo/tests/reference/hardwall_short/input.json" "$scratch/propagation.json"

echo "Running the short propagation ..."
for version in base new; do
  "$scratch/$version-build/spodvr" "$scratch/propagation.json" --no-color --output-dir "$scratch/$version" >/dev/null
done

python3 - "$scratch/base/averages.dat" "$scratch/new/averages.dat" "$tolerance" <<'PY'
import sys
def rows(path):
    return [[float(x) for x in line.split()] for line in open(path) if not line.startswith("#")]
base, new, tolerance = rows(sys.argv[1]), rows(sys.argv[2]), float(sys.argv[3])
if len(base) != len(new):
    sys.exit(f"DIFFERENT: {len(base)} vs {len(new)} rows")
names = ["time", "norm", "<Z>", "<r>", "rest norm"]
worst = 0.0
for column, name in enumerate(names):
    d = max(abs(a[column] - b[column]) / max(abs(a[column]), abs(b[column]), 1e-300) for a, b in zip(base, new))
    worst = max(worst, d)
    print(f"  {name:10s} max relative difference {d:.3e}")
if worst > tolerance:
    sys.exit(f"DIFFERENT: {worst:.3e} > tolerance {tolerance:g}")
print(f"OK: within {tolerance:g} of the base")
PY
