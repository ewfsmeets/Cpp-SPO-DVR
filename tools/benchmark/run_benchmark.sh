#!/usr/bin/env bash
# Benchmark on a realistic-size grid ("medium": the hard wall short case with nx = ny = 16,
# JMax = MJMax = 8, i.e. 45 rotational states and a 1.8 GB wave function; ~5 GB of memory with
# the precomputed exp(-iV dt) table). The tiny test grids are not representative for performance.
#
# Runs the program twice (timings on WSL/laptops vary; compare the printed "Propagation time",
# not the wall time) and the per-operator timing tool once.
# Rebuild the build directory first: comparing a stale binary gives wrong conclusions.
#
# Usage: tools/benchmark/run_benchmark.sh [build-dir] [threads] [steps]
#   build-dir  default: build;  threads: default all;  steps: default 20
set -euo pipefail

repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
build="${1:-build}"
threads="${2:-$(nproc)}"
steps="${3:-20}"
scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT

python3 -c '
import json, sys
config = json.load(open(sys.argv[1]))
steps = int(sys.argv[3])
config.update({"nX": 16, "nY": 16, "JMax": 8, "MJMax": 8,
               "tEnd": 2.0 * steps, "nStep": steps, "dtAve": 10.0})
json.dump(config, open(sys.argv[2], "w"), indent=1)
' "$repo/tests/reference/hardwall_short/input.json" "$scratch/medium.json" "$steps"

echo "Benchmark: $build, $threads threads, $steps time steps (medium grid)"
for run in 1 2; do
  "$repo/$build/spodvr" "$scratch/medium.json" --threads "$threads" --no-color --output-dir "$scratch/out" \
    | grep -E "Propagation time" | sed "s/^ */  run $run: /"
done
if [[ -x "$repo/$build/tools/spodvr_operator_timing" ]]; then
  echo "Per operator:"
  OMP_NUM_THREADS="$threads" "$repo/$build/tools/spodvr_operator_timing" "$scratch/medium.json" 3
fi
