#!/usr/bin/env bash
# Time-loop additions must cost nothing when off (docs/CONTRIBUTING.md, "Additions to the time loop"): runs two
# existing spodvr builds alternately (A B B A A B ..., R times each) on the same input without plot keys and
# compares their "Propagation time" and host peak memory (tools/compare/compare_speed.py: at most
# max(1 %, spread) slower, at most 16 MiB more; GPU memory is not measured). It builds nothing.
#   CPU (default): the medium grid of tools/benchmark/run_benchmark.sh (keep the two in step), 8 threads,
#   50 time steps.
#   --cuda: the DZ input (inputs/thesis_fig5_3_dz) on the GPU, 40 time steps.
#   Both: averages every analysis step; the number of steps must make whole analysis steps (CPU: a multiple
#   of 5, GPU: of 8).
#
# Usage: tools/benchmark/compare_speed.sh --base-build DIR --new-build DIR [--cuda] [--steps N]
#          [--repetitions R] [--logs DIR]
#   DIR: a build directory holding spodvr; R: runs per build, at least 2 (default 5); --logs: an empty or new
#   directory that keeps the run logs (default: a scratch directory, removed at the end).
# Exit status: 0 within the bounds, 1 slower or more memory (or a run failed), 2 bad arguments.
# Environment: none read. Run time: CPU about 10 x 40 s, GPU about 10 x 2 min (measured with R = 5).
# Disk: each run writes its checkpoint (DZ: about 6 GiB) into a scratch directory under /var/tmp (on disk; /tmp
#   may be a RAM-backed tmpfs too small for it), removed right after the run (the log is kept elsewhere).
# Docs: docs/CONTRIBUTING.md, section "Additions to the time loop".
set -euo pipefail

usage_error() { echo "error: $1; see the usage in $0" >&2; exit 2; }
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
base="" new="" cuda=0 steps="" repetitions=5 logs=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cuda) cuda=1; shift; continue ;;
    --base-build|--new-build|--steps|--repetitions|--logs) [[ $# -ge 2 ]] || usage_error "$1 needs a value" ;;
    *) usage_error "unknown argument '$1'" ;;
  esac
  case "$1" in
    --base-build) base="$2" ;;
    --new-build) new="$2" ;;
    --steps) steps="$2" ;;
    --repetitions) repetitions="$2" ;;
    --logs) logs="$2" ;;
  esac
  shift 2
done
[[ -n $base && -n $new ]] || usage_error "give both --base-build and --new-build"
for build in "$base" "$new"; do
  [[ -x "$build/spodvr" ]] || usage_error "'$build/spodvr' is not a program; give a build directory"
done
base="$(realpath "$base")" new="$(realpath "$new")"
[[ $base != "$new" ]] || usage_error "--base-build and --new-build are the same directory"
if cmp -s "$base/spodvr" "$new/spodvr"; then echo "warning: the two spodvr programs are identical" >&2; fi
echo "base: $base/spodvr, $("$base/spodvr" --version 2>&1 || true)"
echo "new:  $new/spodvr, $("$new/spodvr" --version 2>&1 || true)"
[[ $repetitions =~ ^[0-9]+$ && $repetitions -ge 2 ]] || usage_error "--repetitions must be a whole number >= 2"
[[ -z $steps || $steps =~ ^[1-9][0-9]*$ ]] || usage_error "--steps must be a positive whole number"
scratch="$(mktemp -d --tmpdir=/var/tmp compare_speed.XXXXXX)"
trap 'rm -rf "$scratch"' EXIT
if [[ -n $logs ]]; then
  mkdir -p "$logs"
  [[ -z "$(ls -A "$logs")" ]] || usage_error "--logs '$logs' is not empty (older logs would be read too)"
  logs="$(realpath "$logs")"
else
  logs="$scratch/logs"; mkdir -p "$logs"
fi

# The input: the medium grid (CPU) or the DZ grid (GPU), cut to the steps, averages every analysis step.
python3 -I -B -c '
import json, sys
repo, out, cuda, steps = sys.argv[1], sys.argv[2], sys.argv[3] == "1", sys.argv[4]
if cuda:
    config = json.load(open(repo + "/inputs/thesis_fig5_3_dz/input.json"))
    steps = int(steps or 40)
else:  # the medium grid of tools/benchmark/run_benchmark.sh
    config = json.load(open(repo + "/tests/reference/hardwall_short/input.json"))
    config.update({"nX": 16, "nY": 16, "JMax": 8, "MJMax": 8})
    steps = int(steps or 50)
per_analysis = round(config["dtAna"] / config["dtStep"])
if steps % per_analysis != 0:
    sys.exit(f"error: --steps {steps} is not a whole number of analysis steps ({per_analysis} time steps each)")
config.update({"tEnd": config["dtStep"] * steps, "nStep": steps, "dtAve": config["dtAna"]})
json.dump(config, open(out, "w"), indent=1)
' "$repo" "$scratch/input.json" "$cuda" "$steps" || exit 2
device=()
[[ $cuda == 1 ]] && device=(--device cuda)

for k in $(seq 1 "$repetitions"); do
  order="base new"; (( k % 2 == 0 )) && order="new base"  # alternate which build runs first
  for which in $order; do
    build="$base"; [[ $which == new ]] && build="$new"
    log="$logs/$which-$k.log"
    OMP_NUM_THREADS=8 /usr/bin/time -v "$build/spodvr" "$scratch/input.json" --threads 8 "${device[@]}" --no-color \
      --output-dir "$scratch/out-$which" > "$log" 2>&1 || { echo "error: run $which $k failed; see $log" >&2; exit 1; }
    rm -rf -- "${scratch:?}/out-$which"  # only the log is needed; one checkpoint on disk at a time
    echo "$which $k: $(grep -E 'Propagation time' "$log" | sed 's/^ *//')"
  done
done
python3 -I -B "$repo/tools/compare/compare_speed.py" --base "$logs"/base-*.log --new "$logs"/new-*.log
