#!/usr/bin/env bash
# Checks that a refactor did not change any numerical result: builds a base commit and the working
# tree, runs the potential-only hard wall case and a short propagation (50 time steps of the hard
# wall case) with both, and requires the data files to be bitwise identical (and the
# potential-only printed output identical apart from paths, thread counts and "computed in N s"
# timings). The working tree is also run with 4 threads, which must not change the results either.
# Code paths not exercised by these runs need a direct old-vs-new comparison.
#
# Usage: tools/check_bitwise_refactor.sh [base-ref]      (base-ref default: HEAD)
# Output: the differences found, or one "OK" line. Exit 0 identical, 1 different or error, 2 usage.
# Environment: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (nlohmann_json sources for an offline
#   configure; default <repo>/build/_deps/nlohmann_json-src, cloned from GitHub if that is
#   missing), TMPDIR (scratch builds), CMAKE_BUILD_PARALLEL_LEVEL (build jobs, read by cmake);
#   the builds: tools/scratch_build.sh.
# Run time: about 50 s (two Release builds of spodvr and six runs; 48 s measured on the
#   development machine with CMAKE_BUILD_PARALLEL_LEVEL=4, 2026-11-12).
# Docs: docs/CONTRIBUTING.md, section "Testing".
set -euo pipefail

if [[ $# -gt 1 ]]; then
  echo "usage: $0 [base-ref]" >&2
  exit 2
fi
base_ref="${1:-HEAD}"
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
source "$repo/tools/scratch_build.sh"
require_commit "$base_ref"

echo "Building base ($base_ref) and working tree ..."
add_base_worktree "$base_ref"
build_spodvr "$scratch/base-src" "$scratch/base-build"
build_spodvr "$repo" "$scratch/new-build"

# Potential-only variant of the reference input.
python3 -I -B -c '
import json, sys
config = json.load(open(sys.argv[1]))
config["Potential"] = True
json.dump(config, open(sys.argv[2], "w"), indent=2)
' "$repo/tests/reference/hardwall_small/input.json" "$scratch/input.json"

write_short_propagation_input "$scratch/propagation.json"

# Printed output without what legitimately differs between runs: the output path, the thread
# count and loop schedule, and wall-clock timings ("6D potential computed in 0.42 s").
run() {  # <build dir> <output name> <threads>
  "$1/spodvr" "$scratch/input.json" --no-color --threads "$3" --output-dir "$scratch/$2" \
    | grep -vE 'Data files written to| computed in [0-9]+(\.[0-9]+)? s$' \
    | sed -E 's/[0-9]* OpenMP thread\(s\)(, loop schedule [^ ]+ \([^)]*\))?//' > "$scratch/$2.txt"
}
run "$scratch/base-build" base 1
run "$scratch/new-build" new 1
run "$scratch/new-build" new-threads 4

propagate() {  # <build dir> <output name> <threads>
  "$1/spodvr" "$scratch/propagation.json" --no-color --threads "$3" --output-dir "$scratch/$2" > /dev/null
}
echo "Running the short propagation ..."
propagate "$scratch/base-build" prop-base 4
propagate "$scratch/new-build" prop-new 1
propagate "$scratch/new-build" prop-new-threads 4

status=0
for file in asymptotic_potential.dat gas_phase_states.dat; do
  for other in new new-threads; do
    if ! cmp -s "$scratch/base/$file" "$scratch/$other/$file"; then
      echo "DIFFERENT: $file ($other vs base)"; status=1
    fi
  done
done
for other in prop-new prop-new-threads; do
  if ! cmp -s "$scratch/prop-base/averages.dat" "$scratch/$other/averages.dat"; then
    echo "DIFFERENT: averages.dat of the propagation ($other vs base)"; status=1
  fi
done
for other in new new-threads; do
  if ! diff -q "$scratch/base.txt" "$scratch/$other.txt" >/dev/null; then
    echo "DIFFERENT: printed output ($other vs base):"; diff "$scratch/base.txt" "$scratch/$other.txt" | head -20; status=1
  fi
done
[[ $status -eq 0 ]] && echo "OK: potential-only and propagation results are bitwise identical to $base_ref (1 and 4 threads)"
exit $status
