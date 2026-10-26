#!/usr/bin/env bash
# Checks that a refactor did not change any numerical result: builds a base commit (default HEAD)
# and the working tree, runs the potential-only hard wall case with both, and requires the data
# files to be bitwise identical (and the printed output identical apart from paths). The working
# tree is also run with several threads, which must not change the results either.
#
# Only code exercised by the potential-only run is covered; check other code paths (e.g. the
# quadrature with mj > 0) with a direct old-vs-new comparison.
#
# Usage: tools/check_bitwise_refactor.sh [base-ref]
set -euo pipefail

base_ref="${1:-HEAD}"
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
scratch="$(mktemp -d)"
trap 'git -C "$repo" worktree remove --force "$scratch/base-src" >/dev/null 2>&1 || true; rm -rf "$scratch"' EXIT

configure_and_build() {  # <source dir> <build dir>
  local json_src="$repo/build/_deps/nlohmann_json-src"
  local reuse=()
  [[ -d "$json_src" ]] && reuse=(-DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON="$json_src")
  cmake -S "$1" -B "$2" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_BUILD_TESTS=OFF "${reuse[@]}" >/dev/null
  cmake --build "$2" --target spodvr >/dev/null
}

echo "Building base ($base_ref) and working tree ..."
git -C "$repo" worktree add --detach "$scratch/base-src" "$base_ref" >/dev/null 2>&1
configure_and_build "$scratch/base-src" "$scratch/base-build"
configure_and_build "$repo" "$scratch/new-build"

# Potential-only variant of the reference input.
python3 -c '
import json, sys
config = json.load(open(sys.argv[1]))
config["Potential"] = True
json.dump(config, open(sys.argv[2], "w"), indent=2)
' "$repo/tests/reference/hardwall_small/input.json" "$scratch/input.json"

run() {  # <build dir> <output name> <threads>
  "$1/spodvr" "$scratch/input.json" --no-color --threads "$3" --output-dir "$scratch/$2" \
    | grep -v "Data files written to" | sed 's/[0-9]* OpenMP thread(s)//' > "$scratch/$2.txt"
}
run "$scratch/base-build" base 1
run "$scratch/new-build" new 1
run "$scratch/new-build" new-threads 4

status=0
for file in asymptotic_potential.dat gas_phase_states.dat; do
  for other in new new-threads; do
    if ! cmp -s "$scratch/base/$file" "$scratch/$other/$file"; then
      echo "DIFFERENT: $file ($other vs base)"; status=1
    fi
  done
done
for other in new new-threads; do
  if ! diff -q "$scratch/base.txt" "$scratch/$other.txt" >/dev/null; then
    echo "DIFFERENT: printed output ($other vs base):"; diff "$scratch/base.txt" "$scratch/$other.txt" | head -20; status=1
  fi
done
[[ $status -eq 0 ]] && echo "OK: results are bitwise identical to $base_ref (1 and 4 threads)"
exit $status
