#!/usr/bin/env bash
# Checks that a change did not alter any result of a run: builds a base commit and the working tree
# (or takes existing builds of the base), runs every run mode with both and requires every output
# file to be bitwise identical, the printed output identical apart from paths, thread counts and
# timings, and run_info.json identical apart from the git commit and paths. The cases (all run from
# the repository root, so that the CRP and SCM data paths resolve):
#   pot        hard wall, potential only
#   prop       hard wall, 50 time steps
#   first      hard wall, 30 steps; restart: 20 more steps from its checkpoint (run mode 1);
#   analysis   onlyAnalysis of first's BK file (run mode 2)
#   also       prop with alsoAnalysis; both: prop with CheckpointFormat both
#   sym        cu111_dailight_short_sym (UseSymmetry), 50 steps
#   plot       cu111_dailight_short, 50 steps, with the plot data of the base: densities/, potential/ (two sites)
#              and frames/ (every 50 a.u.: tStart, 50, 100)
#   crp, scm   potential-only set-ups of crp_cu111_srp48_small and scm_cu111_srp48_small
# The builds run every case with 4 threads; the working tree also runs pot and prop with 1 thread (for
# prop only the files are compared: the printed memory estimate names the thread count). With --cuda
# the propagating cases (prop, first, restart, also, both, sym, plot) also run on the GPU with CUDA builds of
# both trees. Code paths not exercised by these runs need a direct old-vs-new comparison. Files the
# working tree writes and the base does not are reported, except run_info.json (top level, or in
# fine_analysis/ for the analysis-only run) while the base does not write it yet.
#
# Usage: tools/check_bitwise_refactor.sh [--cuda] [--base-build DIR] [base-ref]   (base-ref default: HEAD)
#   --base-build DIR  use DIR/build/spodvr (and DIR/build-cuda/spodvr with --cuda), builds of base-ref
#                     made with build_spodvr / build_spodvr_cuda (tools/scratch_build.sh), instead of
#                     building the base; that DIR matches base-ref is not checked
# Output: the differences found, or one "OK" line. Exit 0 identical, 1 different or error, 2 usage.
# Environment: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (nlohmann_json sources for an offline
#   configure; default <repo>/build/_deps/nlohmann_json-src, cloned from GitHub if that is
#   missing), TMPDIR (scratch builds), CMAKE_BUILD_PARALLEL_LEVEL (build jobs, read by cmake);
#   the builds: tools/scratch_build.sh.
# Run time: 90 s with --cuda --base-build on the development machine (8 cores, RTX 4070 Ti, 2026-11-19:
#   the working tree's CPU and CUDA builds of spodvr and 34 runs); without --base-build add the base's builds.
# Docs: docs/CONTRIBUTING.md, section "Testing".
set -euo pipefail

usage() {
  echo "usage: $0 [--cuda] [--base-build DIR] [base-ref]" >&2
  exit 2
}
cuda=0
base_build=""
refs=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cuda) cuda=1; shift ;;
    --base-build) [[ $# -ge 2 ]] || usage; base_build="$(cd "$2" && pwd)"; shift 2 ;;
    -*) usage ;;
    *) refs+=("$1"); shift ;;
  esac
done
[[ ${#refs[@]} -le 1 ]] || usage
base_ref="${refs[0]:-HEAD}"
repo="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
source "$repo/tools/scratch_build.sh"
require_commit "$base_ref"

if [[ -n "$base_build" ]]; then
  scratch="$(mktemp -d)"
  trap 'rm -rf "$scratch"' EXIT
  base_cpu="$base_build/build"
  base_gpu="$base_build/build-cuda"
  required=("$base_cpu/spodvr")
  if [[ $cuda -eq 1 ]]; then required+=("$base_gpu/spodvr"); fi
  for exe in "${required[@]}"; do
    [[ -x "$exe" ]] || { echo "error: $exe does not exist (--base-build)" >&2; exit 1; }
  done
  echo "Using the base builds in $base_build; building the working tree ..."
else
  echo "Building base ($base_ref) and working tree ..."
  add_base_worktree "$base_ref"
  base_cpu="$scratch/base-build"
  base_gpu="$scratch/base-build-cuda"
  build_spodvr "$scratch/base-src" "$base_cpu"
  if [[ $cuda -eq 1 ]]; then build_spodvr_cuda "$scratch/base-src" "$base_gpu"; fi
fi
new_cpu="$scratch/new-build"
new_gpu="$scratch/new-build-cuda"
build_spodvr "$repo" "$new_cpu"
if [[ $cuda -eq 1 ]]; then build_spodvr_cuda "$repo" "$new_gpu"; fi

# write_inputs LABEL: the input of every case for the runs named LABEL into $scratch/inputs/LABEL/
# (restart and analysis read the checkpoint and BK file of that label's own first run).
write_inputs() {
  mkdir -p "$scratch/inputs/$1"
  python3 -I -B - "$repo" "$scratch/inputs/$1" "$scratch/$1/first" <<'PY'
import json, sys
repo, out, first_dir = sys.argv[1:4]
def case(name):
    return json.load(open(f"{repo}/tests/reference/{name}/input.json"))
def steps(config, n, t_start=0.0):
    config.update({"tStart": t_start, "tEnd": t_start + n * config["dtStep"], "nStep": n,
                   "dtAve": config["dtAna"]})
    return config
inputs = {
    "pot": {**case("hardwall_small"), "Potential": True},
    "prop": steps(case("hardwall_short"), 50),
    "first": steps(case("hardwall_short"), 30),
    "restart": {**steps(case("hardwall_short"), 20, 60.0), "NewCalc": False, "RestartDir": first_dir},
    "analysis": {**steps(case("hardwall_short"), 30), "onlyAnalysis": True, "RestartDir": first_dir},
    "also": {**steps(case("hardwall_short"), 50), "alsoAnalysis": True},
    "both": {**steps(case("hardwall_short"), 50), "CheckpointFormat": "both"},
    "sym": steps(case("cu111_dailight_short_sym"), 50),
    "plot": {
        **steps(case("cu111_dailight_short"), 50),
        "writeDensities": True,
        "writePotential": True,
        "PlotSites": [{"name": "top", "x": 0.0, "y": 0.0}, {"name": "bridge", "x": 0.5, "y": 0.0}],
        "dtFrames": 50.0,
    },
    "crp": {**case("crp_cu111_srp48_small"), "Potential": True},
    "scm": {**case("scm_cu111_srp48_small"), "Potential": True},
}
for name, config in inputs.items():
    json.dump(config, open(f"{out}/{name}.json", "w"), indent=1)
PY
}

# run_case EXE LABEL CASE THREADS [cuda]: runs one case from the repository root into
# $scratch/LABEL/CASE (stdout in $scratch/LABEL/CASE.txt); a failing run stops the check.
run_case() {
  local exe="$1" label="$2" name="$3" threads="$4" device=()
  if [[ "${5:-}" == cuda ]]; then device=(--device cuda); fi
  mkdir -p "$scratch/$label"
  if ! (cd "$repo" && "$exe" "$scratch/inputs/$label/$name.json" --no-color --threads "$threads" "${device[@]}" \
          --output-dir "$scratch/$label/$name" > "$scratch/$label/$name.txt" 2> "$scratch/$label/$name.err"); then
    echo "error: $label/$name failed:" >&2
    cat "$scratch/$label/$name.err" >&2
    exit 1
  fi
}

cpu_cases=(pot prop first restart analysis also both sym plot crp scm)
gpu_cases=(prop first restart also both sym plot)
echo "Running ${#cpu_cases[@]} cases with each build ..."
for label in base new; do
  write_inputs "$label"
  exe="$base_cpu/spodvr"
  if [[ $label == new ]]; then exe="$new_cpu/spodvr"; fi
  for name in "${cpu_cases[@]}"; do run_case "$exe" "$label" "$name" 4; done
done
write_inputs new1
run_case "$new_cpu/spodvr" new1 pot 1
run_case "$new_cpu/spodvr" new1 prop 1
if [[ $cuda -eq 1 ]]; then
  echo "Running ${#gpu_cases[@]} cases on the GPU with each CUDA build ..."
  for label in base-gpu new-gpu; do
    write_inputs "$label"
    exe="$base_gpu/spodvr"
    if [[ $label == new-gpu ]]; then exe="$new_gpu/spodvr"; fi
    for name in "${gpu_cases[@]}"; do run_case "$exe" "$label" "$name" 4 cuda; done
  done
fi

# The printed output without what legitimately differs between runs: paths, the thread count and
# loop schedule, and wall-clock timings.
filtered() {
  grep -vE 'Data files written to|Plot data written to|Frames \(dtFrames\):|Norm budget \(NormBudget\):|Restarted from|Fine analysis of|Fine-analysis files written to| computed in [0-9]+(\.[0-9]+)? s$|Propagation time:|Fine analysis time:' "$1" \
    | sed -E 's/[0-9]* OpenMP thread\(s\)(, loop schedule [^ ]+ \([^)]*\))?//'
}

# run_info.json without the git commit and the paths, which differ between builds and runs.
run_info_filtered() {
  python3 -I -B -c '
import json, sys
info = json.load(open(sys.argv[1]))
info.get("program", {}).pop("git", None)
for key in ("input_file", "restart_dir"):
    info.pop(key, None)
print(json.dumps(info, indent=1, sort_keys=True))' "$1"
}

status=0
# compare_case BASE_LABEL NEW_LABEL CASE [files-only]: file sets, every file's bytes, run_info.json,
# and (unless files-only) the printed output.
compare_case() {
  local a="$scratch/$1/$3" b="$scratch/$2/$3" file only_a only_b
  only_a="$(comm -23 <(cd "$a" && find . -type f | sort) <(cd "$b" && find . -type f | sort))"
  only_b="$(comm -13 <(cd "$a" && find . -type f | sort) <(cd "$b" && find . -type f | sort) \
              | { grep -vxE '\./(fine_analysis/)?run_info\.json' || true; })"
  if [[ -n "$only_a$only_b" ]]; then
    echo "DIFFERENT: file set of $3 ($2 vs $1): only in $1: ${only_a:-none}; only in $2: ${only_b:-none}"
    status=1
  fi
  while IFS= read -r file; do
    [[ -f "$b/$file" ]] || continue
    if [[ "$file" == */run_info.json ]]; then
      if ! diff -q <(run_info_filtered "$a/$file") <(run_info_filtered "$b/$file") > /dev/null; then
        echo "DIFFERENT: $3/run_info.json ($2 vs $1):"
        diff <(run_info_filtered "$a/$file") <(run_info_filtered "$b/$file") | head -10 || true
        status=1
      fi
    elif ! cmp -s "$a/$file" "$b/$file"; then
      echo "DIFFERENT: $3/${file#./} ($2 vs $1)"
      status=1
    fi
  done < <(cd "$a" && find . -type f | sort)
  [[ "${4:-}" == files-only ]] && return 0
  if ! diff -q <(filtered "$scratch/$1/$3.txt") <(filtered "$scratch/$2/$3.txt") > /dev/null; then
    echo "DIFFERENT: printed output of $3 ($2 vs $1):"
    diff <(filtered "$scratch/$1/$3.txt") <(filtered "$scratch/$2/$3.txt") | head -20 || true
    status=1
  fi
}

for name in "${cpu_cases[@]}"; do compare_case base new "$name"; done
compare_case base new1 pot
compare_case base new1 prop files-only
if [[ $cuda -eq 1 ]]; then
  for name in "${gpu_cases[@]}"; do compare_case base-gpu new-gpu "$name"; done
fi
compared=$(( ${#cpu_cases[@]} + 2 + cuda * ${#gpu_cases[@]} ))
if [[ $status -eq 0 ]]; then
  echo "OK: $compared case comparisons bitwise identical to $base_ref (files, printed output, run_info.json)"
fi
exit $status
