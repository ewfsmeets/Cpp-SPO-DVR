# shellcheck shell=bash
# The helpers of the tools that build trees of this repository in scratch directories:
# tools/check_bitwise_refactor.sh and tools/compare_propagation.sh (a Release build of spodvr from a
# base commit and one from the working tree) and tools/crp_mutation_check/run.sh (offline sources).
# Sourced, never run, after the caller has set repo (its repository root):
#   source "$repo/tools/scratch_build.sh"
# Sourcing sets source_options for nlohmann_json (offline_sources).
# Environment: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON, FETCHCONTENT_SOURCE_DIR_GOOGLETEST (sources
#   for an offline configure; default <repo>/build/_deps/<name>-src, cloned from GitHub if that is
#   missing; a value that is not a directory is refused), CMAKE_BUILD_PARALLEL_LEVEL (build jobs),
#   TMPDIR (add_base_worktree).
# Run time: sourcing, instant (the builds: see the callers). Docs: docs/CONTRIBUTING.md, section "Testing".

# offline_sources NAME...: sets the array source_options to -DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir> for
# each FetchContent dependency NAME (NLOHMANN_JSON, GOOGLETEST) whose sources are at hand: the
# directory FETCHCONTENT_SOURCE_DIR_<NAME> names, else <repo>/build/_deps/<name>-src if it exists.
offline_sources() {
  local name variable dir
  source_options=()
  for name in "$@"; do
    variable="FETCHCONTENT_SOURCE_DIR_$name"
    dir="${!variable:-$repo/build/_deps/${name,,}-src}"
    if [[ -n "${!variable:-}" && ! -d "$dir" ]]; then
      echo "error: $variable='$dir' is not a directory; point it at the ${name,,} sources or unset it" >&2
      exit 1
    fi
    if [[ -d "$dir" ]]; then source_options+=(-D"$variable=$dir"); fi
  done
}

# require_commit REF: refuses (exit 1) unless REF names a commit of the repository.
require_commit() {
  if ! git -C "$repo" rev-parse --verify --quiet "$1^{commit}" >/dev/null; then
    echo "error: '$1' is not a commit in $repo; give a branch, tag or commit hash" >&2
    exit 1
  fi
}

# add_base_worktree REF: sets scratch to a new mktemp -d directory and checks REF out into
# $scratch/base-src (a detached git worktree); both are removed when the script exits.
add_base_worktree() {
  scratch="$(mktemp -d)"
  trap 'git -C "$repo" worktree remove --force "$scratch/base-src" >/dev/null 2>&1 || true; rm -rf "$scratch"' EXIT
  git -C "$repo" worktree add --detach "$scratch/base-src" "$1" >/dev/null 2>&1
}

# build_spodvr SOURCE_DIR BUILD_DIR: configures SOURCE_DIR (Release, no tests, no tools) into
# BUILD_DIR and builds spodvr.
build_spodvr() {
  cmake -S "$1" -B "$2" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_BUILD_TESTS=OFF -DSPODVR_BUILD_TOOLS=OFF \
    "${source_options[@]}" >/dev/null
  cmake --build "$2" --target spodvr >/dev/null
}

# write_short_propagation_input FILE: the input of tests/reference/hardwall_short cut to 50 time
# steps, with the averages after every analysis step.
write_short_propagation_input() {
  python3 -I -B -c '
import json, sys
config = json.load(open(sys.argv[1]))
config.update({"tEnd": 100.0, "nStep": 50, "dtAve": 10.0})
json.dump(config, open(sys.argv[2], "w"), indent=2)
' "$repo/tests/reference/hardwall_short/input.json" "$1"
}

offline_sources NLOHMANN_JSON
