# shellcheck shell=bash
# The Release build of spodvr that tools/check_bitwise_refactor.sh and tools/compare_propagation.sh
# make of each tree they compare. Sourced, never run, after the caller has set repo (its
# repository root): source "$repo/tools/scratch_build.sh"
# Environment: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON (nlohmann_json sources for an offline
#   configure; default <repo>/build/_deps/nlohmann_json-src, cloned from GitHub if that is
#   missing; a value that is not a directory is refused), CMAKE_BUILD_PARALLEL_LEVEL (build jobs).
# Run time: sourcing, instant (the builds: see the callers). Docs: README.md, section "Testing".

json_src="${FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON:-$repo/build/_deps/nlohmann_json-src}"
if [[ -n "${FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON:-}" && ! -d "$json_src" ]]; then
  echo "error: FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON='$json_src' is not a directory; point it at the nlohmann_json sources or unset it" >&2
  exit 1
fi

# build_spodvr SOURCE_DIR BUILD_DIR: configures SOURCE_DIR (Release, no tests, no tools) into
# BUILD_DIR and builds spodvr.
build_spodvr() {
  local reuse=()
  [[ -d "$json_src" ]] && reuse=(-DFETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON="$json_src")
  cmake -S "$1" -B "$2" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_BUILD_TESTS=OFF -DSPODVR_BUILD_TOOLS=OFF \
    "${reuse[@]}" >/dev/null
  cmake --build "$2" --target spodvr >/dev/null
}
