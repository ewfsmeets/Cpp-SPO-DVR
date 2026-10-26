#!/usr/bin/env bash
# clang-tidy on the given files with the build's compile database; called by the pre-commit hook
# (which pins and installs clang-tidy). Usage: tools/clang_tidy.sh FILE... (BUILD_DIR, default build)
# The pinned clang-tidy has no OpenMP header, and GCC's omp.h sits next to GCC's own intrinsics
# headers (which clang must not see), so only omp.h is exposed through a separate directory.
set -euo pipefail

build_dir=${BUILD_DIR:-build}
if [[ ! -f $build_dir/compile_commands.json ]]; then
  echo "clang_tidy.sh: $build_dir/compile_commands.json not found; configure with CMake first" >&2
  exit 1
fi

include_dir=$PWD/$build_dir/clang-tidy-include
mkdir -p "$include_dir"
ln -sf "$(${CXX:-g++} -print-file-name=include/omp.h)" "$include_dir/omp.h"

exec clang-tidy -p "$build_dir" --quiet --extra-arg=-isystem"$include_dir" "$@"
