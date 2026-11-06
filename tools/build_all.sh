#!/usr/bin/env bash
# Configures (when needed) and builds every build variant of the project, see the "Building"
# section of README.md. Usage, from the repository root:
#   tools/build_all.sh                 all variants (the CUDA ones only when nvcc is found)
#   tools/build_all.sh build build-cuda   only these
# build-fast is built without profile-guided optimisation (the training run takes long, and the
# profile must be retrained after every code change); for PGO see docs/performance.md.
set -euo pipefail

cd "$(dirname "$0")/.."

release="-G Ninja -DCMAKE_BUILD_TYPE=Release"
declare -A options=(
  [build]=""
  [build-fftw]="-DSPODVR_FFT_LIBRARY=FFTW"
  [build-nomkl]="-DSPODVR_USE_MKL=OFF"
  [build-poison]="-DSPODVR_POISON_UNINITIALISED=ON"
  [build-fast]="-DSPODVR_FAST_MATH=ON -DSPODVR_ARCH=native -DSPODVR_PGO=OFF"
  [build-cuda]="-DSPODVR_USE_CUDA=ON"
  [build-cuda-fast]="-DSPODVR_USE_CUDA=ON -DSPODVR_FAST_MATH=ON -DSPODVR_ARCH=native -DSPODVR_PGO=OFF"
)
order=(build build-fftw build-nomkl build-poison build-fast build-cuda build-cuda-fast)

variants=("$@")
if [ ${#variants[@]} -eq 0 ]; then variants=("${order[@]}"); fi

# nvcc on PATH, or where NVIDIA's packages put it (CMake looks there too).
have_nvcc=0
if command -v nvcc > /dev/null || [ -x /usr/local/cuda/bin/nvcc ]; then have_nvcc=1; fi

# Configures and builds one variant with the given extra CMake options.
build_variant() {
  local dir=$1
  shift
  # shellcheck disable=SC2086  # the option strings are split on purpose
  cmake -S . -B "$dir" $release "$@" > /dev/null
  cmake --build "$dir"
}

for dir in "${variants[@]}"; do
  if [ -z "${options[$dir]+set}" ]; then
    echo "unknown build variant: $dir (known: ${order[*]})" >&2
    exit 2
  fi
  if [[ $dir == build-cuda* ]] && [ $have_nvcc -eq 0 ]; then
    echo "--- $dir: skipped (nvcc not found)"
    continue
  fi
  echo "--- $dir"
  # shellcheck disable=SC2086
  build_variant "$dir" ${options[$dir]}
done
echo "All requested variants built."
