#pragma once

// Launch configuration shared by the kernels: one thread per work item, 256 threads per block.

#include <cstdint>

namespace spodvr::cuda {

inline constexpr int threads_per_block = 256;

// Number of blocks for n work items.
inline unsigned int blocks_for(std::int64_t n) {
  return static_cast<unsigned int>((n + threads_per_block - 1) / threads_per_block);
}

// A 2D launch for an array of n_inner x n_outer values (inner index fastest): x covers the inner
// index, y the outer one. Used instead of a 1D launch with i % n_inner and i / n_inner: integer
// division is slow on the GPU (64-bit division is emulated in software), and these kernels are
// otherwise pure memory traffic. The outer index runs in a grid-stride loop (gridDim.y <= 65535).
inline dim3 grid_2d(std::int64_t n_inner, std::int64_t n_outer) {
  const std::int64_t y = n_outer < 65535 ? n_outer : 65535;
  return dim3(blocks_for(n_inner), static_cast<unsigned int>(y));
}

// The global index of this thread (64-bit: the wave function can have more than 2^31 elements).
__device__ inline std::int64_t thread_index() {
  return static_cast<std::int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
}

}  // namespace spodvr::cuda
