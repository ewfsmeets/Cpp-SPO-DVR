#pragma once

// Launch configuration shared by the kernels: one thread per work item, 256 threads per block.

#include <cstdint>

namespace spodvr::cuda {

inline constexpr int threads_per_block = 256;

// Number of blocks for n work items.
inline unsigned int blocks_for(std::int64_t n) {
  return static_cast<unsigned int>((n + threads_per_block - 1) / threads_per_block);
}

// The global index of this thread (64-bit: the wave function can have more than 2^31 elements).
__device__ inline std::int64_t thread_index() {
  return static_cast<std::int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
}

}  // namespace spodvr::cuda
