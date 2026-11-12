#pragma once

// Shared set-up of the GPU tests: skip without a GPU, random data, relative differences.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <span>

#include <cuda_runtime_api.h>

#include "comparison.hpp"

namespace spodvr::testing {

// Skips the current test when the machine has no CUDA-capable GPU.
#define SPODVR_REQUIRE_GPU()                                                           \
  do {                                                                                 \
    int spodvr_gpu_count = 0;                                                          \
    if (cudaGetDeviceCount(&spodvr_gpu_count) != cudaSuccess || spodvr_gpu_count == 0) \
      GTEST_SKIP() << "no CUDA-capable GPU";                                           \
  } while (false)

// Reproducible pseudo-random values in [-1, 1) (a 64-bit LCG; no <random> dependence).
inline void fill_pseudo_random(std::span<std::complex<double>> v, std::uint64_t seed = 12345) {
  std::uint64_t x = seed;
  auto next = [&x] {
    x = x * 6364136223846793005ULL + 1442695040888963407ULL;
    return double(x >> 11) / double(1ULL << 53) * 2.0 - 1.0;
  };
  for (auto& c : v) {
    const double re = next();
    c = {re, next()};
  }
}

// max |a - b| / max |a|; NaN if a value of a or b is NaN.
inline double relative_difference(std::span<const std::complex<double>> a, std::span<const std::complex<double>> b) {
  double diff = 0.0;
  double scale = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    diff = nan_aware_max(diff, std::abs(a[i] - b[i]));
    scale = nan_aware_max(scale, std::abs(a[i]));
  }
  return scale > 0.0 ? diff / scale : diff;
}

}  // namespace spodvr::testing
