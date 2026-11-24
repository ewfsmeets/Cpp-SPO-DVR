#pragma once

// Shared set-up of the GPU tests: skip without a GPU (fail on a broken CUDA set-up), random data,
// relative differences.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <span>

#include <cuda_runtime_api.h>

#include "comparison.hpp"

namespace spodvr::testing {

// What a GPU test does on this machine.
enum class GpuCheck { run, skip, fail };

// The decision from what the CUDA runtime reports: driver_version from cudaDriverGetVersion (0: no
// driver installed), count_status and device_count from cudaGetDeviceCount. A machine without a GPU
// (no driver, no device) skips the test; any other error with a driver present fails it (e.g.
// cudaErrorInsufficientDriver: a driver older than the toolkit), so that a broken set-up does not
// pass as skipped tests.
inline GpuCheck gpu_check(int driver_version, cudaError_t count_status, int device_count) {
  if (driver_version == 0 || count_status == cudaErrorNoDevice) return GpuCheck::skip;
  if (count_status != cudaSuccess) return GpuCheck::fail;
  return device_count == 0 ? GpuCheck::skip : GpuCheck::run;
}

// Skips the current test when the machine has no CUDA-capable GPU; fails it when the CUDA runtime
// reports another error (gpu_check).
#define SPODVR_REQUIRE_GPU()                                                                                      \
  do {                                                                                                            \
    int spodvr_driver_version = 0;                                                                                \
    int spodvr_gpu_count = 0;                                                                                     \
    (void)cudaDriverGetVersion(&spodvr_driver_version);                                                           \
    const cudaError_t spodvr_count_status = cudaGetDeviceCount(&spodvr_gpu_count);                                \
    const ::spodvr::testing::GpuCheck spodvr_gpu_check =                                                          \
        ::spodvr::testing::gpu_check(spodvr_driver_version, spodvr_count_status, spodvr_gpu_count);               \
    if (spodvr_gpu_check == ::spodvr::testing::GpuCheck::skip)                                                    \
      GTEST_SKIP() << "no CUDA-capable GPU (driver version " << spodvr_driver_version << ", " << spodvr_gpu_count \
                   << " devices: " << cudaGetErrorString(spodvr_count_status) << ")";                             \
    if (spodvr_gpu_check == ::spodvr::testing::GpuCheck::fail)                                                    \
      GTEST_FAIL() << "CUDA driver " << spodvr_driver_version << " present, but cudaGetDeviceCount "              \
                   << "failed: " << cudaGetErrorString(spodvr_count_status);                                      \
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
