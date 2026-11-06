#pragma once

// cuFFT load callback of the GPU backend: exp(i A) applied while the inverse phi FFT loads the DVR
// values, instead of a separate multiply_by_exp_potential pass over the DVR scratch (one read and
// one write of it saved per slab). The callback (fft_callbacks.cu) is compiled to LTO-IR and
// embedded in the program; cuFFT links it into its FFT kernel at plan creation (nvJitLink).
// Plain C++: included by the callback and by the host code.

#include <cstdint>
#include <span>

#include "cuda/device_complex.hpp"
#include "cuda/kernels.hpp"

namespace spodvr::cuda {

// Everything the callback reads; lives in device memory (one struct), izr_first is updated before
// every slab in stream order. The DVR scratch of a slab as in potential_kernels.cu: element
// o = col + n_columns * (itheta + n_theta * iphi), column col = (izr - izr_first) + n_slab * ixy.
struct ExpPotentialCallbackInfo {
  ExpPotentialTable exp_potential;
  std::int64_t izr_first;  // the current slab
  int n_slab;
  int n_columns;           // n_slab * n_xy
};

// Name of the callback function in the LTO-IR (unmangled; cuFFT allows no namespace or extern "C").
inline constexpr char exp_potential_callback_name[] = "spodvr_exp_potential_load";

// The embedded LTO-IR fatbin of fft_callbacks.cu (generated at build time, fft_callbacks.cpp).
[[nodiscard]] std::span<const unsigned char> exp_potential_callback_fatbin();

}  // namespace spodvr::cuda
