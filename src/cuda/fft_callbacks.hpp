#pragma once

// cuFFT load callback of the GPU backend: exp(i A) applied while the inverse phi FFT loads the DVR
// values, instead of a separate multiply_by_exp_potential pass over the DVR scratch (one read and
// one write of it saved per slab). The callback (fft_callbacks.cu) is compiled to LTO-IR and
// embedded in the program; cuFFT links it into its FFT kernel at plan creation (nvJitLink).
// Plain C++: included by the callback and by the host code.

#include <cstdint>
#include <span>

#include "cuda/device_complex.hpp"

namespace spodvr::cuda {

// Everything the callback reads; lives in device memory (one struct), izr_first is updated before
// every slab in stream order. Same indexing as multiply_by_exp_potential_kernel
// (potential_kernels.cu): scratch element o = col + n_columns * (itheta + n_theta * iphi), column
// col = (izr - izr_first) + n_slab * ixy.
struct ExpPotentialCallbackInfo {
  const Complex* exp_potential;   // exp(i A), (izr, itp, irep) with izr fastest
  const int* representative_of;   // per ixy
  const int* operation_of;        // per ixy
  const int* phi_index;           // per symmetry operation: n_phi values
  std::int64_t n_zr_potential;
  std::int64_t izr_first;         // the current slab
  int n_slab;
  int n_columns;                  // n_slab * n_xy
  int n_theta;
  int n_phi;
};

// Name of the callback function in the LTO-IR (unmangled; cuFFT allows no namespace or extern "C").
inline constexpr char exp_potential_callback_name[] = "spodvr_exp_potential_load";

// The embedded LTO-IR fatbin of fft_callbacks.cu (generated at build time, fft_callbacks.cpp).
[[nodiscard]] std::span<const unsigned char> exp_potential_callback_fatbin();

}  // namespace spodvr::cuda
