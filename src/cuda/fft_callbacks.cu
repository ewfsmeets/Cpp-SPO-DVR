// cuFFT LTO load callback: exp(i A) times the DVR value, applied while the inverse phi FFT loads it.
// CPU twin: PotentialOperator::apply_interaction (the multiplication between the phi transforms).
// Fortran: Potential.f90 PerformEXPInteractionOp
//
// Not part of the spodvr_cuda library: CMake compiles this file to LTO-IR (a fatbin, embedded by
// fft_callbacks.cpp), and cuFFT links it into the FFT kernel when the plan is made. cuFFT finds the
// function by its name, so it is at global scope and not extern "C" (both required).
// cuFFT chooses the link options (nvJitLink contracts to FMA by default), so the product is written
// with explicitly rounded operations: the same product as Complex operator* compiled with
// -fmad=false (multiply_by_exp_potential_kernel). The FFT around it is cuFFT's: for a plan with a
// callback it may use a differently rounding kernel (see DevicePotentialOperator).
#include <cufftXt.h>

#include "cuda/fft_callbacks.hpp"

__device__ cufftDoubleComplex spodvr_exp_potential_load(void* data_in, unsigned long long offset, void* caller_info,
                                                        void* /*shared*/) {
  const auto& info = *static_cast<const spodvr::cuda::ExpPotentialCallbackInfo*>(caller_info);
  // The DVR scratch of the slab, (col, theta, phi): offset = col + n_columns * (itheta + n_theta * iphi).
  // 32-bit index arithmetic (cheaper division; DevicePotentialOperator keeps a slab below 2^31 values).
  const int value = static_cast<int>(offset);
  const int col = value % info.n_columns;
  const int theta_phi = value / info.n_columns;
  const int itheta = theta_phi % info.exp_potential.n_theta;
  const int iphi = theta_phi / info.exp_potential.n_theta;
  // Column col = (izr - izr_first) + n_slab * ixy.
  const std::int64_t izr = info.izr_first + col % info.n_slab;
  const int ixy = col / info.n_slab;

  const spodvr::cuda::Complex f = info.exp_potential.at(izr, ixy, itheta, iphi);
  const cufftDoubleComplex v = static_cast<const cufftDoubleComplex*>(data_in)[offset];
  // f * v as Complex operator* (device_complex.hpp): (f.re v.re - f.im v.im, f.re v.im + f.im v.re)
  return {__dsub_rn(__dmul_rn(f.re, v.x), __dmul_rn(f.im, v.y)), __dadd_rn(__dmul_rn(f.re, v.y), __dmul_rn(f.im, v.x))};
}
