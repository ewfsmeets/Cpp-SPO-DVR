// cuFFT LTO load callback: exp(i A) times the DVR value, applied while the inverse phi FFT loads it.
// CPU twin: PotentialOperator::apply_interaction (the multiplication between the phi transforms).
// Fortran: Potential.f90 PerformEXPInteractionOp
//
// Not part of the spodvr_cuda library: CMake compiles this file to LTO-IR (a fatbin, embedded by
// fft_callbacks.cpp), and cuFFT links it into the FFT kernel when the plan is made. cuFFT finds the
// function by its unmangled name, so it is at global scope and not extern "C" (both required).
// cuFFT chooses the link options (nvJitLink contracts to FMA by default), so the product is written
// with explicitly rounded operations: the same product as Complex operator* compiled with
// -fmad=false (multiply_by_exp_potential_kernel). The FFT around it is cuFFT's: for a plan with a
// callback it may use a differently rounding kernel (see DevicePotentialOperator).
#include <cufftXt.h>

#include "cuda/fft_callbacks.hpp"

__device__ cufftDoubleComplex spodvr_exp_potential_load(void* data_in, unsigned long long offset, void* caller_info,
                                                        void* /*shared*/) {
  const auto& a = *static_cast<const spodvr::cuda::ExpPotentialCallbackInfo*>(caller_info);
  // The scratch of one slab has fewer than 2^31 values (DevicePotentialOperator checks it).
  const int o = static_cast<int>(offset);
  const int col = o % a.n_columns;
  const int tp = o / a.n_columns;  // itheta + n_theta * iphi
  const int itheta = tp % a.n_theta;
  const int iphi = tp / a.n_theta;
  const std::int64_t izr = a.izr_first + col % a.n_slab;
  const int ixy = col / a.n_slab;
  // Surface symmetry: the (x, y) point takes the value of its representative at the mapped phi.
  const std::int64_t itp = a.phi_index[a.n_phi * a.operation_of[ixy] + iphi] + std::int64_t(a.n_phi) * itheta;
  const std::int64_t n_tp = std::int64_t(a.n_theta) * a.n_phi;
  const spodvr::cuda::Complex f = a.exp_potential[izr + a.n_zr_potential * (n_tp * a.representative_of[ixy] + itp)];
  const cufftDoubleComplex v = static_cast<const cufftDoubleComplex*>(data_in)[offset];
  // f * v as Complex operator* (device_complex.hpp): (f.re v.re - f.im v.im, f.re v.im + f.im v.re)
  return {__dsub_rn(__dmul_rn(f.re, v.x), __dmul_rn(f.im, v.y)), __dadd_rn(__dmul_rn(f.re, v.y), __dmul_rn(f.im, v.x))};
}
