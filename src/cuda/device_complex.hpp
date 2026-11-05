#pragma once

// Complex arithmetic and the wave-function layout on the GPU. Included by the CUDA kernels (.cu)
// and by the host code that launches them (.cpp), so everything here is plain C++ that also
// compiles with nvcc.

#include <cstdint>

#ifdef __CUDACC__
#define SPODVR_HOST_DEVICE __host__ __device__
#else
#define SPODVR_HOST_DEVICE
#endif

namespace spodvr::cuda {

// A complex number in device memory. Same memory layout as std::complex<double> (two doubles),
// so host arrays of std::complex<double> are copied to the GPU as they are. (std::complex
// cannot be used in device code.)
struct alignas(16) Complex {
  double re;
  double im;
};

// The products as GCC computes std::complex products with -fcx-fortran-rules (no FMA: the
// kernels are compiled with -fmad=false), so a kernel gives bitwise the same result as the CPU
// code for the same operands in the same order.
SPODVR_HOST_DEVICE inline Complex operator*(Complex a, Complex b) {
  return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
}
// real * complex, as std::complex<double> does it: (r re, r im)
SPODVR_HOST_DEVICE inline Complex operator*(double r, Complex a) { return {r * a.re, r * a.im}; }
SPODVR_HOST_DEVICE inline Complex operator+(Complex a, Complex b) { return {a.re + b.re, a.im + b.im}; }
SPODVR_HOST_DEVICE inline Complex operator-(Complex a, Complex b) { return {a.re - b.re, a.im - b.im}; }

// Layout of the scattering wave function on the GPU: (izr, irot, ixy), izr = ir + nr*iz fastest.
// (The host layout is (irot, ixy, izr).) With izr fastest every operator stage becomes either one
// batched FFT with one batch dimension or a kernel whose neighbouring threads read neighbouring
// memory (coalesced):
//   - each (irot, ixy) is a contiguous (Z, r) block: the 2D (Z, r) FFT needs no gather;
//   - the (x, y) FFT of all (izr, irot) is one strided plan with batches at distance 1;
//   - one thread per (Z, r) point (or per (Z, r, x, y) column) reads consecutive memory.
struct DeviceLayout {
  std::int64_t n_zr;   // nr * nz
  std::int64_t n_rot;  // rotational states
  std::int64_t n_xy;   // nx * ny

  SPODVR_HOST_DEVICE std::int64_t operator()(std::int64_t izr, std::int64_t irot, std::int64_t ixy) const {
    return izr + n_zr * (irot + n_rot * ixy);
  }
  SPODVR_HOST_DEVICE std::int64_t size() const { return n_zr * n_rot * n_xy; }
};

}  // namespace spodvr::cuda
