// Kernels of the kinetic stages: multiplication by exp(-i T dt) in momentum space (between the
// cuFFT transforms) and by the momentum-shift phases. CPU twin: propagation/kinetic_operators.cpp
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// data = (period, n / period): x over the position in the period (consecutive threads,
// consecutive elements: coalesced), y over the periods (grid_2d).
__global__ void multiply_by_periodic_factors_kernel(Complex* data, const Complex* factors, FftScales scales,
                                                    std::int64_t period, std::int64_t n_periods) {
  const std::int64_t k = thread_index();
  if (k >= period) return;
  const Complex factor = factors[k];
  for (std::int64_t p = blockIdx.y; p < n_periods; p += gridDim.y) {
    Complex& v = data[k + period * p];
    v = factor * (scales.backward * (scales.forward * v));
  }
}

// data = (block, n / block): x over the position in the block, y over the blocks.
__global__ void multiply_by_blocked_factors_kernel(Complex* data, const Complex* factors, FftScales scales,
                                                   std::int64_t block, std::int64_t n_blocks) {
  const std::int64_t i = thread_index();
  if (i >= block) return;
  for (std::int64_t b = blockIdx.y; b < n_blocks; b += gridDim.y) {
    Complex& v = data[i + block * b];
    v = factors[b] * (scales.backward * (scales.forward * v));
  }
}

// x over izr (coalesced), y over ixy (grid_2d); each thread runs over all irot of its (izr, ixy),
// n_zr apart (coalesced across the threads), so the factors are loaded once and no index
// division is needed. One pass over psi with both kinetic factors, (Z, r) first, then (x, y), as
// the CPU applies its two stages.
__global__ void multiply_by_kinetic_factors_kernel(Complex* psi, DeviceLayout layout, const Complex* zr_factors,
                                                   FftScales zr_scales, const Complex* xy_factors,
                                                   FftScales xy_scales) {
  const std::int64_t izr = thread_index();
  if (izr >= layout.n_zr) return;
  const Complex zr = zr_factors[izr];
  for (std::int64_t ixy = blockIdx.y; ixy < layout.n_xy; ixy += gridDim.y) {
    const Complex xy = xy_factors[ixy];
    for (std::int64_t irot = 0; irot < layout.n_rot; ++irot) {
      Complex& v = psi[layout(izr, irot, ixy)];
      const Complex after_zr = zr * (zr_scales.backward * (zr_scales.forward * v));
      v = xy * (xy_scales.backward * (xy_scales.forward * after_zr));
    }
  }
}

}  // namespace

// Fortran: KineticOp.f90 MultByEXPKinEnergy (Z, r and x, y in one pass)
void multiply_by_kinetic_factors(Complex* psi, const DeviceLayout& layout, const Complex* zr_factors,
                                 FftScales zr_scales, const Complex* xy_factors, FftScales xy_scales) {
  multiply_by_kinetic_factors_kernel<<<grid_2d(layout.n_zr, layout.n_xy), threads_per_block>>>(
      psi, layout, zr_factors, zr_scales, xy_factors, xy_scales);
  check_kernel_launch("multiply_by_kinetic_factors");
}

// Fortran: KineticOp.f90 MultByEXPKinEnergy (specular grid)
void multiply_by_periodic_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t period,
                                  std::int64_t n) {
  multiply_by_periodic_factors_kernel<<<grid_2d(period, n / period), threads_per_block>>>(data, factors, scales,
                                                                                          period, n / period);
  check_kernel_launch("multiply_by_periodic_factors");
}

// Fortran: SplitOperator.f90 PerformEXPKineticOp (MomentumShiftWF, MultByEXPKinEnergy)
void multiply_by_blocked_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t block,
                                 std::int64_t n) {
  multiply_by_blocked_factors_kernel<<<grid_2d(block, n / block), threads_per_block>>>(data, factors, scales, block,
                                                                                       n / block);
  check_kernel_launch("multiply_by_blocked_factors");
}

}  // namespace spodvr::cuda
