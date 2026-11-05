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

}  // namespace

// Fortran: KineticOp.f90 MultByEXPKinEnergy (Z, r and specular grid)
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
