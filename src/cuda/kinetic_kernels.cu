// Kernels of the kinetic stages: multiplication by exp(-i T dt) in momentum space (between the
// cuFFT transforms) and by the momentum-shift phases. CPU twin: propagation/kinetic_operators.cpp
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// One thread per element; consecutive threads take consecutive elements (coalesced).
__global__ void multiply_by_periodic_factors_kernel(Complex* data, const Complex* factors, FftScales scales,
                                                    std::int64_t period, std::int64_t n) {
  const std::int64_t i = thread_index();
  if (i >= n) return;
  data[i] = factors[i % period] * (scales.backward * (scales.forward * data[i]));
}

__global__ void multiply_by_blocked_factors_kernel(Complex* data, const Complex* factors, FftScales scales,
                                                   std::int64_t block, std::int64_t n) {
  const std::int64_t i = thread_index();
  if (i >= n) return;
  data[i] = factors[i / block] * (scales.backward * (scales.forward * data[i]));
}

}  // namespace

// Fortran: KineticOp.f90 MultByEXPKinEnergy (Z, r and specular grid)
void multiply_by_periodic_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t period,
                                  std::int64_t n) {
  multiply_by_periodic_factors_kernel<<<blocks_for(n), threads_per_block>>>(data, factors, scales, period, n);
  check_kernel_launch("multiply_by_periodic_factors");
}

// Fortran: SplitOperator.f90 PerformEXPKineticOp (MomentumShiftWF, MultByEXPKinEnergy)
void multiply_by_blocked_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t block,
                                 std::int64_t n) {
  multiply_by_blocked_factors_kernel<<<blocks_for(n), threads_per_block>>>(data, factors, scales, block, n);
  check_kernel_launch("multiply_by_blocked_factors");
}

}  // namespace spodvr::cuda
