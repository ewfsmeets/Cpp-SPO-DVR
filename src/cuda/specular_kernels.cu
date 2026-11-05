// Kernel of the specular grid step after the kinetic part. CPU twin: propagation/specular_operator.cpp
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// One thread per specular point (ir, iz), index ir + nr * iz. The three CPU passes (gas phase,
// Z absorber, r absorber) are applied to each point in the same order, so the result is the same.
__global__ void specular_gas_phase_and_absorbers_kernel(SpecularStepArgs a, Complex* spec) {
  const std::int64_t i = thread_index();
  if (i >= std::int64_t(a.nr) * a.nz) return;
  const int ir = static_cast<int>(i % a.nr);
  const int iz = static_cast<int>(i / a.nr);
  Complex v = spec[i];
  if (iz >= a.nz_potential) v = a.gas_phase_step[ir] * v;
  if (iz >= a.z_absorber_first && iz <= a.z_absorber_last) v = a.z_optical[iz] * v;
  if (a.nr > 1 && ir >= a.r_absorber_first && ir <= a.r_absorber_last) v = a.r_optical[ir] * v;
  spec[i] = v;
}

}  // namespace

// Fortran: GasPhase.f90 PerformEXPGasPhaseOnSpecGrid, Absorption.f90 PerformEXPOpticalOpSpec
void apply_specular_gas_phase_and_absorbers(const SpecularStepArgs& args, Complex* specular) {
  const std::int64_t n = std::int64_t(args.nr) * args.nz;
  specular_gas_phase_and_absorbers_kernel<<<blocks_for(n), threads_per_block>>>(args, specular);
  check_kernel_launch("apply_specular_gas_phase_and_absorbers");
}

}  // namespace spodvr::cuda
