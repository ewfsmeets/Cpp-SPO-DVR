// Kernel of the specular grid step after the kinetic part. CPU twin: propagation/specular_operator.cpp
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// One thread per specular point, point = ir + nr * iz. The three CPU passes (gas phase,
// Z absorber, r absorber) are applied to each point in the same order, so the result is the same.
__global__ void specular_gas_phase_and_absorbers_kernel(SpecularStepArgs args, Complex* spec) {
  const std::int64_t point = thread_index();
  if (point >= std::int64_t(args.nr) * args.nz) return;
  const int ir = static_cast<int>(point % args.nr);
  const int iz = static_cast<int>(point / args.nr);
  Complex v = spec[point];
  if (iz >= args.nz_potential) v = args.gas_phase_step[ir] * v;
  v = args.z_optical[iz] * v;
  if (args.nr > 1 && ir >= args.r_absorber_first && ir <= args.r_absorber_last) v = args.r_optical[ir] * v;
  spec[point] = v;
}

}  // namespace

// Fortran: GasPhase.f90 PerformEXPGasPhaseOnSpecGrid, Absorption.f90 PerformEXPOpticalOpSpec
void apply_specular_gas_phase_and_absorbers(const SpecularStepArgs& args, Complex* specular) {
  const std::int64_t n = std::int64_t(args.nr) * args.nz;
  specular_gas_phase_and_absorbers_kernel<<<blocks_for(n), threads_per_block>>>(args, specular);
  check_kernel_launch("apply_specular_gas_phase_and_absorbers");
}

}  // namespace spodvr::cuda
