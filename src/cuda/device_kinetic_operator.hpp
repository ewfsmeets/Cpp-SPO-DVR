#pragma once

#include <memory>
#include <vector>

#include "cuda/cuda_fft.hpp"
#include "cuda/device_buffer.hpp"
#include "cuda/device_wave_packet.hpp"
#include "cuda/kernels.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/specular_operator.hpp"

namespace spodvr::cuda {

// exp(-i T dt) on the scattering grid in GPU memory: T_Z + T_r by one batched 2D (Z, r) cuFFT
// over all (irot, ixy) blocks (contiguous in the device layout, so no gather as on the CPU),
// and T_xy (with the momentum shift) by one strided 2D (x, y) cuFFT over all (izr, irot).
// The FFT scale factors of both directions are applied in the multiplication kernel between the
// transforms (no extra pass over psi; see FftScales).
// CPU twin: ScatteringKineticOperator. Fortran: KineticOp.f90, SplitOperator.f90 PerformEXPKineticOp
class DeviceKineticOperator {
 public:
  DeviceKineticOperator(const ScatteringShape& shape, const PropagationTables& tables, std::size_t max_fft_work);

  // exp(-i (T_Z + T_r) dt). Fortran: KineticOpAlongZrOnScatGrid
  void apply_zr(DeviceWavePacket& psi, TimeStep step) const;
  // exp(-i T_xy dt) at every (Z, r) point. Fortran: PerformEXPKineticOp
  void apply_xy(DeviceWavePacket& psi, TimeStep step) const;

  // The cuFFT plans (to share one work area).
  [[nodiscard]] std::vector<DeviceFftPlan*> fft_plans();

 private:
  ScatteringShape shape_;
  DeviceBuffer<Complex> zr_full_, zr_half_;  // exp(-i (T_Z + T_r) dt), (ir, iz)
  DeviceBuffer<Complex> xy_full_, xy_half_;  // exp(-i T_xy dt), per ixy
  FftScales zr_scales_, xy_scales_;
  bool shift_;                               // momentum shift needed (parallel momentum != 0)
  DeviceBuffer<Complex> shift_forward_, shift_backward_;
  std::unique_ptr<DeviceFftPlan> zr_to_momentum_, zr_to_coordinate_;
  std::unique_ptr<DeviceFftPlan> xy_to_momentum_, xy_to_coordinate_;  // none if nx = ny = 1
};

// The specular (incoming channel) grid in GPU memory: 2D (r, Z) cuFFT kinetic step, then (except
// on the last half step) the gas-phase propagator and the absorbers.
// CPU twin: SpecularOperator. Fortran: SplitOperator.f90 SpecularChannelOps
class DeviceSpecularOperator {
 public:
  DeviceSpecularOperator(int nr, int nz_specular, int nz_potential, int initial_channel, const IndexRange& z_absorber,
                         const IndexRange& r_absorber, const PropagationTables& tables, std::size_t max_fft_work);

  void apply(DeviceWavePacket& psi, TimeStep step, bool kinetic_only) const;
  [[nodiscard]] std::vector<DeviceFftPlan*> fft_plans();

 private:
  int nr_;
  int nz_;
  int nz_potential_;
  IndexRange z_absorber_;
  IndexRange r_absorber_;
  DeviceBuffer<Complex> kinetic_full_, kinetic_half_;
  FftScales scales_;
  DeviceBuffer<Complex> gas_phase_step_;                // per ir
  DeviceBuffer<double> z_optical_;                      // per iz
  DeviceBuffer<double> r_optical_;                      // per ir
  std::unique_ptr<DeviceFftPlan> to_momentum_, to_coordinate_;
};

}  // namespace spodvr::cuda
