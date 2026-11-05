#pragma once

#include <memory>
#include <span>
#include <vector>

#include "cuda/cuda_fft.hpp"
#include "cuda/device_buffer.hpp"
#include "cuda/device_wave_packet.hpp"
#include "cuda/kernels.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "symmetry/grid_symmetry.hpp"

namespace spodvr::cuda {

// exp(-i V dt) on the scattering grid in GPU memory, with absorption and the coupling to the
// specular grid (the same steps as the CPU PotentialOperator::apply_block):
//   beyond the interaction region: one elementwise kernel;
//   in the interaction region, slab by slab of (Z, r) points:
//     absorb_and_transform_to_dvr -> phi FFT -> multiply_by_exp_potential -> inverse phi FFT
//     -> transform_to_fbr_and_rotate -> project_onto_specular.
// The slabs keep the DVR scratch array (n_theta * n_phi values per (Z, r, x, y) point, about
// twice psi) to a fixed size. exp(i A) is computed on the host with the same cos/sin as the CPU
// code (bitwise the same table) and kept on the GPU, re-laid out with izr fastest.
// CPU twin: PotentialOperator. Fortran: Potential.f90 PerformEXPPotentialOp
class DevicePotentialOperator {
 public:
  // Default size limit of the DVR scratch array of one slab.
  static constexpr std::size_t default_scratch_limit = std::size_t(1) << 30;  // 1 GiB

  // premultiplied_potential: A = -V dt on (itp, irep, izr) for the interaction region
  // (PotentialGrid after premultiply_with_time_step).
  DevicePotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                          const PropagationTables& tables, const AngularTransform& transform,
                          std::span<const double> premultiplied_potential, const GridSymmetry& symmetry,
                          std::size_t max_fft_work, std::size_t scratch_limit = default_scratch_limit);

  void apply(DeviceWavePacket& psi) const;
  [[nodiscard]] std::vector<DeviceFftPlan*> fft_plans();
  // (Z, r) points per slab.
  [[nodiscard]] std::int64_t slab_points() const { return args_.n_slab; }
  // The kernels' arguments and the DVR scratch (for the timing tool, which runs the kernels of
  // apply() one by one).
  [[nodiscard]] const PotentialArgs& kernel_args() const { return args_; }
  [[nodiscard]] Complex* scratch() const { return dvr_.data(); }

  // GPU memory of the tables and the scratch array (bytes).
  [[nodiscard]] static std::size_t memory_bytes(const ScatteringShape& shape, int nz_potential, int n_theta, int n_phi,
                                                int n_representatives,
                                                std::size_t scratch_limit = default_scratch_limit);
  // The (Z, r) points per slab: the largest divisor of the interaction-region points whose DVR
  // scratch fits in scratch_limit (at least 1), so that all slabs have the same shape.
  [[nodiscard]] static std::int64_t slab_size(const ScatteringShape& shape, int nz_potential, int n_theta, int n_phi,
                                              std::size_t scratch_limit);

 private:
  PotentialArgs args_{};  // the kernels' view of the tables below (izr_first set per slab)
  DeviceBuffer<double> z_optical_;
  DeviceBuffer<Complex> before_interaction_, gas_phase_step_, rotation_half_;
  DeviceBuffer<Complex> specular_to_scattering_, scattering_to_specular_;
  DeviceBuffer<int> block_first_, block_size_, block_phi_slot_, block_matrix_offset_;
  DeviceBuffer<double> matrix_;
  DeviceBuffer<Complex> exp_potential_;
  DeviceBuffer<int> representative_of_, operation_of_, phi_index_;
  mutable DeviceBuffer<Complex> dvr_;  // scratch of one slab, (col, theta, phi); not part of the state
  std::unique_ptr<DeviceFftPlan> to_dvr_fft_, to_fbr_fft_;  // phi transforms (none if nphi = 1)
};

}  // namespace spodvr::cuda
