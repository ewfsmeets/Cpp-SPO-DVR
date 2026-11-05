// Potential stage of the GPU backend: tables, exp(i A) on the host, the slab loop.
// CPU twin: propagation/potential_operator.cpp
#include "cuda/device_potential_operator.hpp"

#include <complex>
#include <stdexcept>
#include <string>

#include "propagation/potential_operator.hpp"

namespace spodvr::cuda {

namespace {

using cplx = std::complex<double>;

// exp(i A) with izr fastest: (izr, itp, irep) from A on (itp, irep, izr). Same cos/sin as the CPU
// table (PotentialOperator), so the values are bitwise equal; only the order differs.
std::vector<cplx> exp_potential_device_order(std::span<const double> a, std::int64_t n_tp, std::int64_t n_rep,
                                             std::int64_t n_zr) {
  std::vector<cplx> out(a.size());
#pragma omp parallel for schedule(static)
  for (std::int64_t k = 0; k < n_tp * n_rep; ++k)  // k = itp + n_tp * irep
    for (std::int64_t izr = 0; izr < n_zr; ++izr) out[izr + n_zr * k] = exp_i(a[k + n_tp * n_rep * izr]);
  return out;
}

}  // namespace

std::int64_t DevicePotentialOperator::slab_size(const ScatteringShape& s, int nz_potential, int n_theta, int n_phi,
                                                std::size_t scratch_limit) {
  const std::int64_t points = std::int64_t(nz_potential) * s.nr;
  const std::size_t bytes_per_point = sizeof(Complex) * std::size_t(s.nx) * s.ny * n_theta * n_phi;
  for (std::int64_t n = points; n > 1; --n)
    if (points % n == 0 && n * bytes_per_point <= scratch_limit) return n;
  return 1;
}

std::size_t DevicePotentialOperator::memory_bytes(const ScatteringShape& s, int nz_potential, int n_theta, int n_phi,
                                                  int n_representatives, std::size_t scratch_limit) {
  const std::size_t n_tp = std::size_t(n_theta) * n_phi;
  const std::size_t table = sizeof(Complex) * n_tp * n_representatives * std::size_t(nz_potential) * s.nr;
  const std::size_t scratch =
      sizeof(Complex) * std::size_t(s.nx) * s.ny * n_tp * slab_size(s, nz_potential, n_theta, n_phi, scratch_limit);
  return table + scratch;
}

DevicePotentialOperator::DevicePotentialOperator(const ScatteringShape& s, int nz_potential, int initial_channel,
                                                 const PropagationTables& tables, const AngularTransform& transform,
                                                 std::span<const double> premultiplied_potential,
                                                 const GridSymmetry& symmetry, std::size_t max_fft_work,
                                                 std::size_t scratch_limit) {
  const int n_theta = transform.n_theta();
  const int n_phi = transform.n_phi();
  const std::int64_t n_tp = std::int64_t(n_theta) * n_phi;
  const std::int64_t n_rep = symmetry.n_representatives();
  const std::int64_t n_zr_potential = std::int64_t(nz_potential) * s.nr;
  if (premultiplied_potential.size() != std::size_t(n_tp * n_rep * n_zr_potential))
    throw std::invalid_argument("DevicePotentialOperator: potential has the wrong size");

  // Tables shared with the CPU operator.
  const PotentialStepFactors factors = potential_step_factors(s, tables);
  z_optical_ = DeviceBuffer<double>(std::span<const double>(tables.z_optical));
  before_interaction_ = upload_complex(factors.before_interaction);
  gas_phase_step_ = upload_complex(factors.gas_phase_step);
  rotation_half_ = upload_complex(tables.rotation.half);
  specular_to_scattering_ = upload_complex(tables.specular_to_scattering);
  scattering_to_specular_ = upload_complex(tables.scattering_to_specular);

  // The angular transform: mJ blocks and the (theta, J) matrices.
  std::vector<int> first, size, slot, offset;
  int matrix_offset = 0;
  for (const AngularTransform::MjBlock& b : transform.mj_blocks()) {
    if (b.fbr.size() > max_j_per_mj_block)
      throw std::invalid_argument("GPU backend: more than " + std::to_string(max_j_per_mj_block) +
                                  " J values in one mJ block are not supported");
    first.push_back(b.fbr.first);
    size.push_back(b.fbr.size());
    slot.push_back(b.phi_slot);
    offset.push_back(matrix_offset);
    matrix_offset += n_theta * b.fbr.size();
  }
  block_first_ = DeviceBuffer<int>(std::span<const int>(first));
  block_size_ = DeviceBuffer<int>(std::span<const int>(size));
  block_phi_slot_ = DeviceBuffer<int>(std::span<const int>(slot));
  block_matrix_offset_ = DeviceBuffer<int>(std::span<const int>(offset));
  matrix_ = DeviceBuffer<double>(std::span<const double>(transform.matrix()));

  // exp(i A) and the symmetry maps (representative and phi map of every (x, y) point).
  exp_potential_ = upload_complex(exp_potential_device_order(premultiplied_potential, n_tp, n_rep, n_zr_potential));
  std::vector<int> representative, operation, phi_map;
  for (int ixy = 0; ixy < symmetry.n_xy(); ++ixy) {
    representative.push_back(symmetry.representative_of(ixy));
    operation.push_back(symmetry.operation_of(ixy));
  }
  for (std::size_t iop = 0; iop < symmetry.operations().size(); ++iop)
    for (const int phi : symmetry.phi_index(static_cast<int>(iop))) phi_map.push_back(phi);
  representative_of_ = DeviceBuffer<int>(std::span<const int>(representative));
  operation_of_ = DeviceBuffer<int>(std::span<const int>(operation));
  phi_index_ = DeviceBuffer<int>(std::span<const int>(phi_map));

  // Slab, scratch and the phi FFTs: one transform of length n_phi per (column, theta), the phi
  // values n_columns * n_theta apart (scratch layout (col, theta, phi)).
  const std::int64_t n_slab = slab_size(s, nz_potential, n_theta, n_phi, scratch_limit);
  const std::int64_t n_columns = n_slab * s.nx * s.ny;
  if (n_columns * n_tp >= (std::int64_t(1) << 31))  // the kernels index within a slab with int
    throw std::invalid_argument("GPU backend: the DVR scratch of one (Z, r) point exceeds 2^31 values");
  dvr_ = DeviceBuffer<Complex>(std::size_t(n_columns * n_tp));
  if (n_phi > 1) {
    const DeviceFftShape phi{{n_phi}, n_columns * n_theta, 1, n_columns * n_theta};
    to_dvr_fft_ = std::make_unique<DeviceFftPlan>(phi, FftSign::Plus, max_fft_work);  // Fortran: FFTForFBRToDVRInPhi
    to_fbr_fft_ = std::make_unique<DeviceFftPlan>(phi, FftSign::Minus, max_fft_work);  // Fortran: FFTForDVRToFBRInPhi
  }

  args_.psi = DeviceLayout{std::int64_t(s.nr) * s.nz, s.n_rot, std::int64_t(s.nx) * s.ny};
  args_.nr = s.nr;
  args_.n_theta = n_theta;
  args_.n_phi = n_phi;
  args_.initial_channel = initial_channel;
  args_.n_zr_potential = n_zr_potential;
  args_.izr_first = 0;
  args_.n_slab = n_slab;
  args_.z_optical = z_optical_.data();
  args_.before_interaction = before_interaction_.data();
  args_.gas_phase_step = gas_phase_step_.data();
  args_.rotation_half = rotation_half_.data();
  args_.specular_to_scattering = specular_to_scattering_.data();
  args_.scattering_to_specular = scattering_to_specular_.data();
  args_.n_mj_blocks = static_cast<int>(first.size());
  args_.mj_max = transform.mj_max();
  args_.block_first = block_first_.data();
  args_.block_size = block_size_.data();
  args_.block_phi_slot = block_phi_slot_.data();
  args_.block_matrix_offset = block_matrix_offset_.data();
  args_.matrix = matrix_.data();
  args_.exp_potential = exp_potential_.data();
  args_.representative_of = representative_of_.data();
  args_.operation_of = operation_of_.data();
  args_.phi_index = phi_index_.data();
}

std::vector<DeviceFftPlan*> DevicePotentialOperator::fft_plans() {
  if (!to_dvr_fft_) return {};
  return {to_dvr_fft_.get(), to_fbr_fft_.get()};
}

// Fortran: Potential.f90 PerformEXPPotentialOp for every (Z, r) point
void DevicePotentialOperator::apply(DeviceWavePacket& psi) const {
  Complex* scattering = psi.scattering.data();
  Complex* specular = psi.specular.data();
  apply_gas_phase_step(args_, scattering);

  PotentialArgs slab = args_;
  Complex* dvr = dvr_.data();
  for (slab.izr_first = 0; slab.izr_first < args_.n_zr_potential; slab.izr_first += args_.n_slab) {
    absorb_and_transform_to_dvr(slab, scattering, specular, dvr);
    if (to_dvr_fft_) to_dvr_fft_->execute(dvr);
    multiply_by_exp_potential(slab, dvr);
    if (to_fbr_fft_) to_fbr_fft_->execute(dvr);
    transform_to_fbr_and_rotate(slab, dvr, scattering);
    project_onto_specular(slab, scattering, specular);
  }
}

}  // namespace spodvr::cuda
