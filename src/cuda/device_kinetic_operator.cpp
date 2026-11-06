// Kinetic stages of the GPU backend: cuFFT plans and factor tables, launch order.
// CPU twins: propagation/kinetic_operators.cpp, propagation/specular_operator.cpp
#include "cuda/device_kinetic_operator.hpp"

#include <algorithm>
#include <complex>

#include "cuda/kernels.hpp"

namespace spodvr::cuda {

namespace {

using cplx = std::complex<double>;

bool is_one(const cplx& c) { return c.real() == 1.0 && c.imag() == 0.0; }

}  // namespace

DeviceKineticOperator::DeviceKineticOperator(const ScatteringShape& s, const PropagationTables& tables,
                                             std::size_t max_fft_work)
    : shape_(s), shift_(!std::all_of(tables.shift_forward.begin(), tables.shift_forward.end(), is_one)) {
  const std::int64_t n_zr = std::int64_t(s.nr) * s.nz;
  const std::int64_t n_xy = std::int64_t(s.nx) * s.ny;

  // (Z, r): every (irot, ixy) is a contiguous block of n_zr values (Z outer, r inner).
  const double scale_z = tables.fft_scale_z;
  const double scale_r = tables.fft_scale_r;
  zr_scales_ = {scale_z * scale_r, (1.0 / (s.nz * scale_z)) * (1.0 / (s.nr * scale_r))};
  zr_full_ = upload_complex(tables.zr_kinetic.full);
  zr_half_ = upload_complex(tables.zr_kinetic.half);
  const DeviceFftShape zr{{s.nz, s.nr}, 1, n_zr, s.n_rot * n_xy};
  zr_to_momentum_ = std::make_unique<DeviceFftPlan>(zr, FftSign::Minus, max_fft_work);
  zr_to_coordinate_ = std::make_unique<DeviceFftPlan>(zr, FftSign::Plus, max_fft_work);

  // (x, y): element (ix, iy) of batch b = izr + n_zr * irot at b + n_zr*n_rot*(iy + ny*ix). A single
  // grid point in x or y is not transformed (nor scaled), as on the CPU.
  const double scale_x = s.nx > 1 ? tables.fft_scale_x : 1.0;
  const double scale_y = s.ny > 1 ? tables.fft_scale_y : 1.0;
  xy_scales_ = {scale_x * scale_y, (1.0 / (s.nx * scale_x)) * (1.0 / (s.ny * scale_y))};
  xy_full_ = upload_complex(tables.xy_kinetic.full);
  xy_half_ = upload_complex(tables.xy_kinetic.half);
  std::vector<std::int64_t> lengths;
  if (s.nx > 1) lengths.push_back(s.nx);
  if (s.ny > 1) lengths.push_back(s.ny);
  if (!lengths.empty()) {
    const DeviceFftShape xy{lengths, n_zr * s.n_rot, 1, n_zr * s.n_rot};
    xy_to_momentum_ = std::make_unique<DeviceFftPlan>(xy, FftSign::Minus, max_fft_work);
    xy_to_coordinate_ = std::make_unique<DeviceFftPlan>(xy, FftSign::Plus, max_fft_work);
  }
  if (shift_) {
    shift_forward_ = upload_complex(tables.shift_forward);
    shift_backward_ = upload_complex(tables.shift_backward);
  }
}

std::vector<DeviceFftPlan*> DeviceKineticOperator::fft_plans() {
  std::vector<DeviceFftPlan*> plans{zr_to_momentum_.get(), zr_to_coordinate_.get()};
  if (xy_to_momentum_) plans.insert(plans.end(), {xy_to_momentum_.get(), xy_to_coordinate_.get()});
  return plans;
}

// exp(-i (T_Z + T_r) dt) and exp(-i T_xy dt) with one multiplication pass instead of two: both
// forward transforms, both factors, both inverse transforms. T_Z + T_r and T_xy commute, and so
// do the momentum shift (a factor per ixy) and the (Z, r) transform, so this equals the CPU's
// T_zr stage followed by its T_xy stage up to rounding (a methodology exception, docs/cuda.md).
// Fortran: KineticOp.f90 PerformEXPKinOpZScat/RScat, SplitOperator.f90 PerformEXPKineticOp
void DeviceKineticOperator::apply(DeviceWavePacket& psi, TimeStep step) const {
  const bool full = step == TimeStep::Full;
  const DeviceBuffer<Complex>& zr_factors = full ? zr_full_ : zr_half_;
  const DeviceBuffer<Complex>& xy_factors = full ? xy_full_ : xy_half_;
  Complex* data = psi.scattering.data();
  const std::int64_t per_xy = psi.layout.n_zr * psi.layout.n_rot;  // values of one ixy
  const std::int64_t n = psi.layout.size();
  if (shift_) multiply_by_blocked_factors(data, shift_forward_.data(), FftScales{}, per_xy, n);
  zr_to_momentum_->execute(data);
  if (xy_to_momentum_) xy_to_momentum_->execute(data);
  multiply_by_kinetic_factors(data, psi.layout, zr_factors.data(), zr_scales_, xy_factors.data(), xy_scales_);
  if (xy_to_coordinate_) xy_to_coordinate_->execute(data);
  zr_to_coordinate_->execute(data);
  if (shift_) multiply_by_blocked_factors(data, shift_backward_.data(), FftScales{}, per_xy, n);
}

DeviceSpecularOperator::DeviceSpecularOperator(int nr, int nz_specular, int nz_potential, int initial_channel,
                                               const IndexRange& z_absorber, const IndexRange& r_absorber,
                                               const PropagationTables& tables, std::size_t max_fft_work)
    : nr_(nr),
      nz_(nz_specular),
      nz_potential_(nz_potential),
      z_absorber_(z_absorber),
      r_absorber_(r_absorber),
      gas_phase_step_(upload_complex(specular_gas_phase_factors(nr, initial_channel, tables))),
      z_optical_(std::span<const double>(tables.z_specular_optical)),
      r_optical_(std::span<const double>(tables.r_optical)) {
  const double scale_z = tables.fft_scale_z;
  const double scale_r = tables.fft_scale_r;
  scales_ = {scale_z * scale_r, (1.0 / (nz_ * scale_z)) * (1.0 / (nr_ * scale_r))};
  kinetic_full_ = upload_complex(tables.specular_kinetic.full);
  kinetic_half_ = upload_complex(tables.specular_kinetic.half);
  // The wave function (ir, iz) as one 2D transform: Z outer, r inner.
  const DeviceFftShape layout{{nz_, nr_}, 1, std::int64_t(nz_) * nr_, 1};
  to_momentum_ = std::make_unique<DeviceFftPlan>(layout, FftSign::Minus, max_fft_work);
  to_coordinate_ = std::make_unique<DeviceFftPlan>(layout, FftSign::Plus, max_fft_work);
}

std::vector<DeviceFftPlan*> DeviceSpecularOperator::fft_plans() { return {to_momentum_.get(), to_coordinate_.get()}; }

// Fortran: SplitOperator.f90 SpecularChannelOps (kinetic_only on the last half step)
void DeviceSpecularOperator::apply(DeviceWavePacket& psi, TimeStep step, bool kinetic_only) const {
  const DeviceBuffer<Complex>& factors = step == TimeStep::Full ? kinetic_full_ : kinetic_half_;
  Complex* data = psi.specular.data();
  const std::int64_t n = std::int64_t(nr_) * nz_;
  to_momentum_->execute(data);
  multiply_by_periodic_factors(data, factors.data(), scales_, n, n);
  to_coordinate_->execute(data);
  if (kinetic_only) return;
  const SpecularStepArgs args{nr_,
                              nz_,
                              nz_potential_,
                              gas_phase_step_.data(),
                              z_absorber_.first,
                              z_absorber_.last,
                              z_optical_.data(),
                              r_absorber_.first,
                              r_absorber_.last,
                              r_optical_.data()};
  apply_specular_gas_phase_and_absorbers(args, data);
}

}  // namespace spodvr::cuda
