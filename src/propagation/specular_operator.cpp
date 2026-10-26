// Propagation of the incoming channel on the specular grid.
// Fortran: SplitOperator.f90 SpecularChannelOps and the routines it calls
#include "propagation/specular_operator.hpp"

namespace spodvr {

namespace {
using cplx = std::complex<double>;
}  // namespace

// Fortran: KineticOp.f90 SetupFFTs (nproj > 0)
SpecularOperator::SpecularOperator(int nr, int nz_specular, int nz_potential, int initial_channel,
                                   const IndexRange& z_absorber, const IndexRange& r_absorber,
                                   const PropagationTables& tables, const Fft& fft)
    : nr_(nr),
      nz_(nz_specular),
      nz_potential_(nz_potential),
      initial_channel_(initial_channel),
      z_absorber_(z_absorber),
      r_absorber_(r_absorber),
      tables_(tables) {
  // The wave function (ir, iz) as one 2D transform: Z outer (stride nr), r inner (stride 1).
  const FftLayout2D layout{nz_, nr_, nr_, 1, 1, 1};
  const double scale_z = tables.fft_scale_z;
  const double scale_r = tables.fft_scale_r;
  to_momentum_ = fft.plan_complex_2d(layout, FftSign::Minus, scale_z * scale_r);
  to_coordinate_ = fft.plan_complex_2d(layout, FftSign::Plus, (1.0 / (nz_ * scale_z)) * (1.0 / (nr_ * scale_r)));

  const int n_rot = static_cast<int>(tables.rotation.half.size()) / nr_;
  for (int ir = 0; ir < nr_; ++ir) {
    const cplx rot = tables.rotation.half[initial_channel_ + n_rot * ir];
    gas_phase_step_.push_back(rot * tables.gas_phase[ir] * rot);
  }
}

// Fortran: SplitOperator.f90 SpecularChannelOps (kinetic_only on the last half step)
void SpecularOperator::apply(SpecularWave& psi, TimeStep step, bool kinetic_only) const {
  apply_kinetic(psi, step);
  if (kinetic_only) return;
  apply_gas_phase(psi);
  apply_absorbers(psi);
}

// Fortran: KineticOp.f90 PerformEXPKinOpSpec
void SpecularOperator::apply_kinetic(SpecularWave& psi, TimeStep step) const {
  const std::vector<cplx>& factors = step == TimeStep::Full ? tables_.specular_kinetic.full : tables_.specular_kinetic.half;
  to_momentum_->execute(psi.values);
  for (std::size_t i = 0; i < psi.values.size(); ++i) psi.values[i] = factors[i] * psi.values[i];
  to_coordinate_->execute(psi.values);
}

// Free molecule beyond the interaction region: rot(dt/2) exp(-i vasy dt) rot(dt/2) for the
// initial rotational state. Fortran: GasPhase.f90 PerformEXPGasPhaseOnSpecGrid
void SpecularOperator::apply_gas_phase(SpecularWave& psi) const {
  for (int iz = nz_potential_; iz < nz_; ++iz)
    for (int ir = 0; ir < nr_; ++ir) psi(ir, iz) = gas_phase_step_[ir] * psi(ir, iz);
}

// Fortran: Absorption.f90 PerformEXPOpticalOpSpec (nproj = 2)
void SpecularOperator::apply_absorbers(SpecularWave& psi) const {
  for (int iz = z_absorber_.first; iz <= z_absorber_.last; ++iz)
    for (int ir = 0; ir < nr_; ++ir) psi(ir, iz) = tables_.z_specular_optical[iz] * psi(ir, iz);
  if (nr_ == 1) return;
  for (int iz = 0; iz < nz_; ++iz)
    for (int ir = r_absorber_.first; ir <= r_absorber_.last; ++ir) psi(ir, iz) = tables_.r_optical[ir] * psi(ir, iz);
}

}  // namespace spodvr
