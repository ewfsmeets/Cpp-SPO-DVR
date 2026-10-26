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
  // Along Z: nr transforms with stride nr; along r: nzsp contiguous transforms.
  const FftLayout z_layout{nz_, nr_, nr_, 1};
  z_to_momentum_ = fft.plan_complex(z_layout, FftSign::Minus, tables.fft_scale_z);
  z_to_coordinate_ = fft.plan_complex(z_layout, FftSign::Plus, 1.0 / (nz_ * tables.fft_scale_z));
  const FftLayout r_layout{nr_, nz_, 1, nr_};
  r_to_momentum_ = fft.plan_complex(r_layout, FftSign::Minus, tables.fft_scale_r);
  r_to_coordinate_ = fft.plan_complex(r_layout, FftSign::Plus, 1.0 / (nr_ * tables.fft_scale_r));
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
  z_to_momentum_->execute(psi.values);
  if (nr_ != 1) r_to_momentum_->execute(psi.values);
  for (std::size_t i = 0; i < psi.values.size(); ++i) psi.values[i] = factors[i] * psi.values[i];
  z_to_coordinate_->execute(psi.values);
  if (nr_ != 1) r_to_coordinate_->execute(psi.values);
}

// Free molecule beyond the interaction region: rot(dt/2) exp(-i vasy dt) rot(dt/2) for the
// initial rotational state. Fortran: GasPhase.f90 PerformEXPGasPhaseOnSpecGrid
void SpecularOperator::apply_gas_phase(SpecularWave& psi) const {
  const int n_rot = static_cast<int>(tables_.rotation.half.size()) / nr_;
  for (int iz = nz_potential_; iz < nz_; ++iz)
    for (int ir = 0; ir < nr_; ++ir) {
      const cplx rot = tables_.rotation.half[initial_channel_ + n_rot * ir];
      cplx& v = psi(ir, iz);
      v = v * rot;
      v = v * tables_.gas_phase[ir];
      v = v * rot;
    }
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
