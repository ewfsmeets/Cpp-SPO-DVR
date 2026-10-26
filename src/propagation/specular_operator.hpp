#pragma once

#include <memory>

#include "core/index_range.hpp"
#include "fft/fft.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Operations on the specular (incoming channel) grid in one time step: exp(-i T dt) by FFT in
// Z and r (T includes the parallel energy), then, except on the last half step, the gas-phase
// propagator rot(dt/2) . exp(-i vasy dt) . rot(dt/2) beyond the interaction region and the
// absorbing potentials. Fortran: SplitOperator.f90 SpecularChannelOps; KineticOp.f90
// PerformEXPKinOpSpec; GasPhase.f90 PerformEXPGasPhaseOnSpecGrid; Absorption.f90
// PerformEXPOpticalOpSpec
class SpecularOperator {
 public:
  SpecularOperator(int nr, int nz_specular, int nz_potential, int initial_channel, const IndexRange& z_absorber,
                   const IndexRange& r_absorber, const PropagationTables& tables, const Fft& fft);

  void apply(SpecularWave& psi, TimeStep step, bool kinetic_only) const;

 private:
  void apply_kinetic(SpecularWave& psi, TimeStep step) const;
  void apply_gas_phase(SpecularWave& psi) const;
  void apply_absorbers(SpecularWave& psi) const;

  int nr_;
  int nz_;
  int nz_potential_;
  int initial_channel_;
  IndexRange z_absorber_;  // Fortran: nschsp .. nschmx
  IndexRange r_absorber_;  // Fortran: nrchsp .. nrchmx
  const PropagationTables& tables_;
  std::unique_ptr<ComplexFftPlan> z_to_momentum_, z_to_coordinate_;  // Fortran: BackwSpecFFTz / ForwSpecFFTz
  std::unique_ptr<ComplexFftPlan> r_to_momentum_, r_to_coordinate_;  // Fortran: BackwSpecFFTr / ForwSpecFFTr
};

}  // namespace spodvr
