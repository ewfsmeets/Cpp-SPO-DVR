#pragma once

#include <memory>
#include <vector>

#include "core/index_range.hpp"
#include "fft/fft.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// rot(dt/2) exp(-i vasy(r) dt) rot(dt/2) of the initial rotational state, per ir: the free
// molecule beyond the interaction region on the specular grid. Fortran: GasPhase.f90
// PerformEXPGasPhaseOnSpecGrid
[[nodiscard]] std::vector<std::complex<double>> specular_gas_phase_factors(int nr, int initial_channel,
                                                                           const PropagationTables& tables);

// Operations on the specular (incoming channel) grid in one time step: exp(-i T dt) by a 2D FFT
// in (r, Z) (T includes the parallel energy; Fortran: 1D FFTs in Z, then r), then, except on the last half step, the
// gas-phase propagator rot(dt/2) . exp(-i vasy dt) . rot(dt/2) beyond the interaction region and the absorbing
// potentials. Fortran: SplitOperator.f90 SpecularChannelOps; KineticOp.f90 PerformEXPKinOpSpec; GasPhase.f90
// PerformEXPGasPhaseOnSpecGrid; Absorption.f90 PerformEXPOpticalOpSpec
class SpecularOperator {
 public:
  SpecularOperator(int nr, int nz_specular, int nz_potential, int initial_channel, const IndexRange& r_absorber,
                   const PropagationTables& tables, const Fft& fft);

  void apply(SpecularWave& psi, TimeStep step, bool kinetic_only) const;

 private:
  void apply_kinetic(SpecularWave& psi, TimeStep step) const;
  void apply_gas_phase(SpecularWave& psi) const;
  void apply_absorbers(SpecularWave& psi) const;

  int nr_;
  int nz_;
  int nz_potential_;
  int initial_channel_;
  IndexRange r_absorber_;  // Fortran: nrchsp .. nrchmx
  const PropagationTables& tables_;
  std::vector<std::complex<double>> gas_phase_step_;  // per ir: rot(dt/2) exp(-i vasy dt) rot(dt/2)
  // 2D (r, Z) FFTs. Fortran: BackwSpecFFTz/r (to momentum), ForwSpecFFTz/r (back)
  std::unique_ptr<ComplexFftPlan> to_momentum_, to_coordinate_;
};

}  // namespace spodvr
