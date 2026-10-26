// Split-operator propagation over one analysis step. Fortran: SplitOperator.f90 Propagate
#include "propagation/split_operator.hpp"

namespace spodvr {

void SplitOperatorPropagator::propagate_analysis_step(WavePacket& psi) const {
  // First half step of T, then V. Fortran: SpecularChannelOps, KineticOpAlongZrOnScatGrid,
  // KineticAndPotentialOp (HALF_TIME_STEP)
  specular_.apply(psi.specular, TimeStep::Half, false);
  kinetic_.apply_zr(psi.scattering, TimeStep::Half);
  kinetic_.apply_xy(psi.scattering, TimeStep::Half);
  potential_.apply(psi);

  // Merged half steps: full T, then V.
  for (int step = 1; step < steps_per_analysis_; ++step) {
    specular_.apply(psi.specular, TimeStep::Full, false);
    kinetic_.apply_zr(psi.scattering, TimeStep::Full);
    kinetic_.apply_xy(psi.scattering, TimeStep::Full);
    potential_.apply(psi);
  }

  // Final half step of T only (Fortran: lastStep = .TRUE.).
  specular_.apply(psi.specular, TimeStep::Half, true);
  kinetic_.apply_zr(psi.scattering, TimeStep::Half);
  kinetic_.apply_xy(psi.scattering, TimeStep::Half);
}

}  // namespace spodvr
