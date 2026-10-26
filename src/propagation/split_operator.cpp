// Split-operator propagation over one analysis step. Fortran: SplitOperator.f90 Propagate
#include "propagation/split_operator.hpp"

#include <vector>

namespace spodvr {

void SplitOperatorPropagator::propagate_analysis_step(WavePacket& psi) const {
  // First half step of T, then V. Fortran: SpecularChannelOps, KineticOpAlongZrOnScatGrid,
  // KineticAndPotentialOp (HALF_TIME_STEP)
  specular_.apply(psi.specular, TimeStep::Half, false);
  kinetic_.apply_zr(psi.scattering, TimeStep::Half);
  apply_xy_and_potential(psi, TimeStep::Half, false);

  // Merged half steps: full T, then V.
  for (int step = 1; step < steps_per_analysis_; ++step) {
    specular_.apply(psi.specular, TimeStep::Full, false);
    kinetic_.apply_zr(psi.scattering, TimeStep::Full);
    apply_xy_and_potential(psi, TimeStep::Full, false);
  }

  // Final half step of T only (Fortran: lastStep = .TRUE.).
  specular_.apply(psi.specular, TimeStep::Half, true);
  kinetic_.apply_zr(psi.scattering, TimeStep::Half);
  apply_xy_and_potential(psi, TimeStep::Half, true);
}

// Fortran: SplitOperator.f90 KineticAndPotentialOp
void SplitOperatorPropagator::apply_xy_and_potential(WavePacket& psi, TimeStep step, bool kinetic_only) const {
#pragma omp parallel
  {
    std::vector<std::complex<double>> work(potential_.work_size());
#pragma omp for schedule(static)
    for (int izr = 0; izr < psi.scattering.n_blocks(); ++izr) {
      const std::span<std::complex<double>> block = psi.scattering.block(izr);
      kinetic_.apply_xy_block(block, step);
      if (!kinetic_only) potential_.apply_block(izr, block, psi.specular, work);
    }
  }
}

}  // namespace spodvr
