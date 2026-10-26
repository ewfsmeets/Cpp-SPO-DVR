// Split-operator propagation over one analysis step. Fortran: SplitOperator.f90 Propagate
#include "propagation/split_operator.hpp"

#include <vector>

#include "core/aligned_vector.hpp"

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

// Fortran: SplitOperator.f90 KineticAndPotentialOp. Processed per Z value: T_xy for its nr
// points at once, then V point by point, while the data are in cache.
void SplitOperatorPropagator::apply_xy_and_potential(WavePacket& psi, TimeStep step, bool kinetic_only) const {
  const int nr = psi.specular.n_r();
  const std::size_t block_size = psi.scattering.block_size();
#pragma omp parallel
  {
    ComplexVector xy_work(kinetic_.group_work_size());
    ComplexVector potential_work(potential_.work_size());
#pragma omp for schedule(static)
    for (int iz = 0; iz < psi.scattering.n_blocks() / nr; ++iz) {
      kinetic_.apply_xy_group(std::span(psi.scattering.values).subspan(block_size * nr * iz, block_size * nr), step,
                              xy_work);
      if (kinetic_only) continue;
      for (int ir = 0; ir < nr; ++ir) {
        const int izr = ir + nr * iz;
        potential_.apply_block(izr, psi.scattering.block(izr), psi.specular, potential_work);
      }
    }
  }
}

}  // namespace spodvr
