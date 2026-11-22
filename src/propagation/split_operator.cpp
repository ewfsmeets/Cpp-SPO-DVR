// Split-operator propagation over one analysis step. Fortran: SplitOperator.f90 Propagate
#include "propagation/split_operator.hpp"

#include <vector>

#include "core/aligned_vector.hpp"

namespace spodvr {

// T/2 V, (T V) n - 1 times, T/2 (analysis_step_schedule).
void SplitOperatorPropagator::propagate_analysis_step(WavePacket& psi) const {
  for (const SplitOperatorSubStep& sub_step : schedule_)
    kinetic_and_potential_step(psi, sub_step.step, sub_step.kinetic_only);
}

// Fortran: SpecularChannelOps, KineticOpAlongZrOnScatGrid, KineticAndPotentialOp
void SplitOperatorPropagator::kinetic_and_potential_step(WavePacket& psi, TimeStep step, bool kinetic_only) const {
  specular_.apply(psi.specular, step, kinetic_only);
  kinetic_.apply_zr(psi.scattering, step);
  apply_xy_and_potential(psi, step, kinetic_only);
}

// Fortran: SplitOperator.f90 KineticAndPotentialOp. Processed in groups of consecutive (Z, r)
// points: T_xy for the group at once, then V point by point, while the data are in cache.
// Groups with Z < zasymp (angular transform + exp(-iV dt)) cost far more than gas-phase groups
// (one scalar factor); the schedule (dynamic,1 or static,1, core/loop_schedule.hpp) balances
// them (the Fortran uses SCHEDULE(DYNAMIC) here). Results do not depend on the schedule: every
// point is independent and each thread has its own work arrays.
void SplitOperatorPropagator::apply_xy_and_potential(WavePacket& psi, TimeStep step, bool kinetic_only) const {
  const int group_points = kinetic_.group_points();
  const std::size_t block_size = psi.scattering.block_size();
#pragma omp parallel
  {
    ComplexVector xy_work(kinetic_.group_work_size());
    ComplexVector potential_work(potential_.work_size());
#pragma omp for schedule(runtime)
    for (int g = 0; g < psi.scattering.n_blocks() / group_points; ++g) {
      const int first = g * group_points;
      kinetic_.apply_xy_group(std::span(psi.scattering.values).subspan(block_size * first, block_size * group_points),
                              step, xy_work);
      if (kinetic_only) continue;
      for (int izr = first; izr < first + group_points; ++izr)
        potential_.apply_block(izr, psi.scattering.block(izr), psi.specular, potential_work);
    }
  }
}

}  // namespace spodvr
