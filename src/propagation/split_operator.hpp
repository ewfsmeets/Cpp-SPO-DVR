#pragma once

#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"
#include "propagation/specular_operator.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Split-operator propagation over one analysis step of n time steps,
//   [T(dt/2) V(dt) T(dt/2)]^n = T(dt/2) V (T V)^(n-1) T(dt/2),
// with adjacent half steps of T merged. One T block: specular grid operations, T_Z and T_r on
// the scattering grid, then T_xy per (Z, r) point; V follows T except after the final half step.
// Fortran: SplitOperator.f90 Propagate
class SplitOperatorPropagator {
 public:
  SplitOperatorPropagator(int steps_per_analysis, const ScatteringKineticOperator& kinetic,
                          const PotentialOperator& potential, const SpecularOperator& specular)
      : steps_per_analysis_(steps_per_analysis), kinetic_(kinetic), potential_(potential), specular_(specular) {}

  void propagate_analysis_step(WavePacket& psi) const;

 private:
  int steps_per_analysis_;  // Fortran: nanspo
  const ScatteringKineticOperator& kinetic_;
  const PotentialOperator& potential_;
  const SpecularOperator& specular_;
};

}  // namespace spodvr
