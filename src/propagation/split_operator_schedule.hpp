#pragma once

#include <vector>

#include "propagation/kinetic_operators.hpp"

namespace spodvr {

// One T(dt') block of the split-operator step: T over `step` (half or full time step), followed
// by V unless `kinetic_only`. Fortran: SplitOperator.f90 KineticAndPotentialOp (lastStep)
struct SplitOperatorSubStep {
  TimeStep step;
  bool kinetic_only;
};

// The T blocks of one analysis step of n time steps, [T(dt/2) V(dt) T(dt/2)]^n with adjacent half
// steps of T merged: T/2 V, then n - 1 times T V, then T/2 without V. With n <= 1 the two half
// steps alone. Both backends run this list, so the operator order exists once (docs/cuda.md, "One
// method, two executions"). Fortran: SplitOperator.f90 Propagate
[[nodiscard]] std::vector<SplitOperatorSubStep> analysis_step_schedule(int steps_per_analysis);

}  // namespace spodvr
