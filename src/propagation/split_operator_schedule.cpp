// The order of the T and V blocks in one analysis step, shared by the CPU and GPU backends.
// Fortran: SplitOperator.f90 Propagate
#include "propagation/split_operator_schedule.hpp"

namespace spodvr {

std::vector<SplitOperatorSubStep> analysis_step_schedule(int steps_per_analysis) {
  std::vector<SplitOperatorSubStep> schedule;
  schedule.push_back({TimeStep::Half, false});
  for (int step = 1; step < steps_per_analysis; ++step) schedule.push_back({TimeStep::Full, false});
  schedule.push_back({TimeStep::Half, true});
  return schedule;
}

}  // namespace spodvr
