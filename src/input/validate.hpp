#pragma once

#include "input/config.hpp"

namespace spodvr {

// All consistency checks of the Fortran input reader; throws InputError naming the offending
// parameter on the first violation. zasymp is the (possibly switched) asymptotic Z.
// Checks that need the grids (Setup.f90) are done when the grids are built.
// Fortran: InputNew.f90 ValidateInput
void validate_input(const InputParameters& p, const SystemParameters& system, double zasymp);

}  // namespace spodvr
