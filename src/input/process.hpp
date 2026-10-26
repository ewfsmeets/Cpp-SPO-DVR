#pragma once

#include "input/config.hpp"

namespace spodvr {

// Builds the complete configuration from the raw input: obtains the system parameters from the
// PES, applies the switching rule for zasymp, validates (throws InputError) and computes the
// derived parameters.
// Fortran: InputNew.f90 SetupInput (ValidateInput, ProcessInput) and the scalar part of
// Setup.f90 MiscellaneousSetup.
[[nodiscard]] Config build_config(const InputParameters& input);

}  // namespace spodvr
