#pragma once

#include "input/config.hpp"

namespace spodvr {

// Builds the complete configuration from the raw input: obtains the system parameters from the
// PES (built here from the input), applies the switching rule for zasymp, validates (throws InputError) and computes
// the derived parameters. Fortran: InputNew.f90 SetupInput (ValidateInput, ProcessInput) and the scalar part of
// Setup.f90 MiscellaneousSetup.
[[nodiscard]] Config build_config(const InputParameters& input);
// The same with the PES of `input` already built (make_potential), so that a run builds it once.
[[nodiscard]] Config build_config(const InputParameters& input, const Potential& potential);

// The counts process_input derives, truncated as the Fortran INT (nanspo, nave). Fortran: InputNew.f90 ProcessInput
[[nodiscard]] int steps_per_analysis(const InputParameters& p);
[[nodiscard]] int analyses_per_average(const InputParameters& p);

}  // namespace spodvr
