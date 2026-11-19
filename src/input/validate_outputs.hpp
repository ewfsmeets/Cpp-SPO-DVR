#pragma once

#include "input/config.hpp"

namespace spodvr {

// The checks of the plot-data keys (writeDensities, ...; docs/input_schema.md), which need the run mode
// and the derived parameters: called after miscellaneous_setup. C++ only (no Fortran counterpart); kept out
// of validate_input, whose checks mirror InputNew.f90 one by one.
void validate_output_keys(const Config& config);

}  // namespace spodvr
