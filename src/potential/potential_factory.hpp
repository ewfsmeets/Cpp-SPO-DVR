#pragma once

#include <memory>
#include <string_view>

#include "input/config.hpp"
#include "potential/potential.hpp"

namespace spodvr {

// Creates the PES selected by input.system_id (case-insensitive, like the Fortran .SEQ.).
// Throws InputError for an unknown ID. Fortran: H2SurfaceSystem.f90 performActionH2SurfaceSystem
[[nodiscard]] std::unique_ptr<Potential> make_potential(const InputParameters& input);

// True if the system ID names an implemented PES. Fortran: systemIsImplemented
[[nodiscard]] bool potential_is_implemented(std::string_view system_id);

}  // namespace spodvr
