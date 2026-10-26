#pragma once

#include <memory>
#include <string_view>

#include "potential/potential.hpp"

namespace spodvr {

// Creates the PES selected by its system ID (case-insensitive, like the Fortran .SEQ.).
// Throws InputError for an unknown ID. Fortran: H2SurfaceSystem.f90 performActionH2SurfaceSystem
[[nodiscard]] std::unique_ptr<Potential> make_potential(std::string_view system_id);

// True if the system ID names an implemented PES. Fortran: systemIsImplemented
[[nodiscard]] bool potential_is_implemented(std::string_view system_id);

}  // namespace spodvr
