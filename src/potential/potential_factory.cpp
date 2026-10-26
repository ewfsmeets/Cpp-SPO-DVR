// Selection of the PES by SystemID. Fortran: H2SurfaceSystem.f90
#include "potential/potential_factory.hpp"

#include <string>

#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "potential/hard_wall_potential.hpp"

namespace spodvr {

bool potential_is_implemented(std::string_view system_id) {
  return fortran_compare::seq(system_id, "HardWallPES");
}

std::unique_ptr<Potential> make_potential(std::string_view system_id) {
  if (fortran_compare::seq(system_id, "HardWallPES")) return std::make_unique<HardWallPotential>();
  throw InputError("SystemID", "unknown system ID '" + std::string(system_id) +
                                   "' (implemented: HardWallPES)");
}

}  // namespace spodvr
