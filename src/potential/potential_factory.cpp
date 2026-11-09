// Selection of the PES by SystemID. Fortran: H2SurfaceSystem.f90
#include "potential/potential_factory.hpp"

#include <string>

#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "potential/dai_light_cu111_potential.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/persson_cu111_potential.hpp"

namespace spodvr {

namespace {
// The implemented SystemIDs (Fortran names; compared case-insensitively).
constexpr const char* hard_wall_id = "HardWallPES";
constexpr const char* dai_light_id = "H2Cu111DaiLight";
constexpr const char* persson_id = "H2Cu111Persson";
}  // namespace

bool potential_is_implemented(std::string_view system_id) {
  for (const char* id : {hard_wall_id, dai_light_id, persson_id})
    if (fortran_compare::seq(system_id, id)) return true;
  return false;
}

std::unique_ptr<Potential> make_potential(const InputParameters& input) {
  const std::string_view system_id = input.system_id;
  if (fortran_compare::seq(system_id, hard_wall_id)) return std::make_unique<HardWallPotential>();
  if (fortran_compare::seq(system_id, dai_light_id)) return std::make_unique<DaiLightCu111Potential>();
  if (fortran_compare::seq(system_id, persson_id)) return std::make_unique<PerssonCu111Potential>();
  throw InputError("SystemID", "unknown system ID '" + std::string(system_id) + "' (implemented: " + hard_wall_id +
                                   ", " + dai_light_id + ", " + persson_id + ")");
}

}  // namespace spodvr
