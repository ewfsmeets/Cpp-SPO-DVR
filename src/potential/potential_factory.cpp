// Selection of the PES by SystemID. Fortran: H2SurfaceSystem.f90
#include "potential/potential_factory.hpp"

#include <format>
#include <string>
#include <string_view>

#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "potential/crp_potential.hpp"
#include "potential/dai_light_cu111_potential.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/persson_cu111_potential.hpp"
#include "potential/scm/scm_crp_potential.hpp"

namespace spodvr {

namespace {
// The implemented SystemIDs (Fortran names; compared case-insensitively).
constexpr std::string_view hard_wall_id = "HardWallPES";
constexpr std::string_view dai_light_id = "H2Cu111DaiLight";
constexpr std::string_view persson_id = "H2Cu111Persson";
constexpr std::string_view crp_id = "CRPlibrary";  // a CRP data set (CRPlibDataSet)
constexpr std::string_view scm_id = "CRPDistrtSurface";  // a CRP data set plus the SCM correction (upstream 5.24)
}  // namespace

bool potential_is_implemented(std::string_view system_id) {
  for (const std::string_view id : {hard_wall_id, dai_light_id, persson_id, crp_id, scm_id})
    if (fortran_compare::seq(system_id, id)) return true;
  return false;
}

std::unique_ptr<Potential> make_potential(const InputParameters& input) {
  const std::string_view system_id = input.system_id;
  if (fortran_compare::seq(system_id, hard_wall_id)) return std::make_unique<HardWallPotential>();
  if (fortran_compare::seq(system_id, dai_light_id)) return std::make_unique<DaiLightCu111Potential>();
  if (fortran_compare::seq(system_id, persson_id)) return std::make_unique<PerssonCu111Potential>();
  if (fortran_compare::seq(system_id, crp_id)) return std::make_unique<CrpPotential>(input.crp_dataset);
  // Case-insensitive like every SystemID; the Fortran reads the SCM keys only for this exact
  // spelling (SCM D1). Fortran: H2SurfaceSystem.f90 (CRPDistrtSurface branch)
  if (fortran_compare::seq(system_id, scm_id)) return std::make_unique<ScmCrpPotential>(input.crp_dataset, input.scm);
  throw InputError("SystemID", std::format("unknown system ID '{}' (implemented: {}, {}, {}, {}, {})", system_id,
                                           hard_wall_id, dai_light_id, persson_id, crp_id, scm_id));
}

}  // namespace spodvr
