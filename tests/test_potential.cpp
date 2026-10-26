#include <gtest/gtest.h>

#include <cmath>
#include <fstream>

#include <nlohmann/json.hpp>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"

#include <omp.h>

using namespace spodvr;
using Json = nlohmann::json;

namespace {
Json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return Json::parse(file);
}

struct ReferenceSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  explicit ReferenceSetup(const Json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in) {}
};
}  // namespace

TEST(HardWall, GasPhaseZeroAtEquilibrium) { EXPECT_EQ(HardWallPotential::gas_phase(1.40), 0.0); }

TEST(HardWall, Branches) {
  const HardWallPotential v;
  const double r = 1.7;
  const double vr = HardWallPotential::gas_phase(r);
  EXPECT_EQ(v.value(0, 0, 7.0, r, 0, 0, 0, 0), vr);   // Z >= zasymp: gas phase
  EXPECT_EQ(v.value(0, 0, 10.0, r, 0, 0, 0, 0), vr);
  EXPECT_EQ(v.value(0, 0, 0.1, r, 0, 0, 0, 0), vr + 1.0 / (0.1 * 0.1 * 0.1 * 0.1 + 0.001));
  const double z = 3.0;
  const double vz = 1.0 / (z * z * z * z + 0.001);
  const double f = 1.0 - (0.5 - 0.5 * std::cos(Constants::pi * (z - 0.1778) / (7.0 - 0.1778)));
  EXPECT_EQ(v.value(1.0, 2.0, z, r, 0.3, 0.4, 5, 6), vr + f * vz);
}

TEST(HardWall, SystemParameters) {
  const auto p = make_potential("hardwallpes");
  EXPECT_EQ(p->system_parameters().z_asymptotic, 7.0);
  EXPECT_EQ(p->system_parameters().lattice_constant_x, 6.23966);
  EXPECT_EQ(p->system_parameters().lattice_constant_y, 5.23966);
  EXPECT_EQ(p->system_parameters().skewing_angle, Constants::pi / 2.0);
  EXPECT_EQ(p->system_parameters().point_group, "C2v");
  EXPECT_THROW((void)make_potential("H2Pt111"), InputError);
}

TEST(PotentialGrid, LayoutAndValues) {
  const ReferenceSetup s(reference_json());
  const HardWallPotential v;
  const PotentialGrid grid = compute_potential_grid(v, s.config, s.grids, s.angular);
  EXPECT_EQ(grid.index.n0, 30u);
  EXPECT_EQ(grid.index.n1, 16u);
  EXPECT_EQ(grid.index.n2, 113u * 56u);
  const double ev2au = units::legacy.ev_to_hartree(1.0);
  // Spot check (itp, ixy, izr) = (iphi=3 + 10*itheta=1, iy=2 + 4*ix=1, ir=7 + 56*iz=50).
  const double expected =
      v.value(s.grids.x[1], s.grids.y[2], s.grids.z_scattering[50], s.grids.r[7], s.angular.thetas()[1],
              s.angular.phis()[3], 50, 7) * ev2au;
  EXPECT_EQ(grid(3 + 10 * 1, 2 + 4 * 1, 7 + 56 * 50), expected);
  // Asymptotic potential: V/chaev (division, not multiplication by 1/chaev).
  EXPECT_EQ(grid.asymptotic[10], units::legacy.ev_to_hartree(HardWallPotential::gas_phase(s.grids.r[10])));
}

TEST(PotentialGrid, IndependentOfThreadCount) {
  const ReferenceSetup s(reference_json());
  const HardWallPotential v;
  const int threads = omp_get_max_threads();
  omp_set_num_threads(1);
  const PotentialGrid a = compute_potential_grid(v, s.config, s.grids, s.angular);
  omp_set_num_threads(4);
  const PotentialGrid b = compute_potential_grid(v, s.config, s.grids, s.angular);
  omp_set_num_threads(threads);
  EXPECT_EQ(a.values, b.values);
}

TEST(PotentialGrid, SwitchingToGasPhase) {
  Json json = reference_json();
  json["UseSwitching"] = true;
  json["ZminSwitch"] = 3.0;
  json["ZmaxSwitch"] = 5.0;
  const ReferenceSetup s(json);
  const HardWallPotential v;
  const PotentialGrid grid = compute_potential_grid(v, s.config, s.grids, s.angular);
  EXPECT_EQ(s.grids.nz_potential, 88);  // first Z >= 5.0 is iz = 88
  const auto block = grid.index.n0 * grid.index.n1;
  const auto& u = units::legacy;
  for (int iz : {40, 60, 75, 87}) {
    const double z = s.grids.z_scattering[iz];
    for (int ir : {0, 9, 30}) {
      const double vr_unswitched = u.ev_to_hartree(HardWallPotential::gas_phase(s.grids.r[ir]));  // V at zasymp = 7
      const double raw = v.value(0, 0, z, s.grids.r[ir], 0, 0, iz, ir) * u.ev_to_hartree(1.0);
      double expected = raw;
      if (z > 3.0 + 1e-12) {
        const double f = 0.5 + 0.5 * std::cos(Constants::pi * (z - 3.0) / 2.0);
        expected = f * (raw - vr_unswitched) + vr_unswitched;
      }
      EXPECT_EQ(grid.values[(static_cast<std::size_t>(ir) + 56u * static_cast<std::size_t>(iz)) * block + 5], expected);
    }
  }
  // vasy is evaluated at the switched zasymp (5.0), where the hard wall is not yet zero.
  EXPECT_EQ(grid.asymptotic[4], u.ev_to_hartree(v.value(0, 0, 5.0, s.grids.r[4], Constants::pi / 2, 0, 88, 4)));
}

TEST(PotentialGrid, ImposedLimits) {
  Json json = reference_json();
  json["ImposeVMin"] = true;
  json["VMin"] = 0.5;
  json["ImposeVMax"] = true;
  json["VMax"] = 3.0;
  const ReferenceSetup s(json);
  const PotentialGrid grid = compute_potential_grid(HardWallPotential(), s.config, s.grids, s.angular);
  const auto& u = units::legacy;
  for (const double x : grid.values) {
    EXPECT_GE(x, u.ev_to_hartree(0.5));
    EXPECT_LE(x, u.ev_to_hartree(3.0));
  }
  EXPECT_EQ(grid.asymptotic[4], u.ev_to_hartree(0.5));  // V(1.4) = 0 clamped
  EXPECT_EQ(grid.asymptotic[0], u.ev_to_hartree(3.0));  // V(0.8) = 4.3 eV clamped
}

TEST(PotentialGrid, PremultiplyWithTimeStep) {
  std::vector<double> v = {1.0, -2.0, 0.5};
  premultiply_with_time_step(v, 2.0);
  EXPECT_EQ(v, (std::vector<double>{-2.0, 4.0, -1.0}));
}
