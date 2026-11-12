// Tests: hard wall, Dai-Light and Persson potentials and the potential grid (layout, switching, limits, threads, -dt).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "comparison.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/dai_light_cu111_potential.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/persson_cu111_potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"

#include <omp.h>

using namespace spodvr;
using spodvr::testing::nan_aware_max;
using Json = nlohmann::json;

namespace {

// The input of a data-free PES, selected by its SystemID.
InputParameters input_with_id(const std::string& id) {
  InputParameters p;
  p.system_id = id;
  return p;
}
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
        angular(make_angular_grid(config)) {}
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
  const auto p = make_potential(input_with_id("hardwallpes"));
  EXPECT_EQ(p->system_parameters().z_asymptotic, 7.0);
  EXPECT_EQ(p->system_parameters().lattice_constant_x, 6.23966);
  EXPECT_EQ(p->system_parameters().lattice_constant_y, 5.23966);
  EXPECT_EQ(p->system_parameters().skewing_angle, Constants::pi / 2.0);
  EXPECT_EQ(p->system_parameters().point_group, "C2v");
  EXPECT_THROW((void)make_potential(input_with_id("H2Pt111")), InputError);
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
  const double expected = v.value(s.grids.x[1], s.grids.y[2], s.grids.z_scattering[50], s.grids.r[7],
                                  s.angular.thetas()[1], s.angular.phis()[3], 50, 7) *
                          ev2au;
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
      EXPECT_EQ(grid.values[(ir + 56u * iz) * block + 5], expected);
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

// --- Corrugated H2/Cu(111) potentials ------------------------------------------------------------

namespace {

constexpr double pi = Constants::pi;

// A molecular configuration (skewed x, y; bohr, rad).
struct Configuration {
  double x, y, z, r, theta, phi;
};

double value_at(const Potential& v, const Configuration& q) {
  return v.value(q.x, q.y, q.z, q.r, q.theta, q.phi, 0, 0);
}

// Configurations near the surface, spread over the cell and over orientations.
std::vector<Configuration> sample_configurations() {
  std::vector<Configuration> q;
  for (int i = 0; i < 7; ++i)
    for (int k = 0; k < 5; ++k)
      q.push_back({0.37 * i + 0.11 * k, 0.53 * k - 0.21 * i, 2.2 + 0.3 * k, 1.2 + 0.1 * i, 0.3 + 0.27 * k + 0.05 * i,
                   0.41 * i + 0.9 * k});
  return q;
}

double relative(double a, double b) { return std::abs(a - b) / std::max({std::abs(a), std::abs(b), 1e-3}); }

}  // namespace

TEST(Cu111Potentials, SystemParametersAndFactory) {
  const auto dai = make_potential(input_with_id(" h2cu111dailight "));
  EXPECT_EQ(dai->name(), "H2Cu111DaiLight");
  EXPECT_EQ(dai->system_parameters().lattice_constant_x, 4.84);
  EXPECT_EQ(dai->system_parameters().lattice_constant_y, 4.84);
  EXPECT_EQ(dai->system_parameters().z_asymptotic, 15.0);
  EXPECT_EQ(dai->system_parameters().skewing_angle, pi / 3.0);
  EXPECT_EQ(dai->system_parameters().point_group, "C6v");

  const auto persson = make_potential(input_with_id("H2CU111PERSSON"));
  EXPECT_EQ(persson->name(), "H2Cu111Persson");
  EXPECT_EQ(persson->system_parameters().lattice_constant_x, 2.5173 / 0.529177);  // cc / cboan
  EXPECT_EQ(persson->system_parameters().lattice_constant_y, 2.5173 / 0.529177);
  EXPECT_EQ(persson->system_parameters().z_asymptotic, 8.0);
  EXPECT_EQ(persson->system_parameters().skewing_angle, pi / 3.0);
  EXPECT_EQ(persson->system_parameters().point_group, "C6v");

  EXPECT_TRUE(potential_is_implemented("H2Cu111DaiLight"));
  EXPECT_TRUE(potential_is_implemented("h2cu111persson"));
  EXPECT_FALSE(potential_is_implemented("H2Cu111"));
}

// Exchanging the two (identical) atoms: (theta, phi) -> (pi - theta, phi + pi).
TEST(Cu111Potentials, HomonuclearExchange) {
  const DaiLightCu111Potential dai;
  const PerssonCu111Potential persson;
  for (const Potential* v : {static_cast<const Potential*>(&dai), static_cast<const Potential*>(&persson)})
    for (const Configuration& q : sample_configurations()) {
      const Configuration swapped{q.x, q.y, q.z, q.r, pi - q.theta, q.phi + pi};
      EXPECT_LE(relative(value_at(*v, q), value_at(*v, swapped)), 1e-12) << v->name();
    }
}

// C6v about the top site at the origin, in skewed coordinates (a1 = (1, 0), a2 = (1/2, sqrt3/2)):
// rotation by 60 deg (x, y) -> (-y, x + y), phi -> phi + pi/3; mirror x_cart -> -x_cart
// (x, y) -> (-x - y, y), phi -> pi - phi; mirror y_cart -> -y_cart (x, y) -> (x + y, -y), phi -> -phi.
TEST(Cu111Potentials, C6vSymmetryAboutTheTopSite) {
  const DaiLightCu111Potential dai;
  const PerssonCu111Potential persson;
  for (const Potential* v : {static_cast<const Potential*>(&dai), static_cast<const Potential*>(&persson)})
    for (const Configuration& q : sample_configurations()) {
      const double v0 = value_at(*v, q);
      const Configuration rotated{-q.y, q.x + q.y, q.z, q.r, q.theta, q.phi + pi / 3.0};
      const Configuration mirror_x{-q.x - q.y, q.y, q.z, q.r, q.theta, pi - q.phi};
      const Configuration mirror_y{q.x + q.y, -q.y, q.z, q.r, q.theta, -q.phi};
      EXPECT_LE(relative(v0, value_at(*v, rotated)), 1e-12) << v->name();
      EXPECT_LE(relative(v0, value_at(*v, mirror_x)), 1e-12) << v->name();
      EXPECT_LE(relative(v0, value_at(*v, mirror_y)), 1e-12) << v->name();
    }
}

TEST(Cu111Potentials, DaiLightIsPeriodicOnTheLattice) {
  const DaiLightCu111Potential dai;
  const double a = dai.system_parameters().lattice_constant_x;
  for (const Configuration& q : sample_configurations()) {
    const double v0 = value_at(dai, q);
    EXPECT_LE(relative(v0, value_at(dai, {q.x + a, q.y, q.z, q.r, q.theta, q.phi})), 1e-12);
    EXPECT_LE(relative(v0, value_at(dai, {q.x, q.y - a, q.z, q.r, q.theta, q.phi})), 1e-12);
  }
}

// Fortran quirk kept: Persson's coordinates use 0.529167 bohr/Angstrom while the lattice constant
// uses cboan = 0.529177, so the potential has period cc/0.529167, not the grid's alat.
TEST(Cu111Potentials, PerssonPeriodIsNotTheGridLattice) {
  const PerssonCu111Potential persson;
  const double period = 2.5173 / 0.529167;
  const double alat = persson.system_parameters().lattice_constant_x;
  EXPECT_NEAR(period / alat - 1.0, 1.89e-5, 0.01e-5);
  double max_period = 0.0;
  double max_lattice = 0.0;
  for (const Configuration& q : sample_configurations()) {
    const double v0 = value_at(persson, q);
    max_period =
        nan_aware_max(max_period, relative(v0, value_at(persson, {q.x + period, q.y, q.z, q.r, q.theta, q.phi})));
    max_lattice =
        nan_aware_max(max_lattice, relative(v0, value_at(persson, {q.x + alat, q.y, q.z, q.r, q.theta, q.phi})));
  }
  EXPECT_LE(max_period, 1e-12);
  EXPECT_GT(max_lattice, 1e-9);
}

// Dai-Light: for Z >= zasymp (tolerance 1e-12) the gas-phase LEPS potential U_M + A_M, independent
// of x, y and the orientation; continuous with the molecule-surface expression far from the surface.
TEST(Cu111Potentials, DaiLightGasPhase) {
  const DaiLightCu111Potential dai;
  for (const double r : {0.9, 1.4, 2.5}) {
    const double v_ref = dai.value(0.0, 0.0, 15.0, r, 0.0, 0.0, 0, 0);
    for (const Configuration& q : sample_configurations()) {
      EXPECT_EQ(dai.value(q.x, q.y, 15.0, r, q.theta, q.phi, 0, 0), v_ref);
      EXPECT_EQ(dai.value(q.x, q.y, 15.0 - 0.5e-12, r, q.theta, q.phi, 0, 0), v_ref);  // .RGE. tolerance
      EXPECT_EQ(dai.value(q.x, q.y, 20.0, r, q.theta, q.phi, 0, 0), v_ref);
    }
    // Just below zasymp the molecule-surface interaction is still ~1e-4..1e-3 eV: the Fortran
    // cut-off at zasymp is a small step in V (switching to the gas phase removes it).
    const double step = std::abs(dai.value(1.0, 2.0, 14.99, r, 0.7, 0.3, 0, 0) - v_ref);
    EXPECT_GT(step, 1e-5);
    EXPECT_LT(step, 1e-3);
  }
  // The gas-phase minimum: U_M + A_M = -D at r = R0, D = 0.1744 hartree.
  EXPECT_NEAR(dai.value(0.0, 0.0, 15.0, 1.40, 0.0, 0.0, 0, 0), units::legacy.hartree_to_ev(-0.1744), 1e-12);
}

// Persson: far from the surface only the molecular LEPS term is left; it is independent of the
// lateral position and the orientation. Its minimum is near -D = -4.505 eV (shifted by the Sato
// parameter and the switching term) at r ~ 1.4 bohr.
TEST(Cu111Potentials, PerssonFarFromTheSurface) {
  const PerssonCu111Potential persson;
  const double z_far = 40.0;
  const double v_ref = persson.value(0.0, 0.0, z_far, 1.4, 0.0, 0.0, 0, 0);
  for (const Configuration& q : sample_configurations())
    EXPECT_NEAR(persson.value(q.x, q.y, z_far, 1.4, q.theta, q.phi, 0, 0), v_ref, 1e-12);
  double v_min = 0.0;
  double r_min = 0.0;
  for (int i = 0; i <= 400; ++i) {
    const double r = 1.0 + 0.002 * i;
    const double v = persson.value(0.0, 0.0, z_far, r, 0.0, 0.0, 0, 0);
    if (v < v_min) {
      v_min = v;
      r_min = r;
    }
  }
  EXPECT_NEAR(v_min, -4.505, 0.03);
  EXPECT_NEAR(r_min, 1.4, 0.05);
}
