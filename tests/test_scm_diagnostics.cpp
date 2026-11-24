// Tests: the SCM set-up diagnostics (scm_diagnostics.*: SCM D3, D10, D15, the slab statistics) on stub
// potentials with known values and on SRP48 with the example slab (also with an atom raised into the
// gas phase); the scans do not depend on the threads. SCM D9: tests/test_scm_potential.cpp.
#include <gtest/gtest.h>

#include <omp.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "core/error.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "potential/scm/scm_diagnostics.hpp"
#include "scm_example_input.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

// The bit pattern of x, for bitwise comparisons (−0 differs from +0, a NaN equals itself).
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

const fs::path& example_directory = test_support::scm_example_directory;

// A PES given by a function of (u, v, Z, r, θ, φ), in eV, on a C3v cell of lattice constant 5 bohr.
class StubPotential final : public Potential {
 public:
  using Function = std::function<double(double, double, double, double, double, double)>;
  explicit StubPotential(Function f) : f_(std::move(f)) {}
  [[nodiscard]] double value(double u, double v, double z, double r, double theta, double phi, int,
                             int) const override {
    return f_(u, v, z, r, theta, phi);
  }
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "stub"; }

 private:
  Function f_;
  SystemParameters parameters_{5.0, 5.0, 10.0, std::numbers::pi / 3.0, "C3v"};
};

// The hard wall reference case: its grids for the stub scans (Z from −2 bohr, r from 0.8 bohr).
struct HardWallSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;

  HardWallSetup()
      : config(build_config(read_input_file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json"))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)) {}
};

const HardWallSetup& hard_wall() {
  static const HardWallSetup setup;
  return setup;
}

// The SCM example input, potential only.
nlohmann::json scm_input(const fs::path& position_file) {
  nlohmann::json json = test_support::scm_example_json(position_file);
  json["Potential"] = true;
  return json;
}

// The example position file with the distorted z of the atom at the origin (ideal (0, 0, 0)) set to
// `z`: an atom that reaches the gas phase at z = 2 bohr (the pair functions reach 10.46 bohr).
fs::path raised_atom_file(const fs::path& dir, double z) {
  std::ifstream in(example_directory / "examplePosfile.dat");
  std::ofstream out(dir / "raised.dat");
  std::string line;
  int raised = 0;
  for (int n = 0; std::getline(in, line); ++n) {
    std::istringstream fields(line);
    double x, y, zd, xi, yi, zi;
    int layer;
    if (n >= 3 && fields >> x >> y >> zd >> xi >> yi >> zi >> layer && xi == 0.0 && yi == 0.0 && zi == 0.0) {
      line = std::format("{} {} {} {} {} {} {}", x, y, z, xi, yi, zi, layer);
      ++raised;
    }
    out << line << '\n';
  }
  EXPECT_EQ(raised, 1);
  return dir / "raised.dat";
}

// The SCM potential, its config and grids for an input.
struct ScmSetup {
  ScmCrpPotential potential;  // built once, as the program does
  Config config;
  SpatialGrids grids;
  AngularGrid angular;

  explicit ScmSetup(const nlohmann::json& json) : ScmSetup(parse_input(json)) {}

 private:
  explicit ScmSetup(const InputParameters& input)
      : potential(input.crp_dataset, input.scm),
        config(build_config(input, potential)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)) {}
};

// The SCM example (SRP48, examplePosfile.dat), set up once for every test that uses it.
const ScmSetup& example() {
  static const ScmSetup setup(scm_input(example_directory / "examplePosfile.dat"));
  return setup;
}

// Expects two points to be bitwise equal in every coordinate.
void expect_same_point(const ScmPoint& a, const ScmPoint& b) {
  EXPECT_EQ(bits(a.u), bits(b.u));
  EXPECT_EQ(bits(a.v), bits(b.v));
  EXPECT_EQ(bits(a.z), bits(b.z));
  EXPECT_EQ(bits(a.r), bits(b.r));
  EXPECT_EQ(bits(a.theta), bits(b.theta));
  EXPECT_EQ(bits(a.phi), bits(b.phi));
}

// A CRP data set's range: Zmin, Rmin, Rmax (bohr).
CrpDataSetSummary crp_range(double z_min, double r_min, double r_max) {
  CrpDataSetSummary range;
  range.z_min = z_min;
  range.r_min = r_min;
  range.r_max = r_max;
  return range;
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Stub potentials with known values
// ---------------------------------------------------------------------------------------------

// V = V_CRP(u/s, v/s) + ΔV with ΔV = 0.3 eV at r = 2 bohr plus a ramp 1e-4·u: at the V_asy points
// (u = 0) only r = 2 bohr changes, counted, reported with its r, and an error naming PosFilepath;
// without ΔV the check passes. The plane holds the largest |ΔV| at the first grid point where it
// occurs. A CRP evaluated at (u, v) instead of (u/s, v/s) gives another plane maximum (the control).
TEST(ScmDiagnostics, ScmD3ChangedAsymptoticPointsAreAnError) {
  const HardWallSetup& s = hard_wall();
  const double scale = 1.25;
  const StubPotential crp(
      [](double u, double v, double, double r, double, double) { return std::sin(u + 2.0 * v) - r; });
  const auto delta = [](double u, double, double, double r, double, double) {
    return (std::abs(r - 2.0) < 1e-9 ? 0.3 : 0.0) + 1e-4 * u;
  };
  const StubPotential scm([&](double u, double v, double z, double r, double t, double p) {
    return crp.value(u / scale, v / scale, z, r, t, p, 0, 0) + delta(u, v, z, r, t, p);
  });
  const ScmAsymptoticCheck check = check_asymptotic_potential(scm, crp, scale, s.config, s.grids, s.angular, true);
  EXPECT_EQ(check.changed_points, 1);
  EXPECT_EQ(check.largest_at_points.points, s.grids.r.size());
  EXPECT_NEAR(check.largest_at_points.value, 0.3, 1e-12);
  EXPECT_NEAR(check.largest_at_points.at.r, 2.0, 1e-9);
  EXPECT_EQ(check.z_asymptotic, s.config.derived.z_asymptotic);

  // The plane: the largest 1e-4·u at the largest x, first in grid order (y, θ, φ first, r below 4.55).
  ASSERT_TRUE(check.largest_on_plane.has_value());
  const ScmExtreme& plane = *check.largest_on_plane;
  EXPECT_NEAR(plane.value, 0.3 + 1e-4 * s.grids.x[s.grids.x.size() - 1], 1e-12);  // r = 2.0 is below ROptStart
  EXPECT_EQ(plane.at.u, s.grids.x[s.grids.x.size() - 1]);
  EXPECT_EQ(plane.at.v, s.grids.y[0]);
  EXPECT_EQ(plane.at.theta, s.angular.thetas()[0]);
  EXPECT_EQ(plane.at.phi, s.angular.phis()[0]);
  int rows = 0;
  while (s.grids.r[rows] < s.config.input.r_opt_start - 1e-9) ++rows;
  EXPECT_EQ(plane.points, 1LL * rows * s.grids.x.size() * s.grids.y.size() * s.angular.n_theta() * s.angular.n_phi());
  EXPECT_TRUE(check.plane_warning());  // 0.3 eV > 1 meV

  ScmSlab slab;
  slab.configuration = 7;
  try {
    require_unchanged_asymptotic_potential(check, slab, ScmSlabStatistics{});
    ADD_FAILURE() << "no error";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "PosFilepath");
    EXPECT_NE(std::string(e.what()).find("configuration 7: ΔV changes the gas-phase potential V_asy at 1 of its"),
              std::string::npos)
        << e.what();
  }

  const ScmAsymptoticCheck unchanged = check_asymptotic_potential(crp, crp, 1.0, s.config, s.grids, s.angular, false);
  EXPECT_EQ(unchanged.changed_points, 0);
  EXPECT_FALSE(unchanged.largest_on_plane.has_value());
  EXPECT_NO_THROW(require_unchanged_asymptotic_potential(unchanged, slab, ScmSlabStatistics{}));
  const ScmAsymptoticCheck unscaled = check_asymptotic_potential(scm, crp, 1.0, s.config, s.grids, s.angular, true);
  EXPECT_GT(std::abs(unscaled.largest_on_plane.value().value - plane.value), 1e-3);  // the control: the CRP at (u, v)
}

// A NaN ΔV at a V_asy point is an error, and the NaN is the reported maximum.
TEST(ScmDiagnostics, ScmD3NanIsAnError) {
  const HardWallSetup& s = hard_wall();
  const StubPotential crp([](double, double, double, double, double, double) { return -1.0; });
  const StubPotential scm([](double, double, double, double r, double, double) {
    return r > 3.0 ? std::numeric_limits<double>::quiet_NaN() : -1.0;
  });
  const ScmAsymptoticCheck check = check_asymptotic_potential(scm, crp, 1.0, s.config, s.grids, s.angular, true);
  EXPECT_GT(check.changed_points, 0);
  EXPECT_TRUE(std::isnan(check.largest_at_points.value));
  EXPECT_TRUE(check.plane_warning());
  EXPECT_THROW(require_unchanged_asymptotic_potential(check, ScmSlab{}, ScmSlabStatistics{}), InputError);
}

// With UseSwitching the switch takes V_gas at the PES's own zasymp (potential_grid.cpp
// switch_to_gas_phase): those points are checked too. Here ZmaxSwitch lies 3 bohr below zasymp, and ΔV
// is −0.1 eV only at zasymp: every r there is counted, and the error names that Z.
TEST(ScmDiagnostics, ScmD3SwitchGasPhasePointsAreChecked) {
  const HardWallSetup& s = hard_wall();
  Config config = s.config;
  config.input.use_switching = true;
  config.derived.z_asymptotic = config.derived.z_asymptotic_unswitched - 3.0;
  const double z_switch = config.derived.z_asymptotic_unswitched;
  const StubPotential crp([](double, double, double, double r, double, double) { return -r; });
  const StubPotential scm(
      [&](double, double, double z, double r, double, double) { return z == z_switch ? -0.1 - r : -r; });
  const ScmAsymptoticCheck check = check_asymptotic_potential(scm, crp, 1.0, config, s.grids, s.angular, false);
  ASSERT_EQ(check.point_heights.size(), 2u);
  EXPECT_EQ(check.point_heights[0], config.derived.z_asymptotic);
  EXPECT_EQ(check.point_heights[1], z_switch);
  EXPECT_EQ(check.changed_points, s.grids.r.size());
  EXPECT_EQ(check.largest_at_points.points, 2 * s.grids.r.size());
  EXPECT_EQ(check.largest_at_points.at.z, z_switch);
  // At z_asymptotic only: the switch's points, 0.1 eV lower, would give −0.1 − r_last.
  EXPECT_EQ(check.gas_phase_minimum, -s.grids.r[s.grids.r.size() - 1]);
  EXPECT_THROW(require_unchanged_asymptotic_potential(check, ScmSlab{}, ScmSlabStatistics{}), InputError);
  config.input.use_switching = false;
  EXPECT_EQ(check_asymptotic_potential(scm, crp, 1.0, config, s.grids, s.angular, false).changed_points, 0);
}

// V = −5 + 0.3·u/a + 0.1·v/a: the jump is 0.3 eV across u = 0 | a and 0.1 eV across v = 0 | a at
// every (Z, r) and φ. With the limit −4.76 eV only the v-edge pairs with t ≤ 5/12 have both sides
// below it (−4.9 + 0.3·t on the far side): 6 of 12 points × 12 φ.
TEST(ScmDiagnostics, ScmD10JumpOfAKnownRamp) {
  const StubPotential scm(
      [](double u, double v, double, double, double, double) { return -5.0 + 0.3 * u / 5.0 + 0.1 * v / 5.0; });
  const ScmEdgeJumps jumps = check_edge_jumps(scm, -4.76);
  EXPECT_EQ(jumps.accessible_limit, -4.76);
  EXPECT_EQ(jumps.points_per_edge, 12);
  EXPECT_EQ(jumps.phi_points, 12);
  ASSERT_EQ(jumps.rows.size(), 6u);
  const double expected_z[] = {2.25, 2.25, 3.0, 3.0, 4.0, 4.0};
  const double expected_r[] = {1.95, 1.4, 1.95, 1.4, 1.95, 1.4};
  for (int k = 0; k < 6; ++k) {
    const ScmEdgeJump& row = jumps.rows[k];
    EXPECT_EQ(row.z, expected_z[k]);
    EXPECT_EQ(row.r, expected_r[k]);
    EXPECT_NEAR(row.max_jump, 0.3, 1e-14);
    EXPECT_EQ(row.accessible_pairs, 6 * 12);
    EXPECT_NEAR(row.max_accessible_jump, 0.1, 1e-14);
  }
  // Nothing accessible: the jump of the accessible pairs stays 0; a NaN V is the maximum jump.
  const ScmEdgeJumps none = check_edge_jumps(scm, -6.0);
  EXPECT_EQ(none.rows[0].accessible_pairs, 0);
  EXPECT_EQ(none.rows[0].max_accessible_jump, 0.0);
  const StubPotential nan([](double u, double, double, double, double, double) {
    return u > 4.0 ? std::numeric_limits<double>::quiet_NaN() : 0.0;
  });
  EXPECT_TRUE(std::isnan(check_edge_jumps(nan, 1.0).rows[0].max_jump));
}

// ΔV = −10·(Zmin − Z)·r·(1 + u) below Zmin = midway between the second and third Z points: two
// planes. The CRP's r range lies midway between r points 1 and 2 (0-based) and 30 and 31, so r points
// 2 … 30 are not clamped: the lowest ΔV is at the first Z, r point 30 and the last x, and at the first
// y, θ and φ (V does not depend on them, so they tie, and ties keep the first point in grid order).
// Above Zmin nothing is scanned.
TEST(ScmDiagnostics, ScmD15LowestBelowZminOfAKnownHole) {
  const HardWallSetup& s = hard_wall();
  const double z_min = 0.5 * (s.grids.z_scattering[1] + s.grids.z_scattering[2]);
  const CrpDataSetSummary range =
      crp_range(z_min, 0.5 * (s.grids.r[1] + s.grids.r[2]), 0.5 * (s.grids.r[30] + s.grids.r[31]));
  const double scale = 0.9;
  const StubPotential crp([](double u, double, double z, double r, double, double) { return std::cos(u) + z * r; });
  const auto hole = [&](double u, double z, double r) { return z < z_min ? -10.0 * (z_min - z) * r * (1.0 + u) : 0.0; };
  const StubPotential scm([&](double u, double v, double z, double r, double t, double p) {
    return crp.value(u / scale, v / scale, z, r, t, p, 0, 0) + hole(u, z, r);
  });
  const ScmHoleCheck check = check_holes_below_z_min(scm, crp, scale, range, s.grids, s.angular);
  EXPECT_EQ(check.planes, 2);
  EXPECT_TRUE(check.warning());
  EXPECT_EQ(check.z_start, s.grids.z_scattering[0]);
  EXPECT_EQ(check.r_min, range.r_min);
  EXPECT_EQ(check.r_max, range.r_max);
  EXPECT_EQ(check.lowest.points,
            2LL * 29 * s.grids.x.size() * s.grids.y.size() * s.angular.n_theta() * s.angular.n_phi());
  const double x_last = s.grids.x[s.grids.x.size() - 1];
  const double r_30 = s.grids.r[30];
  EXPECT_NEAR(check.lowest.value, hole(x_last, s.grids.z_scattering[0], r_30), 1e-11);
  expect_same_point(check.lowest.at,
                    {x_last, s.grids.y[0], s.grids.z_scattering[0], r_30, s.angular.thetas()[0], s.angular.phis()[0]});

  const ScmHoleCheck above = check_holes_below_z_min(
      scm, crp, scale, crp_range(s.grids.z_scattering[0] - 1.0, range.r_min, range.r_max), s.grids, s.angular);
  EXPECT_EQ(above.planes, 0);
  EXPECT_FALSE(above.warning());
  EXPECT_EQ(above.lowest.points, 0);
}

// The two scans split the points at the library's clamp: with Rmin and Rmax on r points 2 and 30
// themselves, the CRP clamps both (r < Rmin + 1e-12, r > Rmax − 1e-12, crplib QdWrapper::clamped_r),
// so the scan below Zmin takes r points 3 … 29 (27) and the scan outside the r range r points 0 … 2
// and 30 … 55 (29) on every Z plane. A Z on Zmin is clamped too (Z < Zmin + 1e-12).
TEST(ScmDiagnostics, ScmD15ScansSplitAtTheLibrarysClamp) {
  const HardWallSetup& s = hard_wall();
  const CrpDataSetSummary range = crp_range(s.grids.z_scattering[1], s.grids.r[2], s.grids.r[30]);
  EXPECT_TRUE(crp_clamps_r(s.grids.r[2], range));
  EXPECT_FALSE(crp_clamps_r(s.grids.r[3], range));
  EXPECT_FALSE(crp_clamps_r(s.grids.r[29], range));
  EXPECT_TRUE(crp_clamps_r(s.grids.r[30], range));
  EXPECT_TRUE(crp_clamps_z(s.grids.z_scattering[1], range));
  EXPECT_FALSE(crp_clamps_z(s.grids.z_scattering[2], range));
  const StubPotential flat([](double, double, double, double, double, double) { return 0.0; });
  const long long per_r = 1LL * s.grids.x.size() * s.grids.y.size() * s.angular.n_theta();
  EXPECT_EQ(check_holes_below_z_min(flat, flat, 1.0, range, s.grids, s.angular).lowest.points,
            2 * 27 * per_r * s.angular.n_phi());
  EXPECT_EQ(check_holes_outside_r_range(flat, flat, 1.0, range, s.grids, s.angular).lowest.points,
            1LL * s.grids.nz_potential * 29 * per_r);
}

// The r region, every Z plane, with the CRP's r range midway between r points 1 and 2 and 20 and 21
// (2 + 35 r points outside it): a hole −(3 − Z)(1 + u)/(1 + (r − r_c)²) where the CRP clamps r is
// lowest at the first Z, r = r_c, the last x and the first y, θ, φ. Its r lies beyond Rmax inside the
// r absorber (r_c = the last r point, 9.05 bohr), beyond Rmax outside it (r_c = 4.1 bohr, below
// ROptStart = 4.55 bohr), or below Rmin (r_c = the first r point, 0.8 bohr). Only the first φ is
// scanned. The scan below Zmin does not see these holes, and this scan does not see a hole below Zmin
// inside the r range: their lowest ΔV is 0.
TEST(ScmDiagnostics, ScmD15LowestOutsideTheRRangeOfAKnownHole) {
  const HardWallSetup& s = hard_wall();
  const double z_min = 0.5 * (s.grids.z_scattering[1] + s.grids.z_scattering[2]);
  const CrpDataSetSummary range =
      crp_range(z_min, 0.5 * (s.grids.r[1] + s.grids.r[2]), 0.5 * (s.grids.r[20] + s.grids.r[21]));
  const StubPotential crp([](double u, double, double z, double r, double, double) { return std::cos(u) + z * r; });
  const auto clamped_hole = [&](double r_c) {
    return StubPotential([&, r_c](double u, double v, double z, double r, double t, double p) {
      const double hole = crp_clamps_r(r, range) ? -(3.0 - z) * (1.0 + u) / (1.0 + (r - r_c) * (r - r_c)) : 0.0;
      return crp.value(u, v, z, r, t, p, 0, 0) + hole;
    });
  };
  const double x_last = s.grids.x[s.grids.x.size() - 1];
  const double z_first = s.grids.z_scattering[0];
  for (const int ir_c : {s.grids.r.size() - 1, 22, 0}) {
    const double r_c = s.grids.r[ir_c];
    SCOPED_TRACE(r_c);
    const ScmOutsideRRangeCheck check =
        check_holes_outside_r_range(clamped_hole(r_c), crp, 1.0, range, s.grids, s.angular);
    EXPECT_EQ(check.r_min, range.r_min);
    EXPECT_EQ(check.r_max, range.r_max);
    EXPECT_EQ(check.theta_points, s.angular.n_theta());
    EXPECT_EQ(check.phi_points, 1);
    EXPECT_EQ(check.lowest.points,
              1LL * s.grids.nz_potential * (2 + 35) * s.grids.x.size() * s.grids.y.size() * s.angular.n_theta());
    EXPECT_NEAR(check.lowest.value, -(3.0 - z_first) * (1.0 + x_last), 1e-12);
    expect_same_point(check.lowest.at,
                      {x_last, s.grids.y[0], z_first, r_c, s.angular.thetas()[0], s.angular.phis()[0]});
    EXPECT_EQ(check.in_r_absorber, ir_c >= s.grids.r_optical.first);
    EXPECT_EQ(check.r_absorber_points, s.grids.r_optical.size());
    EXPECT_EQ(check.r_absorber_first, s.grids.r[s.grids.r_optical.first]);
    // The scan below Zmin covers the r the CRP does not clamp only: it does not see this hole.
    EXPECT_EQ(check_holes_below_z_min(clamped_hole(r_c), crp, 1.0, range, s.grids, s.angular).lowest.value, 0.0);
  }
  EXPECT_FALSE(s.grids.r[22] >= s.config.input.r_opt_start);  // the outside case lies below ROptStart

  // A hole below Zmin inside the r range is not seen by the scan outside it.
  const StubPotential low_hole([&](double u, double v, double z, double r, double t, double p) {
    return crp.value(u, v, z, r, t, p, 0, 0) + (z < z_min && !crp_clamps_r(r, range) ? -50.0 : 0.0);
  });
  EXPECT_EQ(check_holes_outside_r_range(low_hole, crp, 1.0, range, s.grids, s.angular).lowest.value, 0.0);
  EXPECT_EQ(check_holes_below_z_min(low_hole, crp, 1.0, range, s.grids, s.angular).lowest.value, -50.0);
}

// An r absorber that holds no r grid point (ROptEnd below RStart = 0.8 bohr; ROptStart .. ROptEnd
// between the grid points 5.0 and 5.15 bohr) has no points (r_optical has first > last), and no point
// lies inside it.
TEST(ScmDiagnostics, ScmD15EmptyRAbsorber) {
  const StubPotential crp([](double, double, double, double r, double, double) { return -r; });
  const StubPotential scm([](double, double, double, double r, double, double) { return r > 6.0 ? -1.0 - r : -r; });
  for (const auto& [start, end] : {std::pair{0.0, 0.1}, std::pair{5.05, 5.1}}) {
    SCOPED_TRACE(start);
    std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
    nlohmann::json json = nlohmann::json::parse(file);
    json.update({{"ROptStart", start}, {"ROptEnd", end}});
    const Config config = build_config(parse_input(json));
    const SpatialGrids grids = make_spatial_grids(config);
    ASSERT_GT(grids.r_optical.first, grids.r_optical.last);
    const ScmOutsideRRangeCheck check =
        check_holes_outside_r_range(scm, crp, 1.0, crp_range(-10.0, 0.0, 4.0), grids, make_angular_grid(config));
    EXPECT_GT(check.lowest.at.r, 6.0);
    EXPECT_EQ(check.r_absorber_points, 0);
    EXPECT_FALSE(check.in_r_absorber);
  }
}

// The plane and both hole scans give bitwise the same value and point with 1, 3 and 4 threads, also with
// ties (|ΔV| constant over most of the plane) and with a NaN at a later point.
TEST(ScmDiagnostics, ScansDoNotDependOnTheThreads) {
  const HardWallSetup& s = hard_wall();
  const StubPotential crp([](double, double, double, double, double, double) { return 0.0; });
  const StubPotential ties([](double, double, double z, double r, double, double) {
    return (r > 1.5 ? 1.0 : 0.5) * (z < -1.85 ? -1.0 : 1.0);
  });
  const StubPotential nan([](double u, double, double, double r, double, double) {
    return r > 3.0 && u > 1.0 ? std::numeric_limits<double>::quiet_NaN() : r;
  });
  const int threads = omp_get_max_threads();
  const auto compare_threads = [&](const StubPotential& scm) {
    std::vector<ScmExtreme> planes, holes, outside;
    for (const int n : {1, 3, 4}) {
      omp_set_num_threads(n);
      planes.push_back(
          check_asymptotic_potential(scm, crp, 1.0, s.config, s.grids, s.angular, true).largest_on_plane.value());
      holes.push_back(check_holes_below_z_min(scm, crp, 1.0, crp_range(-1.7, 0.0, 100.0), s.grids, s.angular).lowest);
      outside.push_back(
          check_holes_outside_r_range(scm, crp, 1.0, crp_range(-1.7, 0.5, 2.0), s.grids, s.angular).lowest);
    }
    for (int k = 1; k < 3; ++k) {
      EXPECT_EQ(bits(planes[k].value), bits(planes[0].value));
      EXPECT_EQ(planes[k].points, planes[0].points);
      expect_same_point(planes[k].at, planes[0].at);
      EXPECT_EQ(bits(holes[k].value), bits(holes[0].value));
      expect_same_point(holes[k].at, holes[0].at);
      EXPECT_EQ(bits(outside[k].value), bits(outside[0].value));
      EXPECT_EQ(outside[k].points, outside[0].points);
      expect_same_point(outside[k].at, outside[0].at);
    }
    return outside[0];
  };
  (void)compare_threads(ties);
  EXPECT_TRUE(std::isnan(compare_threads(nan).value));  // a NaN outside the r range is the lowest value
  omp_set_num_threads(threads);
  // The ties: the first point of the first row with |ΔV| = 1 (r > 1.5); a hole of −1 from the first Z.
  omp_set_num_threads(4);
  const ScmExtreme plane =
      check_asymptotic_potential(ties, crp, 1.0, s.config, s.grids, s.angular, true).largest_on_plane.value();
  omp_set_num_threads(threads);
  int first_row = 0;
  while (!(s.grids.r[first_row] > 1.5)) ++first_row;
  expect_same_point(plane.at, {s.grids.x[0], s.grids.y[0], s.config.derived.z_asymptotic, s.grids.r[first_row],
                               s.angular.thetas()[0], s.angular.phis()[0]});
}

// rms and largest displacement, highest atom and the mean shift of layer 1, of a slab built in memory.
TEST(ScmDiagnostics, SlabStatistics) {
  ScmSlab slab;
  slab.atoms = {{{0.3, 0.0, 0.4}, {0.0, 0.0, 0.0}, 1},
                {{5.0, 0.1, 0.2}, {5.0, 0.0, 0.0}, 1},
                {{2.5, 1.4, -3.9}, {2.5, 1.4, -3.9}, 2}};
  const ScmSlabStatistics statistics = slab_statistics(slab);
  EXPECT_EQ(statistics.top_layer_atoms, 2);
  EXPECT_EQ(statistics.second_layer_atoms, 1);
  EXPECT_NEAR(statistics.max_displacement, 0.5, 1e-15);
  EXPECT_NEAR(statistics.rms_displacement, std::sqrt((0.25 + 0.05) / 3.0), 1e-15);
  EXPECT_EQ(statistics.highest_atom, 0.4);
  EXPECT_NEAR(statistics.top_layer_shift.x, 0.15, 1e-15);
  EXPECT_NEAR(statistics.top_layer_shift.y, 0.05, 1e-15);
  EXPECT_NEAR(statistics.top_layer_shift.z, 0.3, 1e-15);
}

// ---------------------------------------------------------------------------------------------
// The SRP48 data set with the example slab
// ---------------------------------------------------------------------------------------------

// At z_asymptotic = ZmaxSwitch = 11.9 bohr and at the data set's zasymp = 12.28 bohr, where the switch
// takes V_gas, ΔV leaves the gas-phase potential unchanged at all 2 × 48 points (the highest kept atom
// is at 0.94 bohr, 0.50 bohr beyond the pair functions' reach at 11.9 bohr); on the plane |ΔV| stays
// below 1 meV. No grid plane lies below the CRP's Zmin = 0.47 bohr.
TEST(ScmDiagnostics, ScmD3ExampleLeavesTheAsymptoticPotentialUnchanged) {
  const ScmSetup& s = example();
  const ScmDiagnostics d = diagnose_scm_potential(s.potential, s.config, s.grids, s.angular);
  EXPECT_EQ(d.asymptotic.z_asymptotic, 11.9);
  EXPECT_EQ(d.asymptotic.changed_points, 0);
  ASSERT_EQ(d.asymptotic.point_heights.size(), 2u);
  EXPECT_EQ(d.asymptotic.point_heights[1], s.config.derived.z_asymptotic_unswitched);
  EXPECT_EQ(d.asymptotic.largest_at_points.points, 2 * 48);
  ASSERT_TRUE(d.asymptotic.largest_on_plane.has_value());
  EXPECT_GT(d.asymptotic.largest_on_plane->value, 0.0);
  EXPECT_FALSE(d.asymptotic.plane_warning());
  EXPECT_NO_THROW(require_unchanged_asymptotic_potential(d.asymptotic, s.potential.slab(), d.slab));
  EXPECT_NEAR(d.slab.highest_atom, 0.936367, 1e-6);
  EXPECT_EQ(d.slab.top_layer_atoms + d.slab.second_layer_atoms, std::ssize(s.potential.slab().atoms));
  ASSERT_TRUE(d.holes.has_value());
  EXPECT_EQ(d.holes->planes, 0);
  // Beyond Rmax (4.35 bohr) the lowest V − V_CRP lies in the r absorber (ROptStart = 4.2 bohr), as the
  // second region of scm_open_questions.md Q3 (r = 6–9.8 bohr) does; the point: fortran_deviations.md SCM D15.
  ASSERT_TRUE(d.outside_r_range.has_value());
  EXPECT_GT(d.outside_r_range->lowest.at.r, d.outside_r_range->r_max);
  EXPECT_TRUE(d.outside_r_range->in_r_absorber);
  EXPECT_LT(d.outside_r_range->lowest.value, -100.0);
  ASSERT_TRUE(d.edge_jumps.has_value());
  EXPECT_EQ(d.edge_jumps->rows.size(), 6u);
}

// The example's origin atom raised to z = 2 bohr reaches the V_asy points at Z = 11.9 bohr (9.9 bohr
// below them, the pair functions reach 10.46 bohr): an error naming PosFilepath, before the 6D
// potential, also in an analysis-only run (V_asy is computed there too, SCM D3).
TEST(ScmDiagnostics, ScmD3RaisedAtomIsAnError) {
  const TemporaryDirectory dir("spodvr_scm_d3");
  const nlohmann::json json = scm_input(raised_atom_file(dir.path(), 2.0));
  const ScmSetup s(json);
  const ScmAsymptoticCheck check = check_asymptotic_potential(
      s.potential, s.potential.crp(), *s.potential.rescaling_factor(), s.config, s.grids, s.angular, false);
  EXPECT_GT(check.changed_points, 0);
  EXPECT_GT(check.largest_at_points.value, 1e-6);
  try {
    require_unchanged_asymptotic_potential(check, s.potential.slab(), slab_statistics(s.potential.slab()));
    ADD_FAILURE() << "no error";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "PosFilepath");
    EXPECT_NE(std::string(e.what()).find("the highest lies at z = 2 bohr"), std::string::npos) << e.what();
  }
}

namespace {

// Runs the program on `json` with 4 threads and output into `dir`; returns the exit code, with the
// output in `out` and the errors in `err`.
int run(const nlohmann::json& json, const fs::path& dir, std::string& out, std::string& err) {
  fs::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir / "output";
  cl.threads = 4;
  std::ostringstream out_stream, err_stream;
  const int code = run_program(cl, out_stream, err_stream, false);
  out = out_stream.str();
  err = err_stream.str();
  return code;
}

}  // namespace

// The raised atom stops a potential-only run and an analysis-only run (onlyAnalysis: V_asy is
// computed there too) with exit code 2, after the set-up report and before the 6D potential or the
// BK file is read.
TEST(ScmDiagnostics, ScmD3RaisedAtomStopsEveryRunMode) {
  const TemporaryDirectory dir("spodvr_scm_d3_modes");
  nlohmann::json json = scm_input(raised_atom_file(dir.path(), 2.0));
  json.update({{"nX", 2}, {"nY", 2}});
  std::string out, err;
  EXPECT_EQ(run(json, dir.path() / "potential_only", out, err), 2);
  EXPECT_NE(err.find("PosFilepath: configuration 1: ΔV changes the gas-phase potential V_asy"), std::string::npos)
      << err;
  EXPECT_NE(out.find("SCM DISTORTED SURFACE"), std::string::npos);
  EXPECT_EQ(out.find("6D potential computed"), std::string::npos);

  json.erase("Potential");
  json.update({{"onlyAnalysis", true}, {"RestartDir", (dir.path() / "earlier").string()}});
  fs::create_directories(dir.path() / "earlier");
  std::ofstream(dir.path() / "earlier" / "bk_coefficients.bin") << "not read";
  std::ofstream(dir.path() / "earlier" / "checkpoint.bin") << "not read";  // the specular analysis line
  EXPECT_EQ(run(json, dir.path() / "analysis_only", out, err), 2);
  EXPECT_NE(err.find("ΔV changes the gas-phase potential V_asy"), std::string::npos) << err;
  EXPECT_NE(out.find("Analysis-only run: no 6D potential, so only the gas-phase points are checked."),
            std::string::npos);
}

// SCM D11: the C++ has no .pes cache. A potential-only SCM run writes only the gas-phase files and the
// run metadata (run_info.json) into its output directory, nothing it would read back, and prints the
// SCM block after the CRP data set.
TEST(ScmDiagnostics, ScmD11PotentialOnlyRunWritesNoPesCache) {
  const TemporaryDirectory dir("spodvr_scm_d11");
  nlohmann::json json = scm_input(example_directory / "examplePosfile.dat");
  json.update({{"nX", 2}, {"nY", 2}});
  std::string out, err;
  ASSERT_EQ(run(json, dir.path(), out, err), 0) << err;
  std::vector<std::string> files;
  for (const fs::directory_entry& entry : fs::directory_iterator(dir.path() / "output"))
    files.push_back(entry.path().filename().string());
  std::sort(files.begin(), files.end());
  EXPECT_EQ(files, (std::vector<std::string>{"asymptotic_potential.dat", "gas_phase_states.dat", "run_info.json"}));
  const nlohmann::json info = nlohmann::json::parse(std::ifstream(dir.path() / "output" / "run_info.json"));
  EXPECT_EQ(info["system_id"], "CRPDistrtSurface");
  EXPECT_LT(out.find("CRP DATA SET"), out.find("SCM DISTORTED SURFACE"));
  EXPECT_EQ(out.find("Warning: an ideal atom"), std::string::npos);
}

// SCM D10 on the example near the barrier, (Z, r) = (2.25, 1.95) bohr: at (u, v) = (2a/3, 0) | (2a/3,
// a), φ = 0, the jump is 0.68 eV in docs/scm_open_questions.md (an independent Python transcription
// of ΔV, its figure edge_jump.png). That pair is sampled (t = 8/12, φ = 0), so the maximum is at least
// its jump; the jump where both sides are accessible is at most the unrestricted one. At Z = 4 bohr,
// r = 1.4 bohr the maximum along the edges is 0.13 eV there. Tolerance 0.005 eV: half a unit of the
// last digit printed there (the transcription itself is accurate to well below 1 meV).
TEST(ScmDiagnostics, ScmD10ExampleJumpNearTheBarrier) {
  const ScmSetup& s = example();
  const ScmAsymptoticCheck asymptotic = check_asymptotic_potential(
      s.potential, s.potential.crp(), *s.potential.rescaling_factor(), s.config, s.grids, s.angular, false);
  const ScmEdgeJumps jumps = check_edge_jumps(s.potential, asymptotic.gas_phase_minimum + 1.5);
  const double a = s.potential.system_parameters().lattice_constant_x;
  const double inside = s.potential.value(a * 8 / 12, 0.0, 2.25, 1.95, std::numbers::pi / 2.0, 0.0, 0, 0);
  const double image = s.potential.value(a * 8 / 12, a, 2.25, 1.95, std::numbers::pi / 2.0, 0.0, 0, 0);
  EXPECT_NEAR(std::abs(inside - image), 0.68, 0.005);
  EXPECT_GE(jumps.rows[0].max_jump, std::abs(inside - image));
  EXPECT_GT(jumps.rows[0].accessible_pairs, 0);
  EXPECT_LE(jumps.rows[0].max_accessible_jump, jumps.rows[0].max_jump);
  EXPECT_NEAR(jumps.rows[5].max_jump, 0.13, 0.005);
  // fortran_deviations.md SCM D10 quotes 0.700 eV for this margin: half a unit of its last digit.
  EXPECT_NEAR(jumps.rows[0].max_accessible_jump, 0.700, 0.0005);
}

// SCM D15 on the example: with ZStart = −0.55 bohr the 7 Z planes below Zmin hold V − V_CRP of
// −100 eV and lower inside the CRP's r range (docs/scm_open_questions.md Q3: V_SCM down to −255 eV at
// Z = −0.5 bohr, r = 3.0 bohr); the reported value is V − V_CRP at the reported point. The lowest
// value where the CRP clamps r lies beyond Rmax (Q3's second region). The clamp tests count the
// points the CRP library reports as clamped: 7 Z planes and 30 of the 48 r points.
TEST(ScmDiagnostics, ScmD15ExampleHolesBelowZmin) {
  nlohmann::json json = scm_input(example_directory / "examplePosfile.dat");
  json.update({{"ZStart", -0.55}, {"nZ", 120}, {"nZSpecular", 210}, {"nX", 2}, {"nY", 2}});
  const ScmSetup s(json);
  const double scale = *s.potential.rescaling_factor();
  const CrpDataSetSummary& data_set = s.potential.crp().summary();
  const ScmHoleCheck check =
      check_holes_below_z_min(s.potential, s.potential.crp(), scale, data_set, s.grids, s.angular);
  EXPECT_EQ(check.planes, 7);  // Z = −0.55 … 0.35
  EXPECT_TRUE(check.warning());
  EXPECT_LT(check.lowest.value, -100.0);
  const ScmPoint& p = check.lowest.at;
  const double v = s.potential.value(p.u, p.v, p.z, p.r, p.theta, p.phi, 0, 0);
  const double v_crp = s.potential.crp().value(p.u / scale, p.v / scale, p.z, p.r, p.theta, p.phi, 0, 0);
  EXPECT_EQ(bits(check.lowest.value), bits(v - v_crp));
  EXPECT_FALSE(crp_clamps_r(p.r, data_set));
  const ScmOutsideRRangeCheck outside =
      check_holes_outside_r_range(s.potential, s.potential.crp(), scale, data_set, s.grids, s.angular);
  EXPECT_GT(outside.lowest.at.r, data_set.r_max);
  EXPECT_LT(outside.lowest.value, check.lowest.value);

  const int z_6d = s.grids.nz_potential;
  const CrpRangeReport library =
      s.potential.crp().range_report(s.grids.z_scattering.points().first(z_6d), s.grids.r.points());
  int clamped_z = 0, clamped_r = 0;
  for (int iz = 0; iz < z_6d; ++iz) clamped_z += crp_clamps_z(s.grids.z_scattering[iz], data_set) ? 1 : 0;
  for (int ir = 0; ir < s.grids.r.size(); ++ir) clamped_r += crp_clamps_r(s.grids.r[ir], data_set) ? 1 : 0;
  EXPECT_EQ(clamped_z, library.z_clamped);
  EXPECT_EQ(clamped_r, library.r_clamped);
  EXPECT_EQ(clamped_z, 7);
  EXPECT_EQ(clamped_r, 30);
}
