// Regression gates 1 and 2 (docs/validation.md) for tests/reference/hardwall_small against the Fortran
// reference files produced by tools/fortran_reference/run_setup_reference.sh:
//   fortran_setup_stdout.txt  vasy table (gate 1), list-directed full precision
//   pes_reference.txt         the full potential V(iz, ir) and the angular grid (gate 1)
//   debug_asyp.dat            rovibrational energies and eigenvectors (gate 2)
// Tests are skipped when a reference file is missing.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "comparison.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/potential_grid.hpp"
#include "setup_reference.hpp"

using namespace spodvr;
using spodvr::test::parse_fortran_real;
using spodvr::test::relative_difference;
using spodvr::testing::nan_aware_max;
namespace fs = std::filesystem;

namespace {

// The Fortran grids (Setup.f90 `rmin + dfloat(i-1)*delr`) are evaluated in quad precision and
// differ from the C++ double grids by up to 1 ulp. Accepted as within numerical limits; the
// tolerances below allow for it and its propagation through the steep hard wall potential.
constexpr double grid_tolerance = 4e-16;       // relative, ~2 ulp
constexpr double potential_tolerance = 1e-13;  // relative, full 6D potential

const fs::path reference_dir = fs::path(SPODVR_REFERENCE_DIR) / "hardwall_small";

struct Calculation {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  GasPhase gas;
};

const Calculation& calculation() {
  static const Calculation calc = [] {
    std::ifstream file(reference_dir / "input.json");
    Config config = build_config(parse_input(nlohmann::json::parse(file)));
    SpatialGrids grids = make_spatial_grids(config);
    AngularGrid angular = make_angular_grid(config);
    PotentialGrid potential = compute_potential_grid(HardWallPotential(), config, grids, angular);
    GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, FftwFft(), LapackEigensolver());
    return Calculation{std::move(config), std::move(grids), std::move(angular), std::move(potential), std::move(gas)};
  }();
  return calc;
}

#define REQUIRE_REFERENCE(file)            \
  if (!fs::exists(reference_dir / (file))) \
  GTEST_SKIP() << "reference file " << (file) << " missing; run tools/fortran_reference/run_setup_reference.sh"

}  // namespace

// Gate 1a: the Fortran "INFO FROM SETTING UP THE ROVIBARIONAL ASYMPTOTIC POTENTIAL" table,
// which prints r and vasy*chaev.
TEST(RegressionHardWall, Gate1AsymptoticPotentialTable) {
  REQUIRE_REFERENCE("fortran_setup_stdout.txt");
  const auto table = test::read_vasy_table(reference_dir / "fortran_setup_stdout.txt");
  ASSERT_FALSE(table.empty()) << "vasy table not found";
  const Calculation& c = calculation();
  ASSERT_EQ(table.size(), (c.grids.r.size()));
  const auto& u = units::legacy;
  double max_rel = 0.0;
  for (int ir = 0; ir < c.grids.r.size(); ++ir) {
    const auto [r_ref, v_ref] = table[ir];
    EXPECT_LE(relative_difference(c.grids.r[ir], r_ref), grid_tolerance) << "ir = " << ir;
    const double v = u.hartree_to_ev(c.potential.asymptotic[ir]);
    max_rel = nan_aware_max(max_rel, relative_difference(v, v_ref));
    EXPECT_LE(relative_difference(v, v_ref), 1e-14) << "ir = " << ir << ": " << v << " vs " << v_ref;
  }
  RecordProperty("max_relative_difference", std::to_string(max_rel));
  EXPECT_EQ(table[4].second, 0.0);  // V(r = 1.4) = 0 exactly
}

// Gate 1b: the complete 6D potential (hartree, before the -dt premultiplication) and the
// angular grid, extracted from the Fortran .pes file.
TEST(RegressionHardWall, Gate1FullPotentialAndAngularGrid) {
  REQUIRE_REFERENCE("pes_reference.txt");
  std::ifstream file(reference_dir / "pes_reference.txt");
  std::vector<std::string> lines;
  for (std::string line; std::getline(file, line);)
    if (!line.empty() && line[0] != '#') lines.push_back(line);
  std::istringstream header(lines.at(0));
  int nzpot = 0, nr = 0, nx = 0, ny = 0, ntheta = 0, nphi = 0;
  header >> nzpot >> nr >> nx >> ny >> ntheta >> nphi;

  const Calculation& c = calculation();
  ASSERT_EQ(nzpot, c.grids.nz_potential);
  ASSERT_EQ(nr, c.grids.r.size());
  ASSERT_EQ(nx, c.grids.x.size());
  ASSERT_EQ(ny, c.grids.y.size());
  ASSERT_EQ(ntheta, c.angular.n_theta());
  ASSERT_EQ(nphi, c.angular.n_phi());

  std::size_t pos = 1;
  for (int t = 0; t < ntheta; ++t, ++pos)
    EXPECT_NEAR(c.angular.thetas()[t], parse_fortran_real(lines.at(pos)), 1e-15) << t;
  for (int p = 0; p < nphi; ++p, ++pos) EXPECT_EQ(c.angular.phis()[p], parse_fortran_real(lines.at(pos))) << p;

  double max_rel = 0.0;
  const std::size_t block = c.potential.index.n0 * c.potential.index.n1;
  for (int izr = 0; izr < nzpot * nr; ++izr, ++pos) {
    std::istringstream in(lines.at(pos));
    int iz = 0, ir = 0;
    std::string v_text;
    in >> iz >> ir >> v_text;
    ASSERT_EQ((iz - 1) * nr + (ir - 1), izr);
    const double v_ref = parse_fortran_real(v_text);
    for (std::size_t k = 0; k < block; ++k) {
      const double v = c.potential.values[izr * block + k];
      max_rel = nan_aware_max(max_rel, relative_difference(v, v_ref));
    }
  }
  RecordProperty("max_relative_difference", std::to_string(max_rel));
  EXPECT_LE(max_rel, potential_tolerance);
}

// Gate 2: rovibrational energies and eigenvectors vs debug_asyp.dat (DoDump).
TEST(RegressionHardWall, Gate2RovibrationalStates) {
  REQUIRE_REFERENCE("debug_asyp.dat");
  const auto reference = test::read_rovibrational_states(reference_dir / "debug_asyp.dat");
  const Calculation& c = calculation();
  const auto& states = c.gas.states;
  ASSERT_EQ(reference.size(), (states.n_v() * (states.j_max() + 1)));

  // The analysed (open) levels are the regression gate; all levels are reported.
  double max_rel_energy_open = 0.0;
  double max_abs_psi_open = 0.0;
  double max_rel_energy_all = 0.0;
  for (const auto& [key, ref] : reference) {
    const auto [v, j] = key;
    const double e = states.energy(v, j);
    const auto psi = states.wave_function(v, j);
    ASSERT_EQ(ref.psi.size(), psi.size());
    double max_abs_psi = 0.0;
    for (std::size_t i = 0; i < psi.size(); ++i)
      max_abs_psi = nan_aware_max(max_abs_psi, std::abs(psi[i] - ref.psi[i]));
    max_rel_energy_all = nan_aware_max(max_rel_energy_all, relative_difference(e, ref.energy));
    if (v < c.gas.n_v_analysis) {
      max_rel_energy_open = nan_aware_max(max_rel_energy_open, relative_difference(e, ref.energy));
      max_abs_psi_open = nan_aware_max(max_abs_psi_open, max_abs_psi);
      EXPECT_LE(relative_difference(e, ref.energy), 1e-12) << "v = " << v << ", j = " << j;
      EXPECT_LE(max_abs_psi, 1e-10) << "v = " << v << ", j = " << j;
    }
    // Higher levels: LAPACK vs Jacobi; only a loose check.
    EXPECT_LE(relative_difference(e, ref.energy), 1e-10) << "v = " << v << ", j = " << j;
  }
  RecordProperty("max_relative_energy_difference_open", std::to_string(max_rel_energy_open));
  RecordProperty("max_abs_psi_difference_open", std::to_string(max_abs_psi_open));
  RecordProperty("max_relative_energy_difference_all", std::to_string(max_rel_energy_all));
  std::cout << "  open levels: max rel dE = " << max_rel_energy_open << ", max |dpsi| = " << max_abs_psi_open
            << "; all levels: max rel dE = " << max_rel_energy_all << '\n';
}
