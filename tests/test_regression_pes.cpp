// Regression gates 1 and 2 (docs/validation.md) for the corrugated PES cases against the Fortran reference
// files produced by tools/fortran_reference/run_setup_reference.sh <case>:
//   fortran_setup_stdout.txt  vasy table (gate 1)
//   pes_reference.txt         6D subsample of the potential and the angular grid (gate 1)
//   debug_asyp.dat            rovibrational energies and eigenvectors (gate 2)
// Tests are skipped when a reference file is missing.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include <nlohmann/json.hpp>

#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "setup_reference.hpp"

using namespace spodvr;
using spodvr::test::relative_difference;
namespace fs = std::filesystem;

namespace {

// The Fortran grids are evaluated in quad precision and differ from the C++ grids by up to 1 ulp
// (accepted, see docs/fortran_deviations.md). The potential is compared with an absolute floor of 1 eV: near
// V = 0 the LEPS cancellation makes a pure relative test meaningless.
constexpr double grid_tolerance = 4e-16;        // relative, ~2 ulp
constexpr double vasy_tolerance = 1e-14;        // relative to max(|V|, 1 eV)
constexpr double potential_tolerance = 1e-13;   // relative to max(|V|, 1 eV), full 6D potential

double scaled_difference(double a, double b, double floor) {
  return std::abs(a - b) / std::max({std::abs(a), std::abs(b), floor});
}

struct Calculation {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  GasPhase gas;
};

fs::path case_dir(const std::string& name) { return fs::path(SPODVR_REFERENCE_DIR) / name; }

const Calculation& calculation(const std::string& name) {
  static std::map<std::string, std::unique_ptr<Calculation>> cache;
  auto& entry = cache[name];
  if (!entry) {
    std::ifstream file(case_dir(name) / "input.json");
    Config config = build_config(parse_input(nlohmann::json::parse(file)));
    SpatialGrids grids = make_spatial_grids(config);
    AngularGrid angular = make_angular_grid(config);
    const std::unique_ptr<Potential> pes = make_potential(config.input.system_id);
    PotentialGrid potential = compute_potential_grid(*pes, config, grids, angular);
    GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, FftwFft(), LapackEigensolver());
    entry = std::make_unique<Calculation>(
        Calculation{std::move(config), std::move(grids), std::move(angular), std::move(potential), std::move(gas)});
  }
  return *entry;
}

class RegressionPes : public ::testing::TestWithParam<std::string> {
 protected:
  [[nodiscard]] fs::path reference(const std::string& file) const { return case_dir(GetParam()) / file; }
};

#define REQUIRE_REFERENCE(file)                                                                                  \
  if (!fs::exists(reference(file)))                                                                              \
  GTEST_SKIP() << "reference file " << (file) << " missing; run tools/fortran_reference/run_setup_reference.sh " \
               << GetParam()

}  // namespace

// Gate 1a: the Fortran vasy table (r, V in eV).
TEST_P(RegressionPes, Gate1AsymptoticPotentialTable) {
  REQUIRE_REFERENCE("fortran_setup_stdout.txt");
  const auto table = test::read_vasy_table(reference("fortran_setup_stdout.txt"));
  ASSERT_FALSE(table.empty()) << "vasy table not found";
  const Calculation& c = calculation(GetParam());
  ASSERT_EQ(table.size(), c.grids.r.size());
  double max_difference = 0.0;
  for (int ir = 0; ir < c.grids.r.size(); ++ir) {
    const auto [r_ref, v_ref] = table[ir];
    EXPECT_LE(relative_difference(c.grids.r[ir], r_ref), grid_tolerance) << "ir = " << ir;
    const double v = units::legacy.hartree_to_ev(c.potential.asymptotic[ir]);
    const double difference = scaled_difference(v, v_ref, 1.0);
    max_difference = std::max(max_difference, difference);
    EXPECT_LE(difference, vasy_tolerance) << "ir = " << ir << ": " << v << " vs " << v_ref;
  }
  std::cout << "  " << GetParam() << ": vasy max scaled difference " << max_difference << '\n';
}

// Gate 1b: the 6D potential (hartree, before the -dt premultiplication) at the subsampled points
// of the Fortran .pes file, and the angular grid.
TEST_P(RegressionPes, Gate1FullPotentialAndAngularGrid) {
  REQUIRE_REFERENCE("pes_reference.txt");
  const test::PesReference ref = test::read_pes_reference_6d(reference("pes_reference.txt"));
  const Calculation& c = calculation(GetParam());
  ASSERT_EQ(ref.nzpot, c.grids.nz_potential);
  ASSERT_EQ(ref.nr, c.grids.r.size());
  ASSERT_EQ(ref.nx, c.grids.x.size());
  ASSERT_EQ(ref.ny, c.grids.y.size());
  ASSERT_EQ(ref.ntheta, c.angular.n_theta());
  ASSERT_EQ(ref.nphi, c.angular.n_phi());
  for (int t = 0; t < ref.ntheta; ++t) EXPECT_NEAR(c.angular.thetas()[t], ref.thetas[t], 1e-15) << t;
  for (int p = 0; p < ref.nphi; ++p) EXPECT_EQ(c.angular.phis()[p], ref.phis[p]) << p;

  // With surface symmetry the Fortran computed V on its triangle 1 only (C6v, nphi a multiple of
  // 6): the same number of (x, y) points as the C++ representatives (one per orbit).
  if (ref.use_symmetry) {
    EXPECT_EQ(ref.n_xy_computed, c.potential.symmetry.n_representatives());
    EXPECT_FALSE(c.potential.symmetry.is_trivial());
  } else {
    EXPECT_EQ(ref.n_xy_computed, ref.nx * ref.ny);
  }
  ASSERT_GT(ref.points.size(), 10000U);
  const double one_ev = units::legacy.ev_to_hartree(1.0);
  double max_difference = 0.0;
  for (const test::PesReference::Point& p : ref.points) {
    const std::size_t itp = p.iphi + ref.nphi * p.itheta;
    const std::size_t ixy = p.iy + ref.ny * p.ix;
    const std::size_t izr = p.ir + ref.nr * p.iz;
    const double v = c.potential(itp, ixy, izr);
    const double difference = scaled_difference(v, p.v, one_ev);
    max_difference = std::max(max_difference, difference);
    EXPECT_LE(difference, potential_tolerance)
        << "itheta " << p.itheta << " iphi " << p.iphi << " ix " << p.ix << " iy " << p.iy << " iz " << p.iz << " ir "
        << p.ir << ": " << v << " vs " << p.v;
  }
  std::cout << "  " << GetParam() << ": " << ref.points.size() << " points, max scaled difference " << max_difference
            << '\n';
}

// Gate 2: rovibrational energies and eigenvectors vs debug_asyp.dat.
TEST_P(RegressionPes, Gate2RovibrationalStates) {
  REQUIRE_REFERENCE("debug_asyp.dat");
  const auto reference_states = test::read_rovibrational_states(reference("debug_asyp.dat"));
  const Calculation& c = calculation(GetParam());
  const auto& states = c.gas.states;
  ASSERT_EQ(reference_states.size(), (states.n_v() * (states.j_max() + 1)));
  for (const auto& [key, ref] : reference_states) {
    const auto [v, j] = key;
    const double e = states.energy(v, j);
    const auto psi = states.wave_function(v, j);
    ASSERT_EQ(ref.psi.size(), psi.size());
    double max_abs_psi = 0.0;
    for (std::size_t i = 0; i < psi.size(); ++i) max_abs_psi = std::max(max_abs_psi, std::abs(psi[i] - ref.psi[i]));
    if (v < c.gas.n_v_analysis) {
      EXPECT_LE(relative_difference(e, ref.energy), 1e-12) << "v = " << v << ", j = " << j;
      EXPECT_LE(max_abs_psi, 1e-10) << "v = " << v << ", j = " << j;
    }
    // Higher levels: LAPACK vs Jacobi; only a loose check.
    EXPECT_LE(relative_difference(e, ref.energy), 1e-10) << "v = " << v << ", j = " << j;
  }
}

// cu111_dailight_short_sym: C6v surface symmetry (UseSymmetry, nPhisOptimal); the Fortran .pes
// stores the reduced potential and symmetryIndex, expanded by extract_pes.py, and the C++
// potential is gathered through its own maps: the first Fortran gate of the symmetry maps.
INSTANTIATE_TEST_SUITE_P(Cu111, RegressionPes,
                         ::testing::Values("cu111_dailight_small", "cu111_persson_small", "cu111_dailight_short_sym"));
