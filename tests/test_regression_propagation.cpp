// Regression gates 3, 4 and 6 (docs/validation.md): the propagation, the on-the-fly analysis and (a case
// with doFluxAnalysis) the flux analysis of tests/reference/hardwall_short (1000 steps) and the
// corrugated cases against the Fortran output of the same case (fortran_stdout.txt, fortran_averages.out,
// fortran_fort21.bin if present, fortran_flux_*, written by
// tools/fortran_reference/run_propagation_reference.sh). Skipped when missing (gate 6: fails).
// The full 7000-step cases (hardwall_small, cu111_dailight_small) run only with
// SPODVR_LONG_TESTS=1. Each case skips while its Fortran reference files are missing.
// LowZGateFailsWithoutTheScatteringAbsorberLong (also SPODVR_LONG_TESTS=1) is the negative test of
// the low-Z case.
#include <algorithm>
#include <cstdlib>
#include <fstream>

#include "regression_gate.hpp"
#include "symmetry/grid_symmetry.hpp"
#include "symmetry/point_group.hpp"
#include "symmetry/surface_lattice.hpp"

using namespace spodvr;
using namespace spodvr::testing;
namespace fs = std::filesystem;

TEST(RegressionPropagation, Gate3ShortRun) { run_gate("hardwall_short"); }

TEST(RegressionPropagation, Gate3FullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("hardwall_small");
}

// The corrugated Dai-Light H2/Cu(111) cases (60 deg cell, C6v): the first gates with theta/phi/x/y
// coupling. _sym: C6v surface symmetry in the potential (UseSymmetry, nPhisOptimal) in both codes.
TEST(RegressionPropagation, Gate3DaiLightShort) { run_gate("cu111_dailight_short"); }
TEST(RegressionPropagation, Gate3DaiLightShortSymmetric) { run_gate("cu111_dailight_short_sym"); }
TEST(RegressionPropagation, Gate3DaiLightFullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("cu111_dailight_small");
}

// The low-Z absorbers of upstream 5.24: cu111_dailight_short with ZOptLow… (energy form) and
// ZSpOptLow… (factor form), both from ZStart to Zana.
TEST(RegressionPropagation, Gate3DaiLightShortLowZ) { run_gate("cu111_dailight_short_lowz"); }

// Negative test of that gate: the C++ run without the scattering low-Z absorber must fail gates 3
// and 4 by at least 100 times their tolerance, so the gate sees that absorber. It is the weaker of
// the two (removing the specular one changes the averages about 100 times more: docs/validation.md).
// Long suite: it cannot catch a code regression the gate misses (a broken absorber fails the gate),
// its outcome changes only with the input or the reference, and a shorter run does not pay: the
// deviation passes 100 tolerances in both gates only from t = 1500 of 2000 (measured 2026-11-12).
TEST(RegressionPropagation, LowZGateFailsWithoutTheScatteringAbsorberLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / "cu111_dailight_short_lowz";
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt", dir / "fortran_averages.out");
  if (!ref || !fs::exists(dir / "fortran_fort21.bin"))
    GTEST_SKIP() << "no Fortran reference in " << dir << "; run "
                 << "tools/fortran_reference/run_propagation_reference.sh cu111_dailight_short_lowz";
  std::ifstream file(dir / "input.json");
  nlohmann::json input = nlohmann::json::parse(file);
  for (const char* key : {"ZOptLowStart", "ZOptLowEnd", "ZOptLowEnergy"}) ASSERT_EQ(input.erase(key), 1u) << key;
  const PropagationRun cpp = run_cpp_propagation(input);
  const auto bk_ref = read_fortran_bk_records(dir / "fortran_fort21.bin");
  ASSERT_EQ(cpp.averages.size(), ref->averages.size());
  ASSERT_EQ(cpp.bk.size(), bk_ref.size());

  // In the units of the gates: printed half units (gate 3 allows 1) and the BK tolerance (gate 4).
  const double averages = averages_deviation(cpp, *ref);
  const BkDeviation bk = bk_deviation(cpp, bk_ref);
  const double bk_in_tolerances = bk.worst / (bk_gate_tolerance * bk.scale);
  EXPECT_GT(averages, 100.0);
  EXPECT_GT(bk_in_tolerances, 100.0);
  std::cout << "  without the scattering low-Z absorber: averages " << averages << " x (half printed unit), BK "
            << bk_in_tolerances << " x the gate-4 tolerance\n";
}

TEST(RegressionPropagation, ParsesTheFortranOutputFormats) {
  // The relevant lines in the formats of WaveFunction.f90 CoherentStateFunc (stdout) and
  // Averages.f90 (averages.out). Five-column lines in the stdout must not count as averages.
  const fs::path stdout_file = fs::temp_directory_path() / "spodvr_test_fortran_stdout.txt";
  const fs::path averages_file = fs::temp_directory_path() / "spodvr_test_fortran_averages.out";
  {
    std::ofstream out(stdout_file);
    out << " * Norm initial wave packet (Z,r) =  0.10000000E+01\n"
        << "   - at k = 0                                              0.48657986E-06\n"
        << "   - at the end/beginning of the grid                      0.79631124E-16\n"
        << "   - maximum value =                                       0.54436230E+00\n"
        << "   - Low boundary (asymptotic value of Z), iLow =                                               126\n"
        << "   - High boundary (start of opt. pot on spec. grid OR a specifically set ZAnaSp ), iHigh =     235\n"
        << " * Norm of the initial wave packet between iLow and iHigh  =  0.10000000E+01\n"
        << "  0    7     0    0  0.4317480E-32\n";
    std::ofstream avg(averages_file);
    avg << "********************************************************************************\n"
        << "               INFO ON PROPAGATING THE WAVE PACKET\n"
        << "********************************************************************************\n\n"
        << "      time             norm           z              r        rest norm\n\n"
        << "  0.00000000E+00  0.10000000E+01  0.11440000E+02  0.14708348E+01  0.10000000E+01\n\n"
        << "  0.10000000E+03  0.99999946E+00  0.11166768E+02  0.14708477E+01  0.99999945E+00\n";
  }
  const auto ref = read_fortran_reference(stdout_file, averages_file);
  fs::remove(stdout_file);
  fs::remove(averages_file);
  ASSERT_TRUE(ref.has_value());
  EXPECT_EQ(ref->norm, 1.0);
  EXPECT_EQ(ref->k_at_zero, 0.48657986E-06);
  EXPECT_EQ(ref->k_maximum, 0.54436230E+00);
  EXPECT_EQ(ref->i_low, 126);
  EXPECT_EQ(ref->i_high, 235);
  ASSERT_EQ(ref->averages.size(), 2u);
  EXPECT_EQ(ref->averages[1].time, 100.0);
  EXPECT_EQ(ref->averages[1].rest_norm, 0.99999945);
}

// The Fortran's own P(x, y; E) of the gate-6 case is C2-symmetric only to about 1e-4 of its largest
// value, not to rounding: the discrete T_xy of the 60-degree grid is not C2-invariant (docs/validation.md,
// gate 6, which the C++ matches; with a C2-invariant T_xy the C++ density is symmetric to rounding:
// FluxPropagation.DensityIsC2InvariantWithAC2SymmetricKineticOperatorLong). This checks that premise
// on the reference itself (measured 2026-11-13: 1.7e-4 for the C++).
TEST(RegressionPropagation, FortranFluxDensityIsC2SymmetricOnlyToTheGridsAccuracy) {
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / "cu111_dailight_short";
  ASSERT_TRUE(fs::exists(dir / "fortran_flux_density.dat"))
      << "run tools/fortran_reference/run_propagation_reference.sh cu111_dailight_short";
  std::ifstream file(dir / "input.json");
  const Config config = build_config(parse_input(nlohmann::json::parse(file)));
  const FluxGrid grid = make_flux_grid(config, make_spatial_grids(config), 1);  // n_rot 1: a block is an ixy
  const int n_energies = flux_energy_count(config.input);
  const GridSymmetry symmetry(SurfaceLattice(config.derived.system), PointGroup(config.derived.system.point_group),
                              grid.nx, grid.ny, make_angular_grid(config).n_phi());
  const auto c2 = std::ranges::find(symmetry.operations(), SurfaceOperation{false, 6});
  ASSERT_NE(c2, symmetry.operations().end());
  const int op = static_cast<int>(c2 - symmetry.operations().begin());

  const std::vector<std::vector<std::string>> lines = read_tokens(dir / "fortran_flux_density.dat");
  auto density = [&](int n, int ix, int iy) {
    return parse_fortran_real(lines.at(fortran_density_line(n, ix, iy, grid, n_energies)).at(3));
  };
  double worst = 0.0;
  double largest = 0.0;
  for (int n = 0; n < n_energies; ++n)
    for (int ix = 0; ix < grid.nx; ++ix)
      for (int iy = 0; iy < grid.ny; ++iy) {
        const FluxGridXY image = flux_block_xy(symmetry.apply_to_xy(op, ix, iy), grid);
        worst = nan_aware_max(worst, std::abs(density(n, image.ix, image.iy) - density(n, ix, iy)));
        largest = nan_aware_max(largest, std::abs(density(n, ix, iy)));
      }
  const double asymmetry = worst / largest;
  EXPECT_GT(asymmetry, 1e-5);
  EXPECT_LT(asymmetry, 1e-3);
  std::cout << "  the Fortran's P(x, y; E): max |P(-x, -y) - P(x, y)| / max |P| = " << asymmetry << '\n';
}
