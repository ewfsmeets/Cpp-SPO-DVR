// Tests: eigensolver, Fourier grid Hamiltonian and gas-phase rovibrational states.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/fourier_grid_hamiltonian.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/hard_wall_potential.hpp"
#include "potential/potential_grid.hpp"

using namespace spodvr;

TEST(LapackEigensolver, KnownMatrix) {
  // [[2, 1], [1, 2]] -> 1, 3
  const std::vector<double> m = {2.0, 1.0, 1.0, 2.0};
  const auto e = LapackEigensolver().solve(m, 2);
  EXPECT_NEAR(e.values[0], 1.0, 1e-15);
  EXPECT_NEAR(e.values[1], 3.0, 1e-15);
  EXPECT_NEAR(std::abs(e.vector(0)[0]), std::sqrt(0.5), 1e-15);
  EXPECT_NEAR(e.vector(0)[0], -e.vector(0)[1], 1e-15);
}

TEST(FourierGridHamiltonian, KineticMatrixMatchesAnalyticFgh) {
  // Even-N Fourier grid kinetic matrix (Colbert-Miller style, periodic):
  // T_ii = (pi^2 / (6 m dx^2)) (1 + 2/N^2); T_ij = (-1)^(i-j) pi^2 / (m dx^2 N^2 sin^2(pi (i-j)/N))
  // (the Nyquist term k = N/2 is included, as in the Fortran construction via the half spectrum).
  const FftwFft fft;
  const int n = 16;
  const double mass = 3.0;
  const FourierGridHamiltonian fgh(n, 0.0, 4.0, mass, fft);
  const double dx = fgh.spacing();
  const double c = Constants::pi * Constants::pi / (mass * dx * dx);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double expected;
      if (i == j) {
        expected = c / 6.0 * (1.0 + 2.0 / (n * n));
      } else {
        const double s = std::sin(Constants::pi * (i - j) / n);
        expected = ((i - j) % 2 == 0 ? 1.0 : -1.0) * c / (n * n * s * s);
      }
      EXPECT_NEAR(fgh.kinetic()[(i + n * j)], expected, 1e-12 * c) << i << ' ' << j;
    }
  }
}

TEST(FourierGridHamiltonian, HarmonicOscillator) {
  const FftwFft fft;
  const int n = 64;
  const double mass = 1.0;
  const double omega = 1.0;
  const FourierGridHamiltonian fgh(n, -8.0, 8.0, mass, fft);
  std::vector<double> v(n);
  for (int i = 0; i < n; ++i) v[i] = 0.5 * mass * omega * omega * std::pow(fgh.grid()[i], 2);
  const auto e = fgh.eigenstates(v, LapackEigensolver());
  for (int k = 0; k < 6; ++k) EXPECT_NEAR(e.values[k], (k + 0.5) * omega, 1e-10);
  // Normalisation with dx and sign rule (ground state positive everywhere).
  double norm = 0.0;
  for (const double y : e.vector(0)) norm += y * y;
  EXPECT_NEAR(norm * fgh.spacing(), 1.0, 1e-14);
  EXPECT_GT(e.vector(0)[(n / 2)], 0.0);
}

TEST(GasPhase, ReferenceCaseSixDigitValues) {
  // Fortran "INFO FROM SETTING UP THE ASYMPTOTIC ROVIBRATIONAL STATES" (6 printed digits).
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  const Config config = build_config(parse_input(nlohmann::json::parse(file)));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular(4, 4, true, 0, 0);
  const PotentialGrid pot = compute_potential_grid(HardWallPotential(), config, grids, angular);
  const GasPhase gas = compute_gas_phase(config, grids, angular, pot.asymptotic, FftwFft(), LapackEigensolver());
  const auto& u = units::legacy;
  EXPECT_NEAR(u.hartree_to_ev(gas.states.energy(0, 0)), 0.258378, 5e-7);
  EXPECT_NEAR(u.hartree_to_ev(gas.states.energy(1, 0)), 0.762842, 5e-7);
  EXPECT_EQ(gas.n_v_analysis, 2);
  EXPECT_EQ(gas.j_max_analysis, 4);
  for (const auto& s : gas.summary) EXPECT_NEAR(s.norm, 1.0, 1e-12);
}

namespace {

// hardwall_small with the given input changes (validated).
Config hard_wall_config(const nlohmann::json& changes) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json.update(changes);
  return build_config(parse_input(json));
}

// The gas phase of a hard-wall configuration (the hard wall's asymptotic curve).
GasPhase hard_wall_gas_phase(const Config& config) {
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const PotentialGrid pot = compute_potential_grid(HardWallPotential(), config, grids, angular);
  return compute_gas_phase(config, grids, angular, pot.asymptotic, FftwFft(), LapackEigensolver());
}

bool has_warning(const GasPhase& gas, const std::string& text) {
  for (const std::string& w : gas.warnings)
    if (w == text) return true;
  return false;
}

// The bug report's example 1 (docs/fortran_bugs/analysis_limits.md): odd JIn, homonuclear, levels j = 1, 3
// open (E(0, j) - E(0, 0) = 0.01424, 0.08482 eV; available 0.12 + 0.01424), j = 5 closed (0.20923).
const nlohmann::json odd_jin_example = {{"JMax", 5}, {"MJMax", 5},     {"JIn", 1},
                                        {"MJIn", 0}, {"EMinWP", 0.06}, {"EMaxWP", 0.12}};

}  // namespace

// Upstream r1727: the search starts at JIn, so it tests the odd j of the basis (from 0 it gave 2).
TEST(GasPhase, OddJInSearchReachesJ3) {
  const GasPhase gas = hard_wall_gas_phase(hard_wall_config(odd_jin_example));
  EXPECT_EQ(gas.j_max_open, 3);
  EXPECT_EQ(gas.j_max_analysis, 3);
  EXPECT_FALSE(has_warning(gas, "Not all open j channels analysed!"));
}

// JMax 7 and EMaxWP 0.3 eV: j = 5 is open (0.209 <= 0.3 + 0.014), j = 7 closed. The search from 0 gave 6, a
// level of the wrong parity, and the "Rotational basis too small" warning.
TEST(GasPhase, JMax7OpensJ5) {
  nlohmann::json changes = odd_jin_example;
  changes.update({{"JMax", 7}, {"MJMax", 7}, {"EMinWP", 0.2}, {"EMaxWP", 0.3}});
  const GasPhase gas = hard_wall_gas_phase(hard_wall_config(changes));
  EXPECT_EQ(gas.j_max_analysis, 5);
  EXPECT_FALSE(has_warning(gas, "Rotational basis too small to describe all open channels!"));
}

// With JAna 3 the open level j = 5 is not analysed: the warning tests the next level of the basis, j = 5. With
// EMaxWP 0.17 eV j = 4 (0.142 eV) is open but j = 5 (0.209) closed: no warning, so the warning tests j = 5,
// not j = 4.
TEST(GasPhase, OddJInWarnsForTheNextOddLevel) {
  nlohmann::json changes = odd_jin_example;
  changes.update({{"JMax", 7}, {"MJMax", 7}, {"EMinWP", 0.2}, {"EMaxWP", 0.3}, {"JAna", 3}});
  const GasPhase gas = hard_wall_gas_phase(hard_wall_config(changes));
  EXPECT_EQ(gas.j_max_analysis, 3);
  EXPECT_TRUE(has_warning(gas, "Not all open j channels analysed!"));
  changes.update({{"EMinWP", 0.1}, {"EMaxWP", 0.17}});
  const GasPhase even_open = hard_wall_gas_phase(hard_wall_config(changes));
  EXPECT_EQ(even_open.j_max_analysis, 3);
  EXPECT_FALSE(has_warning(even_open, "Not all open j channels analysed!"));
}

// The message of the std::logic_error that f throws ("" if it throws none).
template <class F>
std::string logic_error_message(F f) {
  try {
    f();
  } catch (const std::logic_error& error) {
    return error.what();
  }
  return "";
}

// The post-conditions, reached only through a configuration changed after validation.
TEST(GasPhase, PostConditionsThrow) {
  Config config = hard_wall_config({{"JIn", 2}});
  config.input.j_analysis = 0;
  EXPECT_NE(logic_error_message([&] { (void)hard_wall_gas_phase(config); }).find("is above the highest analysed j"),
            std::string::npos);
  config = hard_wall_config({{"VIn", 1}});
  config.input.v_analysis = 1;
  EXPECT_NE(logic_error_message([&] { (void)hard_wall_gas_phase(config); }).find("is not among the 1 analysed"),
            std::string::npos);
}
