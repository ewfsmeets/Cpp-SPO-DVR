// Tests: eigensolver, Fourier grid Hamiltonian and gas-phase rovibrational states.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>

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
  // T_ii = (pi^2 / (6 m dx^2)) (1 + 2/N^2); T_ij = (-1)^(i-j) pi^2 / (2 m dx^2 N^2 sin^2(pi (i-j)/N))
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
