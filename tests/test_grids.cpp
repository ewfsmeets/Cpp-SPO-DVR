#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <numeric>

#include <nlohmann/json.hpp>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "grids/angular_grid.hpp"
#include "grids/gauss_legendre.hpp"
#include "grids/momentum_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"

using namespace spodvr;

namespace {
Config reference_config() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return build_config(parse_input(nlohmann::json::parse(file)));
}
}  // namespace

TEST(UniformGrid, PointsAsFortran) {
  const UniformGrid g(0.8, 0.15, 56);
  EXPECT_EQ(g.size(), 56);
  EXPECT_EQ(g[4], 0.8 + 4.0 * 0.15);
  EXPECT_EQ(g.back(), 0.8 + 55.0 * 0.15);
}

TEST(MomentumGrid, FftOrderingAndEnergies) {
  const MomentumGrid g(8, 0.5, 2.0);
  const double dk = 2.0 * Constants::pi / (8 * 0.5);
  EXPECT_EQ(g.spacing(), dk);
  const double expected_k[] = {0, dk, 2 * dk, 3 * dk, 4 * dk, -3 * dk, -2 * dk, -dk};
  for (int i = 0; i < 8; ++i) {
    EXPECT_DOUBLE_EQ(g.k()[static_cast<std::size_t>(i)], expected_k[i]);
    EXPECT_DOUBLE_EQ(g.kinetic_energies()[static_cast<std::size_t>(i)], 0.25 * expected_k[i] * expected_k[i]);
  }
  EXPECT_DOUBLE_EQ(g.max_kinetic_energy(), 0.25 * 16 * dk * dk);
}

TEST(GaussLegendre, IntegratesPolynomialsExactly) {
  for (int order = 1; order <= 9; ++order) {
    const auto rule = gauss_associated_legendre(-1.0, 1.0, order, 0);
    for (int p = 0; p <= 2 * order - 1; ++p) {
      double sum = 0.0;
      for (int i = 0; i < order; ++i)
        sum += rule.weights[static_cast<std::size_t>(i)] * std::pow(rule.roots[static_cast<std::size_t>(i)], p);
      const double exact = p % 2 == 1 ? 0.0 : 2.0 / (p + 1);
      EXPECT_NEAR(sum, exact, 1e-14) << "order " << order << " power " << p;
    }
  }
}

TEST(GaussLegendre, FortranIntervalGivesDescendingRoots) {
  const auto rule = gauss_associated_legendre(1.0, -1.0, 5, 0);
  for (std::size_t i = 1; i < 5; ++i) EXPECT_LT(rule.roots[i], rule.roots[i - 1]);
  EXPECT_EQ(rule.roots[2], 0.0);
  EXPECT_NEAR(std::accumulate(rule.weights.begin(), rule.weights.end(), 0.0), -2.0, 1e-15);
}

TEST(AngularGrid, HomonuclearEvenJ) {
  const AngularGrid g(4, 4, true, 0, 0);
  EXPECT_EQ(g.n_theta(), 3);
  EXPECT_EQ(g.n_phi(), 10);
  EXPECT_EQ(g.delta_j(), 2);
  EXPECT_EQ(g.thetas()[2], Constants::pi / 2.0);
  double sum = 0.0;
  for (const double w : g.theta_weights()) sum += w;
  EXPECT_NEAR(sum, 2.0, 1e-15);  // integral of sin(theta) over the full sphere half grid doubled
  EXPECT_EQ(g.phis()[1], Constants::two_pi / 10.0);
}

TEST(AngularGrid, HomonuclearOddJ) {
  const AngularGrid g(5, 5, true, 1, 0);
  EXPECT_EQ(g.j_min(), 1);
  EXPECT_EQ(g.j_max(), 5);
  EXPECT_EQ(g.n_theta(), 3);  // half of 6 roots, none at pi/2
  EXPECT_EQ(g.n_phi(), 12);
  double sum = 0.0;
  for (const double w : g.theta_weights()) sum += w;
  EXPECT_NEAR(sum, 2.0, 1e-15);
}

TEST(AngularGrid, Heteronuclear) {
  const AngularGrid g(3, 2, false, 0, 0);
  EXPECT_EQ(g.n_theta(), 4);
  EXPECT_EQ(g.delta_j(), 1);
  EXPECT_EQ(g.n_phi(), 6);
  double sum = 0.0;
  for (const double w : g.theta_weights()) sum += w;
  EXPECT_NEAR(sum, 2.0, 1e-15);
}

TEST(AngularGrid, InvalidInput) {
  EXPECT_THROW(AngularGrid(4, 4, true, 1, 0), InputError);  // parity
  EXPECT_THROW(AngularGrid(4, 6, true, 0, 0), InputError);  // MJMax > JMax
  EXPECT_THROW(AngularGrid(2, 2, true, 4, 0), InputError);  // JIn > JMax (parity ok)
}

TEST(SpatialGrids, ReferenceCase) {
  const Config config = reference_config();
  const SpatialGrids g = make_spatial_grids(config);
  EXPECT_EQ(g.z_scattering.size(), 176);
  EXPECT_EQ(g.z_specular.size(), 280);
  EXPECT_EQ(g.nz_potential, 113);       // Z = -2 + 0.08*113 = 7.04 is the first gas-phase point
  EXPECT_EQ(g.iz_analysis, 125);
  EXPECT_EQ(g.z_analysis, g.z_scattering[125]);
  EXPECT_EQ(g.z_optical.first, 136);    // Z = 8.88
  EXPECT_EQ(g.z_optical.last, 175);
  ASSERT_TRUE(g.z_specular_optical.has_value());
  EXPECT_EQ(g.z_specular_optical->first, 235);  // Z = 16.8
  EXPECT_EQ(g.z_specular_optical->last, 279);
  EXPECT_EQ(g.r_optical.first, 25);     // r = 4.55
  EXPECT_EQ(g.r_optical.last, 55);
  EXPECT_EQ(g.x.spacing(), 6.239660 / 4.0);
}
