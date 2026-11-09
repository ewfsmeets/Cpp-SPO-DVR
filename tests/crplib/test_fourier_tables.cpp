// Tests: the Fourier coefficient tables (src/crplib/fourier_tables.*): bitwise equal to the tables
// of the Fortran library, and inverse to the Fourier basis at the reference sites (Appendix A of
// docs/crp_library.md), which pins down the site geometry the tables assume.
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "crplib/fourier_tables.hpp"
#include "fortran_reference.hpp"

using namespace crplib::detail;

namespace {

constexpr double two_pi = 2.0 * M_PI;

// Site positions in skewed lattice units (u, v): x = u + v / 2, y = v sqrt(3) / 2.
struct Site {
  double u, v;
};

// The 12 sites of the hexagonal 6D interpolation (Psite index order).
const std::array<Site, 12> molecule_sites = {{{0, 0},
                                              {1. / 2, 0},
                                              {1. / 3, 1. / 3},
                                              {1. / 6, 1. / 6},
                                              {1. / 2, 1. / 2},
                                              {0, 1. / 2},
                                              {2. / 3, 1. / 6},
                                              {1. / 6, 2. / 3},
                                              {2. / 3, 2. / 3},
                                              {5. / 6, 5. / 6},
                                              {1. / 3, 5. / 6},
                                              {5. / 6, 1. / 3}}};

// FOUR6II basis in crystal coordinates (X = 2 pi u, Y = 2 pi v).
std::array<double, 12> basis_6d(Site s) {
  const double x = two_pi * s.u, y = two_pi * s.v;
  return {1.0,
          std::cos(x),
          std::cos(y),
          std::cos(x + y),
          std::sin(x),
          std::sin(y),
          std::sin(x + y),
          std::cos(2 * x + y),
          std::cos(x + 2 * y),
          std::cos(x - y),
          std::cos(2 * x) + std::cos(2 * y) + std::cos(2 * (x + y)),
          -std::sin(2 * x) - std::sin(2 * y) + std::sin(2 * (x + y))};
}

// FOUR3II basis (C6v: the 7 cosine stars; C3v adds 3 sine stars).
std::vector<double> basis_3d(Site s, bool c3v) {
  const double x = two_pi * s.u, y = two_pi * s.v;
  auto c = [](double a) { return std::cos(a); };
  auto sn = [](double a) { return std::sin(a); };
  const double star1 = c(x) + c(y) + c(x + y);
  const double sine1 = -sn(x) - sn(y) + sn(x + y);
  const double star2 = c(2 * x + y) + c(x + 2 * y) + c(x - y);
  const double star3 = c(2 * x) + c(2 * y) + c(2 * (x + y));
  const double sine3 = -sn(2 * x) - sn(2 * y) + sn(2 * (x + y));
  const double star4 = c(3 * x + y) + c(x + 3 * y) + c(2 * x + 3 * y) + c(3 * x + 2 * y) + c(-x + 2 * y) + c(2 * x - y);
  const double sine4 =
      -sn(3 * x + y) - sn(x + 3 * y) + sn(2 * x + 3 * y) + sn(3 * x + 2 * y) - sn(-x + 2 * y) - sn(2 * x - y);
  const double star5 = c(3 * x) + c(3 * y) + c(3 * (x + y));
  const double star6 = c(2 * (2 * x + y)) + c(2 * (x + 2 * y)) + c(2 * (x - y));
  if (!c3v) return {1.0, star1, star2, star3, star4, star5, star6};
  return {1.0, star1, sine1, star2, star3, sine3, star4, sine4, star5, star6};
}

// max |B T^T - I| with B(k, j) = basis_j(site_k) and T(ksit, j) = table[ksit + n j].
template <std::size_t N>
double inverse_error(const std::array<double, N * N>& table, const std::vector<std::vector<double>>& basis) {
  double error = 0.0;
  for (std::size_t k = 0; k < N; ++k) {
    for (std::size_t ksit = 0; ksit < N; ++ksit) {
      double sum = 0.0;
      for (std::size_t j = 0; j < N; ++j) sum += basis[k][j] * table[ksit + N * j];
      error = std::max(error, std::abs(sum - (k == ksit ? 1.0 : 0.0)));
    }
  }
  return error;
}

}  // namespace

TEST(FourierTables, HexagonalMoleculeTableIsTheInverseAtTheTwelveSites) {
  std::vector<std::vector<double>> basis;
  for (const Site& s : molecule_sites) {
    const auto b = basis_6d(s);
    basis.emplace_back(b.begin(), b.end());
  }
  EXPECT_LT(inverse_error<12>(fourier_table_molecule_hexagonal(), basis), 2e-15);
}

TEST(FourierTables, AtomTablesAreTheInverseAtTheirSites) {
  const std::vector<Site> c6v_sites = {{0, 0},           {0, 1. / 2}, {1. / 3, 1. / 3}, {1. / 6, 1. / 6},
                                       {1. / 6, 1. / 3}, {0, 1. / 6}, {0, 2. / 3}};
  const std::vector<Site> c3v_sites = {{0, 0},           {1. / 2, 1. / 2}, {1. / 3, 1. / 3}, {1. / 6, 2. / 3},
                                       {1. / 6, 1. / 2}, {0, 1. / 6},      {0, 1. / 3},      {1. / 3, 5. / 6},
                                       {1. / 2, 2. / 3}, {2. / 3, 2. / 3}};
  std::vector<std::vector<double>> b6, b3;
  for (const Site& s : c6v_sites) b6.push_back(basis_3d(s, false));
  for (const Site& s : c3v_sites) b3.push_back(basis_3d(s, true));
  EXPECT_LT(inverse_error<10>(fourier_table_atom_c3v(), b3), 1e-14);
  // The Fortran DATA of the C6v table has truncated literals (0.666666666666667, 1.333333333333):
  // only ~1e-13 (kept as in the Fortran, docs/fortran_deviations.md).
  const double c6v_error = inverse_error<7>(fourier_table_atom_c6v(), b6);
  EXPECT_LT(c6v_error, 1e-13);
  EXPECT_GT(c6v_error, 1e-14);
}

TEST(FourierTables, BitwiseEqualToTheFortranTables) {
  const auto file = crplib_test::reference_file("kernels/tables.out");
  REQUIRE_CRP_REFERENCE(file);
  // tables.out lines: "<name> <index> <hex> <value>" (written by the C shim, index 0-based).
  const auto records = crplib_test::read_driver_output(file);
  int checked = 0;
  for (const auto& r : records) {
    ASSERT_TRUE(r.has_value) << r.text;
    const auto check = [&](const auto& table) {
      ASSERT_LT(r.index, static_cast<int>(table.size()));
      EXPECT_BITS(table[r.index], r.bits) << " " << r.label << " " << r.index;
      ++checked;
    };
    if (r.label == "fouri6d_c6v" || r.label == "fouri6d_c3v") check(fourier_table_molecule_hexagonal());
    if (r.label == "fouri3d_c6v") check(fourier_table_atom_c6v());
    if (r.label == "fouri3d_c3v") check(fourier_table_atom_c3v());
  }
  EXPECT_EQ(checked, 144 + 144 + 49 + 100);
}
