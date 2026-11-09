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
#include "reference_basis.hpp"

using namespace crplib::detail;

namespace {

using crplib_test::hexagonal_atom_basis;
using crplib_test::hexagonal_molecule_basis;
using crplib_test::hexagonal_molecule_sites;
using crplib_test::Site;

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
  for (const Site& s : hexagonal_molecule_sites) {
    const auto b = hexagonal_molecule_basis(s.u, s.v);
    basis.emplace_back(b.begin(), b.end());
  }
  EXPECT_LT(inverse_error<12>(fourier_table_molecule_hexagonal(), basis), 2e-15);
}

TEST(FourierTables, AtomTablesAreTheInverseAtTheirSites) {
  std::vector<std::vector<double>> b6, b3;
  for (const Site& s : crplib_test::c6v_atom_sites) b6.push_back(hexagonal_atom_basis(s.u, s.v, false));
  for (const Site& s : crplib_test::c3v_atom_sites) b3.push_back(hexagonal_atom_basis(s.u, s.v, true));
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
