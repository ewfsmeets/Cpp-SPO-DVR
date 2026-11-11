// Tests: the real CRP data set of the repository, potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v
// (H2/Cu(111) SRP48, C3v, Angstrom and eV; Mondal et al., J. Phys. Chem. A 117, 8770 (2013)):
// reading (incl. the two listed extra cuts without a file), the ranges, the data reproduced at the
// cut knots, the C3v symmetry, the gradient against finite differences, and the speed-ups. Its
// comparison with the Fortran library is the reference case cu111_srp48 (test_atom_surface,
// test_molecule_surface).
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "crplib/dataset.hpp"
#include "crplib/molecule_surface.hpp"
#include "crplib/qd_wrapper.hpp"
#include "crplib/units.hpp"

using namespace crplib;

namespace {

constexpr double degree = std::numbers::pi / 180.0;

const DataSet& real_data_set() {
  static const DataSet data = read_data_set(CRPLIB_REAL_DATA_SET);
  return data;
}

double bohr(double angstrom) { return convert_length(angstrom, LengthUnit::angstrom, LengthUnit::bohr); }
double hartree(double ev) { return convert_energy(ev, EnergyUnit::ev, EnergyUnit::hartree); }

// The atoms of a molecule with centre (x, y, z), bond length r and orientation (theta, phi), as
// the library's reduction of the cuts places them (molecule_surface.cpp).
std::pair<Vec3, Vec3> atoms(double x, double y, double z, double r, double theta, double phi) {
  const Vec3 a{x - 0.5 * r * std::sin(theta) * std::cos(phi), y - 0.5 * r * std::sin(theta) * std::sin(phi),
               z - 0.5 * r * std::cos(theta)};
  const Vec3 b{x + 0.5 * r * std::sin(theta) * std::cos(phi), y + 0.5 * r * std::sin(theta) * std::sin(phi),
               z + 0.5 * r * std::cos(theta)};
  return {a, b};
}

}  // namespace

// The files as gfortran's list-directed input reads them: the trailing comments ("<-- ...",
// "# ...", "ZCM SWITCH OFF MINIMO ...") after the items are skipped.
TEST(RealDataSet, Reading) {
  const DataSet& d = real_data_set();
  EXPECT_EQ(d.symmetry, Symmetry::c3v);
  EXPECT_EQ(d.units, (Units{LengthUnit::angstrom, EnergyUnit::ev}));
  EXPECT_EQ(d.lattice_constant, 2.601516559);
  EXPECT_EQ(d.interlayer_spacings, (std::vector<double>{2.09713, 2.10385, 2.09713}));
  EXPECT_EQ(d.displacement_1, 0.0);
  ASSERT_TRUE(d.info.has_value());
  EXPECT_NE(d.info->short_info.find("H2/Cu(111)"), std::string::npos);
  EXPECT_EQ(d.atom_cuts.size(), 10U);
  EXPECT_EQ(d.corrugation.x.back(), 5.6);
  EXPECT_FALSE(d.cuts_interpolation);
  EXPECT_EQ(d.zso_min, 3.4);  // "3.4d0"
  EXPECT_EQ(d.zso_max, 4.0);
  EXPECT_EQ(d.z_dependence.x.back(), 6.5);

  ASSERT_EQ(d.molecule_cuts.size(), 31U);
  for (int k = 0; k < 29; ++k) EXPECT_TRUE(d.molecule_cuts[k].present) << k;
  EXPECT_FALSE(d.molecule_cuts[29].present);  // Cu111_t2fAT90P30.2d: listed, no file
  EXPECT_FALSE(d.molecule_cuts[30].present);  // Cu111_t2fBT90P90.2d
  const MoleculeCut& top = d.molecule_cuts[0];
  EXPECT_EQ(top.file_name, "Cu111_topT0.2d");
  ASSERT_EQ(top.r.size(), 14U);
  ASSERT_EQ(top.z.size(), 15U);
  EXPECT_EQ(top.r.front(), 0.4);
  EXPECT_EQ(top.r.back(), 2.3);
  EXPECT_EQ(top.z.front(), 0.25);
  EXPECT_EQ(top.z.back(), 4.0);
  EXPECT_EQ(top.v.front(), 1950.0);  // "0.250000  0.400000   1950.0   #3018.9095661653"
  EXPECT_EQ(top.v.back(), 4.8212461653);

  const std::vector<std::string> warnings = data_set_warnings(d);
  ASSERT_EQ(warnings.size(), 2U);
  EXPECT_NE(warnings[1].find("Cu111_t2fAT90P30.2d, Cu111_t2fBT90P90.2d"), std::string::npos) << warnings[1];
}

// Rmin, Rmax, Zmin: the intersection of the 29 cut grids (the 2 cuts without a file do not count);
// the atom cuts reach every atom of an evaluated molecule; the corrugation function stays inside
// the C3v cluster's coverage (no warning beyond the skipped cuts).
TEST(RealDataSet, Ranges) {
  const DataSet& d = real_data_set();
  const MoleculeSurface pes(d);
  double r_min = 0.0, r_max = 100.0, z_min = -100.0;
  for (const MoleculeCut& cut : d.molecule_cuts) {
    if (!cut.present) continue;
    r_min = std::max(r_min, cut.r.front());
    r_max = std::min(r_max, cut.r.back());
    z_min = std::max(z_min, cut.z.front());
  }
  EXPECT_EQ(r_min, 0.4);
  EXPECT_EQ(r_max, 2.3);
  EXPECT_EQ(z_min, 0.25);
  EXPECT_EQ(pes.r_min(), bohr(0.4));
  EXPECT_EQ(pes.r_max(), bohr(2.3));
  EXPECT_EQ(pes.z_min(), bohr(0.25));
  EXPECT_EQ(pes.zso_min(), bohr(3.4));
  EXPECT_EQ(pes.zso_max(), bohr(4.0));
  EXPECT_EQ(pes.z_asymptotic(), bohr(6.5));
  EXPECT_EQ(pes.warnings().size(), 2U);
  EXPECT_EQ(pes.cut_grids(), 1);  // S4: one shared grid

  const std::vector<double>& atom_z = d.atom_cuts.front().potential.x;  // Angstrom
  EXPECT_LE(atom_z.front(), 0.25 - 0.5 * 2.3);
  EXPECT_GE(atom_z.back(), 4.0 + 0.5 * 2.3);

  const ClusterSize cluster = pes.atom_surface().cluster_size();
  EXPECT_EQ(cluster.atoms, 256);
  std::cout << "  cluster atoms within reach: " << cluster.within_reach << " of " << cluster.atoms << '\n';
}

// The data are total energies (CutsInterpolation = F): at a cut's site and orientation, at its knots
// below ZsoMin, V = data. Exact up to the rounding of the reduction and of the interpolation. Knots
// on the grid edges are left out: the centre recomputed from the atoms (CAMBCORD) can fall a
// rounding outside the grid there, which is the silent V = 0 branch.
TEST(RealDataSet, ReproducesTheDataAtTheCutKnots) {
  const DataSet& d = real_data_set();
  const MoleculeSurface pes(d);
  const double a = bohr(d.lattice_constant);
  double worst = 0.0;
  int n = 0;
  for (int k = 0; k < 29; ++k) {
    const MoleculeCut& cut = d.molecule_cuts[k];
    const double ucm = cut.x * a, vcm = cut.y * a;  // skewed lattice units, as the reduction
    const double x = ucm + vcm / 2.0, y = std::sqrt(3.0) * vcm / 2.0;
    for (std::size_t j = 1; j < cut.z.size(); ++j) {
      if (!(cut.z[j] < d.zso_min)) continue;
      for (std::size_t i = 1; i + 1 < cut.r.size(); ++i) {
        const auto [atom_a, atom_b] = atoms(x, y, bohr(cut.z[j]), bohr(cut.r[i]), cut.theta * degree, cut.phi * degree);
        const double data = hartree(cut.v[i + cut.r.size() * j]);
        const std::optional<MoleculeSurfaceValue> value = pes.evaluate(atom_a, atom_b);
        ASSERT_TRUE(value.has_value()) << cut.file_name << " at Z " << cut.z[j] << ", r " << cut.r[i];
        const double v = value->v;
        const double difference = std::abs(v - data) / std::max(std::abs(data), hartree(1.0));
        worst = std::max(worst, difference);
        EXPECT_LT(difference, 1e-11) << cut.file_name << " at Z " << cut.z[j] << ", r " << cut.r[i] << ": " << v
                                     << " vs " << data;
        ++n;
      }
    }
  }
  EXPECT_GT(n, 2000);
  std::cout << "  " << n << " knots, max |V - data| / max(|data|, 1 eV) = " << worst << '\n';
}

// V is invariant under the C3v operations about the top site (rotations by 120 degrees, the mirror
// x -> -x) and under the exchange of the atoms, to rounding.
TEST(RealDataSet, SymmetricToRounding) {
  const MoleculeSurface pes(real_data_set());
  const double c = std::cos(2.0 * std::numbers::pi / 3.0), s = std::sin(2.0 * std::numbers::pi / 3.0);
  auto rotate = [&](const Vec3& p) { return Vec3{c * p.x - s * p.y, s * p.x + c * p.y, p.z}; };
  auto mirror = [](const Vec3& p) { return Vec3{-p.x, p.y, p.z}; };
  double worst = 0.0;
  for (double x : {0.3, 1.7, 3.1})
    for (double y : {0.2, 1.1, 2.6})
      for (double z : {1.2, 2.5, 4.0, 6.9, 8.0})
        for (double r : {1.0, 1.4, 2.4, 3.6})
          for (const auto& [theta, phi] : {std::pair{0.0, 0.0}, {0.7, 0.4}, {1.5707963267948966, 2.2}, {2.4, -1.1}}) {
            const auto [a0, b0] = atoms(x, y, z, r, theta, phi);
            const auto v0 = pes.evaluate(a0, b0);
            ASSERT_TRUE(v0.has_value());
            std::vector<std::pair<Vec3, Vec3>> images = {{b0, a0},
                                                         {rotate(a0), rotate(b0)},
                                                         {rotate(rotate(a0)), rotate(rotate(b0))},
                                                         {mirror(a0), mirror(b0)},
                                                         {mirror(rotate(a0)), mirror(rotate(b0))},
                                                         {mirror(rotate(rotate(a0))), mirror(rotate(rotate(b0)))}};
            for (const auto& [a1, b1] : images) {
              const double difference = std::abs(pes.evaluate(a1, b1)->v - v0->v) / std::max(std::abs(v0->v), 1e-3);
              worst = std::max(worst, difference);
            }
          }
  EXPECT_LT(worst, 1e-11);
  std::cout << "  max |V(g q) - V(q)| / max(|V|, 1e-3 Ha) = " << worst << '\n';
}

// The analytic gradient against central differences of V, below ZsoMin, in the switching region
// and above ZsoMax (where V is the gas-phase potential).
TEST(RealDataSet, GradientMatchesFiniteDifferences) {
  const MoleculeSurface pes(real_data_set());
  const double h = 1e-5;
  for (double z : {2.1, 4.4, 6.9, 8.2})
    for (const auto& [x, y] : {std::pair{0.37, 0.61}, {2.9, 1.3}})
      for (const auto& [theta, phi] : {std::pair{0.4, 0.3}, {1.3, 2.1}}) {
        const auto [a, b] = atoms(x, y, z, 1.45, theta, phi);
        const MoleculeSurfaceValue v = *pes.evaluate(a, b);
        const std::array<double, 6> analytic = {v.gradient_a.x, v.gradient_a.y, v.gradient_a.z,
                                                v.gradient_b.x, v.gradient_b.y, v.gradient_b.z};
        for (int i = 0; i < 6; ++i) {
          Vec3 ap = a, am = a, bp = b, bm = b;
          double& plus = i < 3 ? (i == 0 ? ap.x : i == 1 ? ap.y : ap.z) : (i == 3 ? bp.x : i == 4 ? bp.y : bp.z);
          double& minus = i < 3 ? (i == 0 ? am.x : i == 1 ? am.y : am.z) : (i == 3 ? bm.x : i == 4 ? bm.y : bm.z);
          plus += h;
          minus -= h;
          const double numeric = (pes.evaluate(ap, bp)->v - pes.evaluate(am, bm)->v) / (2.0 * h);
          EXPECT_NEAR(analytic[i], numeric, 1e-6 * std::max(1.0, std::abs(numeric)))
              << "component " << i << " at Z " << z << ", theta " << theta;
        }
      }
}

// Every combination of the speed-ups (S2, S3, S4) gives bitwise the plain QD-wrapper value on a
// dynamics-like grid (also where r and Z are clamped and above Zasy).
TEST(RealDataSet, SpeedUpsAreBitwise) {
  const DataSet& d = real_data_set();
  const QdWrapper plain(d, EvaluationOptions{.prune_cluster = false, .share_grid_lookup = false, .value_only = false});
  const QdWrapper fast(d);
  const double a = plain.parameters().lattice_constant;
  int differences = 0, n = 0;
  for (int iu = 0; iu < 6; ++iu)
    for (int iv = 0; iv < 6; ++iv)
      for (double z = -1.0; z < 15.0; z += 0.7)
        for (double r : {0.4, 0.9, 1.4, 2.2, 3.3, 4.6})
          for (double theta : {0.0, 0.8, std::numbers::pi / 2})
            for (double phi : {-1.0, 0.6, 2.4}) {
              const double u = iu * a / 6.0, v = iv * a / 6.0;
              const double p = plain.potential_ev(u, v, z, r, theta, phi),
                           f = fast.potential_ev(u, v, z, r, theta, phi);
              if (std::bit_cast<std::uint64_t>(p) != std::bit_cast<std::uint64_t>(f)) ++differences;
              ++n;
            }
  EXPECT_EQ(differences, 0) << " of " << n;
}
