// Tests: surface lattice (reciprocal metric, diffraction order, lattice action), surface
// operations and point groups, and the point group on the (x, y, phi) grid (orbits, maps).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fft_sizes.hpp"
#include "symmetry/grid_symmetry.hpp"
#include "symmetry/point_group.hpp"
#include "symmetry/surface_lattice.hpp"

using namespace spodvr;

namespace {

// Test-local port of the Fortran table, the reference for the lattice-based order.
// Fortran: RotDifBasis.f90 getDiffractionOrder
int fortran_diffraction_order(const std::string& point_group, int n, int m) {
  const int abs_n = std::abs(n);
  const int abs_m = std::abs(m);
  if (point_group == "C2v" || point_group == "C4v") return abs_n + abs_m;
  if (n * m > 0) return std::max(abs_n, abs_m);  // C3, C3v, C6v
  return abs_n + abs_m;
}

struct LatticeCase {
  std::string name;
  std::string point_group;  // the Fortran point group of such a cell
  SurfaceLattice lattice;
};

// The cells of the Fortran systems (hard wall; Cu(100)-like square; the two Cu(111) PESs).
std::vector<LatticeCase> fortran_lattices() {
  const double pi = Constants::pi;
  const double persson = units::legacy.angstrom_to_bohr(2.5173);
  return {
      {"hard wall (rectangular)", "C2v", SurfaceLattice(6.239660, 5.239660, pi / 2.0)},
      {"rectangular, a < b", "C2v", SurfaceLattice(4.0, 7.5, pi / 2.0)},
      {"square", "C4v", SurfaceLattice(4.82, 4.82, pi / 2.0)},
      {"Dai-Light Cu(111)", "C6v", SurfaceLattice(4.84, 4.84, pi / 3.0)},
      {"Persson Cu(111)", "C6v", SurfaceLattice(persson, persson, pi / 3.0)},
      {"hexagonal C3v", "C3v", SurfaceLattice(5.0, 5.0, pi / 3.0)},
  };
}

}  // namespace

TEST(SurfaceLattice, Accessors) {
  const SurfaceLattice lattice(4.0, 5.0, 1.2);
  EXPECT_EQ(lattice.lattice_constant_x(), 4.0);
  EXPECT_EQ(lattice.lattice_constant_y(), 5.0);
  EXPECT_EQ(lattice.skewing_angle(), 1.2);
}

TEST(SurfaceLattice, FromSystemParameters) {
  const SystemParameters system{6.0, 7.0, 8.0, Constants::pi / 3.0, "C6v"};
  const SurfaceLattice lattice(system);
  EXPECT_EQ(lattice.lattice_constant_x(), 6.0);
  EXPECT_EQ(lattice.lattice_constant_y(), 7.0);
  EXPECT_EQ(lattice.skewing_angle(), Constants::pi / 3.0);
}

// |b1|^2 = 1/(a sin)^2, |b2|^2 = 1/(b sin)^2 and b1.b2 = -cos/(a b sin^2) (a_i . b_j = delta_ij).
TEST(SurfaceLattice, ReciprocalMetric) {
  const double a = 4.0;
  const double b = 5.5;
  const double skew = 1.1;
  const SurfaceLattice lattice(a, b, skew);
  const double s = std::sin(skew);
  const double c = std::cos(skew);
  EXPECT_NEAR(lattice.reciprocal_length_squared(1, 0), 1.0 / (a * a * s * s), 1e-15);
  EXPECT_NEAR(lattice.reciprocal_length_squared(0, 1), 1.0 / (b * b * s * s), 1e-15);
  const double b1_dot_b2 = -c / (a * b * s * s);
  for (int n = -3; n <= 3; ++n)
    for (int m = -3; m <= 3; ++m) {
      const double expected = n * n * lattice.reciprocal_length_squared(1, 0) +
                              m * m * lattice.reciprocal_length_squared(0, 1) + 2.0 * n * m * b1_dot_b2;
      EXPECT_NEAR(lattice.reciprocal_length_squared(n, m), expected, 1e-14) << n << " " << m;
    }
  EXPECT_EQ(lattice.reciprocal_length_squared(0, 0), 0.0);
}

// Same metric as the parallel kinetic energy ekinxy = (kx^2 + ky^2 - 2 cos kx ky)/2M with
// kx = n 2 pi/(a sin), ky = m 2 pi/(b sin) (propagation_tables.cpp, Fortran DifChanEnergies).
TEST(SurfaceLattice, MetricMatchesParallelKineticEnergy) {
  const double a = 4.84;
  const double b = 4.84;
  const double skew = Constants::pi / 3.0;
  const SurfaceLattice lattice(a, b, skew);
  const double dkx = Constants::two_pi / (a * std::sin(skew));
  const double dky = Constants::two_pi / (b * std::sin(skew));
  for (int n = -4; n <= 4; ++n)
    for (int m = -4; m <= 4; ++m) {
      const double kx = n * dkx;
      const double ky = m * dky;
      const double k2 = kx * kx + ky * ky - 2.0 * std::cos(skew) * kx * ky;
      const double expected = lattice.reciprocal_length_squared(n, m) * Constants::two_pi * Constants::two_pi;
      EXPECT_NEAR(k2, expected, 1e-12 * std::max(1.0, expected)) << n << " " << m;
    }
}

TEST(SurfaceLattice, ShortDiagonalSteps) {
  const double pi = Constants::pi;
  // Rectangular and square: neither diagonal is short.
  for (const SurfaceLattice& l : {SurfaceLattice(6.2, 5.2, pi / 2.0), SurfaceLattice(4.0, 4.0, pi / 2.0)}) {
    EXPECT_FALSE(l.has_step_sum());
    EXPECT_FALSE(l.has_step_difference());
  }
  // Hexagonal described with a 60 deg cell: b1 + b2 is as short as b1 (exact tie).
  const SurfaceLattice hex60(4.84, 4.84, pi / 3.0);
  EXPECT_TRUE(hex60.has_step_sum());
  EXPECT_FALSE(hex60.has_step_difference());
  // The same lattice described with a 120 deg cell: now b1 - b2 is the short diagonal.
  const SurfaceLattice hex120(4.84, 4.84, 2.0 * pi / 3.0);
  EXPECT_FALSE(hex120.has_step_sum());
  EXPECT_TRUE(hex120.has_step_difference());
  // An oblique cell (75 deg): no diagonal is short.
  const SurfaceLattice oblique(4.0, 4.0, 75.0 * pi / 180.0);
  EXPECT_FALSE(oblique.has_step_sum());
  EXPECT_FALSE(oblique.has_step_difference());
}

// The lattice-based order equals the Fortran per-point-group table for every channel of grids
// up to 64 x 64 on the cells of the Fortran systems.
TEST(SurfaceLattice, DiffractionOrderEqualsFortranTable) {
  for (const LatticeCase& c : fortran_lattices()) {
    for (int n = -31; n <= 32; ++n)
      for (int m = -31; m <= 32; ++m)
        ASSERT_EQ(c.lattice.diffraction_order(n, m), fortran_diffraction_order(c.point_group, n, m))
            << c.name << ": n = " << n << ", m = " << m;
  }
}

TEST(SurfaceLattice, DiffractionOrderExamples) {
  const double pi = Constants::pi;
  const SurfaceLattice rectangular(6.2, 5.2, pi / 2.0);
  EXPECT_EQ(rectangular.diffraction_order(0, 0), 0);
  EXPECT_EQ(rectangular.diffraction_order(2, -1), 3);
  EXPECT_EQ(rectangular.diffraction_order(-2, -2), 4);
  const SurfaceLattice hex(4.84, 4.84, pi / 3.0);
  EXPECT_EQ(hex.diffraction_order(2, 1), 2);
  EXPECT_EQ(hex.diffraction_order(-1, -3), 3);
  EXPECT_EQ(hex.diffraction_order(2, -1), 3);
  EXPECT_EQ(hex.diffraction_order(0, -2), 2);
  // 120 deg description of the same lattice: the max rule applies where n m < 0.
  const SurfaceLattice hex120(4.84, 4.84, 2.0 * pi / 3.0);
  EXPECT_EQ(hex120.diffraction_order(2, -1), 2);
  EXPECT_EQ(hex120.diffraction_order(2, 1), 3);
}

// The order is a step distance: 0 only at (0, 0), symmetric under (n, m) -> (-n, -m), and it
// changes by at most one between channels one shortest step apart (triangle inequality).
TEST(SurfaceLattice, DiffractionOrderIsAStepDistance) {
  const double pi = Constants::pi;
  for (const SurfaceLattice& l : {SurfaceLattice(6.2, 5.2, pi / 2.0), SurfaceLattice(4.84, 4.84, pi / 3.0),
                                  SurfaceLattice(4.84, 4.84, 2.0 * pi / 3.0)}) {
    std::vector<std::pair<int, int>> steps = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    if (l.has_step_sum()) steps.insert(steps.end(), {{1, 1}, {-1, -1}});
    if (l.has_step_difference()) steps.insert(steps.end(), {{1, -1}, {-1, 1}});
    for (int n = -10; n <= 10; ++n)
      for (int m = -10; m <= 10; ++m) {
        const int order = l.diffraction_order(n, m);
        EXPECT_EQ(order == 0, n == 0 && m == 0);
        EXPECT_EQ(order, l.diffraction_order(-n, -m));
        bool has_closer_neighbour = order == 0;
        for (const auto& [dn, dm] : steps) {
          const int neighbour = l.diffraction_order(n + dn, m + dm);
          EXPECT_LE(std::abs(neighbour - order), 1) << n << " " << m;
          if (neighbour == order - 1) has_closer_neighbour = true;
        }
        EXPECT_TRUE(has_closer_neighbour) << n << " " << m;
      }
  }
}

// Orders are shells, not symmetry orbits: on a square lattice (2, 0) and (1, 1) are not
// equivalent but both have order 2 (as in the Fortran, whose sums per order rely on it).
TEST(SurfaceLattice, OrdersAreShellsNotOrbits) {
  const SurfaceLattice square(4.0, 4.0, Constants::pi / 2.0);
  EXPECT_EQ(square.diffraction_order(2, 0), 2);
  EXPECT_EQ(square.diffraction_order(1, 1), 2);
  EXPECT_NE(square.reciprocal_length_squared(2, 0), square.reciprocal_length_squared(1, 1));
}

// =================================================================================================
// Surface operations
// =================================================================================================

namespace {

// All 24 operations (rotations and mirrors in steps of 30 deg).
std::vector<SurfaceOperation> all_operations() {
  std::vector<SurfaceOperation> ops;
  for (const bool mirror : {false, true})
    for (int k = 0; k < 12; ++k) ops.push_back({mirror, k});
  return ops;
}

// Angle difference reduced to (-pi, pi].
double angle_difference(double a, double b) {
  const double d = std::remainder(a - b, Constants::two_pi);
  return std::abs(d);
}

}  // namespace

TEST(SurfaceOperation, ComposeMatchesActionOnAzimuthAndVectors) {
  const std::vector<double> phis = {0.0, 0.3, 1.7, 3.1, 4.4, 6.0};
  const CartesianVector v{0.7, -1.3};
  for (const SurfaceOperation& a : all_operations())
    for (const SurfaceOperation& b : all_operations()) {
      const SurfaceOperation ab = compose(a, b);
      EXPECT_GE(ab.angle_twelfths, 0);
      EXPECT_LT(ab.angle_twelfths, 12);
      for (const double phi : phis)
        EXPECT_LT(angle_difference(apply_to_azimuth(ab, phi), apply_to_azimuth(a, apply_to_azimuth(b, phi))), 1e-12);
      const CartesianVector direct = apply_to_cartesian(ab, v);
      const CartesianVector nested = apply_to_cartesian(a, apply_to_cartesian(b, v));
      EXPECT_NEAR(direct.x, nested.x, 1e-14);
      EXPECT_NEAR(direct.y, nested.y, 1e-14);
    }
}

TEST(SurfaceOperation, GroupAxioms) {
  const std::vector<SurfaceOperation> ops = all_operations();
  for (const SurfaceOperation& g : ops) {
    EXPECT_EQ(compose(g, identity_operation()), g);
    EXPECT_EQ(compose(identity_operation(), g), g);
    EXPECT_EQ(compose(g, inverse(g)), identity_operation());
    EXPECT_EQ(compose(inverse(g), g), identity_operation());
    if (g.mirror) {
      EXPECT_EQ(compose(g, g), identity_operation());  // reflections are involutions
    }
  }
  for (const SurfaceOperation& a : ops)
    for (const SurfaceOperation& b : ops)
      for (const SurfaceOperation& c : ops) ASSERT_EQ(compose(compose(a, b), c), compose(a, compose(b, c)));
}

// The dihedral group D12 in closed form: r^i r^j = r^(i+j), r^i m_j = m_(i+j), m_i r^j = m_(i-j),
// m_i m_j = r^(i-j) (indices mod 12).
TEST(SurfaceOperation, MultiplicationTableIsDihedral) {
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j) {
      const int sum = (i + j) % 12;
      const int difference = ((i - j) % 12 + 12) % 12;
      EXPECT_EQ(compose({false, i}, {false, j}), (SurfaceOperation{false, sum}));
      EXPECT_EQ(compose({false, i}, {true, j}), (SurfaceOperation{true, sum}));
      EXPECT_EQ(compose({true, i}, {false, j}), (SurfaceOperation{true, difference}));
      EXPECT_EQ(compose({true, i}, {true, j}), (SurfaceOperation{false, difference}));
    }
}

// Mirror (true, k) reflects in the line at 15 k degrees: vectors on that line are fixed.
TEST(SurfaceOperation, MirrorLine) {
  for (int k = 0; k < 12; ++k) {
    const double line = k * Constants::pi / 12.0;
    const CartesianVector on_line{std::cos(line), std::sin(line)};
    const CartesianVector image = apply_to_cartesian({true, k}, on_line);
    EXPECT_NEAR(image.x, on_line.x, 1e-15);
    EXPECT_NEAR(image.y, on_line.y, 1e-15);
    const CartesianVector normal{-std::sin(line), std::cos(line)};
    const CartesianVector flipped = apply_to_cartesian({true, k}, normal);
    EXPECT_NEAR(flipped.x, -normal.x, 1e-15);
    EXPECT_NEAR(flipped.y, -normal.y, 1e-15);
  }
}

TEST(SurfaceOperation, Describe) {
  EXPECT_EQ(describe(identity_operation()), "identity");
  EXPECT_EQ(describe({false, 4}), "rotation 120 deg");
  EXPECT_EQ(describe({true, 6}), "mirror at 90 deg");
  EXPECT_EQ(describe({true, 0}), "mirror at 0 deg");
}

// =================================================================================================
// Point groups
// =================================================================================================

TEST(PointGroup, OrdersAndElements) {
  const std::map<std::string, std::vector<SurfaceOperation>> expected = {
      {"C2v", {{false, 0}, {false, 6}, {true, 0}, {true, 6}}},
      {"C4v", {{false, 0}, {false, 3}, {false, 6}, {false, 9}, {true, 0}, {true, 3}, {true, 6}, {true, 9}}},
      {"C6v",
       {{false, 0},
        {false, 2},
        {false, 4},
        {false, 6},
        {false, 8},
        {false, 10},
        {true, 0},
        {true, 2},
        {true, 4},
        {true, 6},
        {true, 8},
        {true, 10}}},
      {"C3v", {{false, 0}, {false, 4}, {false, 8}, {true, 2}, {true, 6}, {true, 10}}},
      {"C3", {{false, 0}, {false, 4}, {false, 8}}},
      {"none", {{false, 0}}},
  };
  for (const auto& [name, elements] : expected) {
    const PointGroup group(name);
    EXPECT_EQ(group.name(), name);
    EXPECT_EQ(group.operations(), elements) << name;
    EXPECT_EQ(group.order(), static_cast<int>(elements.size())) << name;
    EXPECT_TRUE(is_group(group.operations())) << name;
    EXPECT_EQ(group_closure(group.generators()), group.operations()) << name;
  }
}

// C3v about a top site: the mirrors through the top site that keep both hollows (1/3, 1/3) and
// (2/3, 2/3) of the 60 deg cell in place (lines at 30, 90, 150 deg); those at 0, 60, 120 deg would
// exchange the fcc and hcp hollows.
TEST(PointGroup, C3vMirrorsKeepTheHollows) {
  const SurfaceLattice hex(1.0, 1.0, Constants::pi / 3.0);
  const PointGroup c3v("C3v");
  for (const SurfaceOperation& g : c3v.operations())
    for (const double hollow : {1.0 / 3.0, 2.0 / 3.0}) {
      const auto [x, y] = hex.to_skewed(apply_to_cartesian(g, hex.to_cartesian(hollow, hollow)));
      EXPECT_NEAR(std::remainder(x - hollow, 1.0), 0.0, 1e-12) << describe(g);
      EXPECT_NEAR(std::remainder(y - hollow, 1.0), 0.0, 1e-12) << describe(g);
    }
  for (const int k : {0, 4, 8}) {
    const auto [x, y] = hex.to_skewed(apply_to_cartesian({true, k}, hex.to_cartesian(1.0 / 3.0, 1.0 / 3.0)));
    EXPECT_NEAR(std::remainder(x - 2.0 / 3.0, 1.0), 0.0, 1e-12) << k;
    EXPECT_NEAR(std::remainder(y - 2.0 / 3.0, 1.0), 0.0, 1e-12) << k;
  }
}

TEST(PointGroup, NamesAndErrors) {
  EXPECT_EQ(PointGroup(" c6V ").name(), "C6v");
  EXPECT_EQ(PointGroup("NONE").name(), "none");
  EXPECT_THROW(PointGroup("D4h"), InputError);
  EXPECT_THROW(PointGroup(""), InputError);
}

// Fortran getCxvNumPhisSmallestMultiple: 2, 4, 6 for C2v, C4v, C6v.
TEST(PointGroup, PhiPointsMultiple) {
  EXPECT_EQ(phi_points_multiple(PointGroup("C2v").operations()), 2);
  EXPECT_EQ(phi_points_multiple(PointGroup("C4v").operations()), 4);
  EXPECT_EQ(phi_points_multiple(PointGroup("C6v").operations()), 6);
  EXPECT_EQ(phi_points_multiple(PointGroup("C3v").operations()), 6);
  EXPECT_EQ(phi_points_multiple(PointGroup("C3").operations()), 3);
  EXPECT_EQ(phi_points_multiple(PointGroup("none").operations()), 1);
}

TEST(PointGroup, IsGroupDetectsNonGroups) {
  EXPECT_FALSE(is_group({{false, 4}}));                         // no identity
  EXPECT_FALSE(is_group({{false, 0}, {false, 4}}));             // not closed
  EXPECT_FALSE(is_group({{false, 0}, {true, 0}, {true, 2}}));  // m0 m2 is a rotation
  EXPECT_TRUE(is_group({{false, 0}, {true, 5}}));
}

// =================================================================================================
// Lattice action
// =================================================================================================

namespace {

struct GroupCase {
  std::string group;
  SurfaceLattice lattice;
};

std::vector<GroupCase> group_cases() {
  const double pi = Constants::pi;
  return {
      {"C2v", SurfaceLattice(6.2, 5.2, pi / 2.0)},   {"C4v", SurfaceLattice(4.82, 4.82, pi / 2.0)},
      {"C6v", SurfaceLattice(4.84, 4.84, pi / 3.0)}, {"C3v", SurfaceLattice(5.0, 5.0, pi / 3.0)},
      {"C3", SurfaceLattice(5.0, 5.0, pi / 3.0)},    {"none", SurfaceLattice(4.0, 5.0, 1.2)},
  };
}

IntegerMatrix2 multiply(const IntegerMatrix2& a, const IntegerMatrix2& b) {
  return {a.a11 * b.a11 + a.a12 * b.a21, a.a11 * b.a12 + a.a12 * b.a22, a.a21 * b.a11 + a.a22 * b.a21,
          a.a21 * b.a12 + a.a22 * b.a22};
}

}  // namespace

TEST(LatticeAction, MatricesOfTheGroups) {
  for (const GroupCase& c : group_cases()) {
    const PointGroup group(c.group);
    const SurfaceLattice& l = c.lattice;
    // Gram matrix of (a1, a2) in fractional coordinates.
    const double a = l.lattice_constant_x();
    const double b = l.lattice_constant_y();
    const double g11 = a * a;
    const double g12 = a * b * std::cos(l.skewing_angle());
    const double g22 = b * b;
    for (const SurfaceOperation& g : group.operations()) {
      const std::optional<IntegerMatrix2> m = l.lattice_matrix(g);
      ASSERT_TRUE(m.has_value()) << c.group << " " << describe(g);
      EXPECT_EQ(m->a11 * m->a22 - m->a12 * m->a21, g.mirror ? -1 : 1) << c.group << " " << describe(g);
      // The metric is preserved: A^T G A = G.
      const double t11 = m->a11 * (g11 * m->a11 + g12 * m->a21) + m->a21 * (g12 * m->a11 + g22 * m->a21);
      const double t12 = m->a11 * (g11 * m->a12 + g12 * m->a22) + m->a21 * (g12 * m->a12 + g22 * m->a22);
      const double t22 = m->a12 * (g11 * m->a12 + g12 * m->a22) + m->a22 * (g12 * m->a12 + g22 * m->a22);
      EXPECT_NEAR(t11, g11, 1e-12 * g11);
      EXPECT_NEAR(t12, g12, 1e-12 * g11);
      EXPECT_NEAR(t22, g22, 1e-12 * g22);
    }
    // Homomorphism: A(g h) = A(g) A(h).
    for (const SurfaceOperation& g : group.operations())
      for (const SurfaceOperation& h : group.operations())
        EXPECT_EQ(*l.lattice_matrix(compose(g, h)), multiply(*l.lattice_matrix(g), *l.lattice_matrix(h)));
  }
}

TEST(LatticeAction, IncompatibleLattices) {
  const double pi = Constants::pi;
  const SurfaceLattice rectangular(6.2, 5.2, pi / 2.0);
  const SurfaceLattice square(4.0, 4.0, pi / 2.0);
  const SurfaceLattice hexagonal(4.0, 4.0, pi / 3.0);
  EXPECT_FALSE(rectangular.lattice_matrix({false, 3}).has_value());  // 90 deg needs a = b
  EXPECT_FALSE(rectangular.lattice_matrix({true, 3}).has_value());   // diagonal mirror
  EXPECT_FALSE(square.lattice_matrix({false, 2}).has_value());       // 60 deg
  EXPECT_FALSE(hexagonal.lattice_matrix({false, 3}).has_value());    // 90 deg
  EXPECT_TRUE(hexagonal.lattice_matrix({true, 0}).has_value());      // C2v is a subgroup of C6v
  // Examples: the 60 deg rotation of the hexagonal cell (x, y) -> (-y, x + y).
  EXPECT_EQ(*hexagonal.lattice_matrix({false, 2}), (IntegerMatrix2{0, -1, 1, 1}));
  EXPECT_EQ(*hexagonal.lattice_matrix({true, 6}), (IntegerMatrix2{-1, -1, 0, 1}));  // x' = -x - y
  EXPECT_EQ(*square.lattice_matrix({true, 3}), (IntegerMatrix2{0, 1, 1, 0}));       // x <-> y
  EXPECT_THROW(GridSymmetry(rectangular, PointGroup("C4v"), 8, 8, 8), InputError);
  EXPECT_THROW(GridSymmetry(square, PointGroup("C6v"), 8, 8, 12), InputError);
  EXPECT_THROW(GridSymmetry(hexagonal, PointGroup("C6v"), 1, 8, 12), InputError);
  EXPECT_THROW(GridSymmetry(hexagonal, PointGroup("C6v"), 8, 1, 12), InputError);
}

TEST(LatticeAction, CartesianRoundTrip) {
  const SurfaceLattice l(4.0, 5.5, 1.1);
  for (const double x : {0.0, 1.3, -2.1})
    for (const double y : {0.0, 0.7, 4.9}) {
      const auto [x2, y2] = l.to_skewed(l.to_cartesian(x, y));
      EXPECT_NEAR(x2, x, 1e-14);
      EXPECT_NEAR(y2, y, 1e-14);
    }
}

// The reciprocal lattice transforms with A^-T: symmetry-equivalent channels have the same
// length (parallel energy at normal incidence) and the same diffraction order.
TEST(LatticeAction, DiffractionChannelOrbits) {
  for (const GroupCase& c : group_cases()) {
    const PointGroup group(c.group);
    for (const SurfaceOperation& g : group.operations()) {
      const IntegerMatrix2 a = *c.lattice.lattice_matrix(g);
      const int det = a.a11 * a.a22 - a.a12 * a.a21;
      for (int n = -6; n <= 6; ++n)
        for (int m = -6; m <= 6; ++m) {
          // (n', m') = A^-T (n, m); A^-1 = [[a22, -a12], [-a21, a11]] / det.
          const int n2 = (a.a22 * n - a.a21 * m) / det;
          const int m2 = (-a.a12 * n + a.a11 * m) / det;
          EXPECT_EQ(c.lattice.diffraction_order(n2, m2), c.lattice.diffraction_order(n, m))
              << c.group << " " << describe(g) << " " << n << " " << m;
          const double l0 = c.lattice.reciprocal_length_squared(n, m);
          EXPECT_NEAR(c.lattice.reciprocal_length_squared(n2, m2), l0, 1e-12 * std::max(1.0, l0));
        }
    }
  }
}

// =================================================================================================
// The point group on the (x, y, phi) grid
// =================================================================================================

namespace {

const std::vector<int> grid_sizes = {2, 4, 6, 8, 10, 12, 16, 18, 20, 24};
const std::vector<int> phi_sizes = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 30};

// All (group, nx, ny, nphi) combinations of the tests: nx = ny for the hexagonal and square
// cells and a few nx != ny ones (where operations mixing x and y must be dropped).
struct GridCase {
  GroupCase group;
  int nx, ny, nphi;
};

std::vector<GridCase> grid_cases() {
  std::vector<GridCase> cases;
  for (const GroupCase& c : group_cases())
    for (const int n : grid_sizes)
      for (const int nphi : phi_sizes) {
        cases.push_back({c, n, n, nphi});
        if (n == 8) {
          cases.push_back({c, 8, 6, nphi});
          cases.push_back({c, 4, 12, nphi});
        }
      }
  return cases;
}

int index_of(const std::vector<SurfaceOperation>& ops, const SurfaceOperation& g) {
  for (std::size_t i = 0; i < ops.size(); ++i)
    if (ops[i] == g) return static_cast<int>(i);
  return -1;
}

// Fractional coordinate difference reduced to (-1/2, 1/2].
double wrapped(double d) { return std::remainder(d, 1.0); }

}  // namespace

TEST(GridSymmetry, UsableOperationsFormASubgroup) {
  for (const GridCase& c : grid_cases()) {
    const PointGroup group(c.group.group);
    const GridSymmetry s(c.group.lattice, group, c.nx, c.ny, c.nphi);
    EXPECT_TRUE(is_group(s.operations())) << c.group.group << " " << c.nx << "x" << c.ny << " nphi " << c.nphi;
    EXPECT_EQ(s.operations().size() + s.dropped_operations().size(), group.operations().size());
    EXPECT_EQ(s.group_order(), group.order());
    EXPECT_EQ(s.point_group(), group.name());
    EXPECT_EQ(s.operations().front(), identity_operation());
    for (const auto& d : s.dropped_operations()) EXPECT_FALSE(d.reason.empty());
    // With nphi a multiple of the group's multiple and nx = ny every operation is usable.
    if (c.nphi % s.phi_points_multiple() == 0 && c.nx == c.ny) {
      EXPECT_TRUE(s.dropped_operations().empty());
    }
  }
}

TEST(GridSymmetry, UsableSubgroupSizes) {
  const double pi = Constants::pi;
  const SurfaceLattice hex(4.84, 4.84, pi / 3.0);
  const SurfaceLattice square(4.0, 4.0, pi / 2.0);
  auto usable = [](const SurfaceLattice& l, const char* g, int nx, int ny, int nphi) {
    return GridSymmetry(l, PointGroup(g), nx, ny, nphi).operations().size();
  };
  EXPECT_EQ(usable(hex, "C6v", 8, 8, 12), 12U);
  EXPECT_EQ(usable(hex, "C6v", 8, 8, 6), 12U);
  EXPECT_EQ(usable(hex, "C6v", 8, 8, 4), 4U);   // E, C2 and the mirrors at 0 and 90 deg
  EXPECT_EQ(usable(hex, "C6v", 8, 8, 10), 4U);  // the hardwall/Dai-Light nphi without nPhisOptimal
  EXPECT_EQ(usable(hex, "C3v", 8, 8, 10), 2U);  // E and the mirror at 90 deg
  EXPECT_EQ(usable(hex, "C3v", 8, 8, 12), 6U);
  EXPECT_EQ(usable(hex, "C3", 8, 8, 4), 1U);    // trivial
  EXPECT_EQ(usable(hex, "C3", 8, 8, 6), 3U);
  EXPECT_EQ(usable(square, "C4v", 8, 8, 4), 8U);
  EXPECT_EQ(usable(square, "C4v", 8, 8, 6), 4U);
  EXPECT_EQ(usable(square, "C4v", 8, 6, 8), 4U);  // nx != ny: no 90 deg rotations, no diagonal mirrors
  EXPECT_EQ(usable(hex, "C6v", 8, 6, 12), 2U);    // nx != ny: only E and C2 (diagonal lattice matrices)
  EXPECT_TRUE(GridSymmetry(hex, PointGroup("C3"), 8, 8, 4).is_trivial());
}

// Every usable operation permutes the (x, y) points and the phi points; the grid action is a
// homomorphism; images agree with the geometric action (Cartesian rotation/reflection modulo the
// lattice, azimuth map modulo 2 pi).
TEST(GridSymmetry, GridActionIsTheGeometricAction) {
  for (const GridCase& c : grid_cases()) {
    const GridSymmetry s(c.group.lattice, PointGroup(c.group.group), c.nx, c.ny, c.nphi);
    const SurfaceLattice& l = c.group.lattice;
    const int n_ops = static_cast<int>(s.operations().size());
    const double dx = l.lattice_constant_x() / c.nx;
    const double dy = l.lattice_constant_y() / c.ny;
    const double dphi = Constants::two_pi / c.nphi;
    for (int iop = 0; iop < n_ops; ++iop) {
      const SurfaceOperation& g = s.operations()[iop];
      std::set<int> xy_images;
      for (int ix = 0; ix < c.nx; ++ix)
        for (int iy = 0; iy < c.ny; ++iy) {
          const int image = s.apply_to_xy(iop, ix, iy);
          xy_images.insert(image);
          const auto [x, y] = l.to_skewed(apply_to_cartesian(g, l.to_cartesian(ix * dx, iy * dy)));
          const int ix_image = image / c.ny;
          const int iy_image = image % c.ny;
          ASSERT_NEAR(wrapped(x / l.lattice_constant_x() - static_cast<double>(ix_image) / c.nx), 0.0, 1e-12);
          ASSERT_NEAR(wrapped(y / l.lattice_constant_y() - static_cast<double>(iy_image) / c.ny), 0.0, 1e-12);
          // Homomorphism on the grid.
          for (int jop = 0; jop < n_ops; ++jop) {
            const int composed = index_of(s.operations(), compose(s.operations()[jop], g));
            ASSERT_GE(composed, 0);
            ASSERT_EQ(s.apply_to_xy(composed, ix, iy), s.apply_to_xy(jop, ix_image, iy_image));
          }
        }
      EXPECT_EQ(static_cast<int>(xy_images.size()), c.nx * c.ny);
      const std::vector<int>& phi_index = s.phi_index(iop);
      EXPECT_EQ(std::set<int>(phi_index.begin(), phi_index.end()).size(), static_cast<std::size_t>(c.nphi));
      for (int i = 0; i < c.nphi; ++i)
        ASSERT_LT(angle_difference(phi_index[i] * dphi, apply_to_azimuth(g, i * dphi)), 1e-12);
    }
  }
}

// Orbits: a partition of the (x, y) points; the stored operation maps every point onto its
// representative; orbit sizes obey the orbit-stabiliser theorem; the number of orbits equals
// Burnside's count (1/|G|) sum_g fix(g).
TEST(GridSymmetry, OrbitsAndRepresentatives) {
  for (const GridCase& c : grid_cases()) {
    const GridSymmetry s(c.group.lattice, PointGroup(c.group.group), c.nx, c.ny, c.nphi);
    const int n_ops = static_cast<int>(s.operations().size());
    std::vector<int> orbit_size(s.n_representatives(), 0);
    for (int ixy = 0; ixy < s.n_xy(); ++ixy) {
      const int irep = s.representative_of(ixy);
      ASSERT_GE(irep, 0);
      ASSERT_LT(irep, s.n_representatives());
      ASSERT_EQ(s.apply_to_xy(s.operation_of(ixy), ixy / c.ny, ixy % c.ny), s.representative_xy(irep));
      ++orbit_size[irep];
    }
    for (int irep = 0; irep < s.n_representatives(); ++irep) {
      const int rep = s.representative_xy(irep);
      EXPECT_EQ(s.representative_of(rep), irep);
      EXPECT_EQ(s.operation_of(rep), 0);  // the identity
      if (irep > 0) {
        EXPECT_GT(rep, s.representative_xy(irep - 1));  // first of its orbit in memory order
      }
      int stabiliser = 0;
      for (int iop = 0; iop < n_ops; ++iop)
        if (s.apply_to_xy(iop, rep / c.ny, rep % c.ny) == rep) ++stabiliser;
      EXPECT_EQ(orbit_size[irep] * stabiliser, n_ops);
    }
    int fixed_points = 0;
    for (int iop = 0; iop < n_ops; ++iop)
      for (int ixy = 0; ixy < s.n_xy(); ++ixy)
        if (s.apply_to_xy(iop, ixy / c.ny, ixy % c.ny) == ixy) ++fixed_points;
    EXPECT_EQ(fixed_points % n_ops, 0);
    EXPECT_EQ(s.n_representatives(), fixed_points / n_ops);
  }
}

// Number of (x, y) points of the Fortran regions (1-based loops of setupReducedXYGrids*).
// Rectangle A (C2v) / square A (C4v, nphi not a multiple of 4), triangle 1 (C4v), triangle 1 and
// triangle A (C6v). Triangle 1 and the rectangles hold one point per orbit; triangle A holds the
// whole x = 0 column although (0, y) and (0, -y) are equivalent, so it has more points.
TEST(GridSymmetry, RepresentativeCountsAgainstTheFortranRegions) {
  const double pi = Constants::pi;
  for (const int n : {4, 6, 8, 10, 12, 16, 20, 24, 32}) {
    const int rectangle = (n / 2 + 1) * (n / 2 + 1);
    int square_triangle = 0;
    for (int ix = 1; ix <= n / 2 + 1; ++ix) square_triangle += ix;
    int hex_triangle_1 = 0;
    int hex_triangle_a = 0;
    for (int ix = 1; ix <= n / 2 + 1; ++ix) {
      hex_triangle_1 += ix <= n / 3 + 1 ? ix : n - 2 * ix + 3;
      hex_triangle_a += ix == 1 ? n : n - 2 * ix + 3;
    }
    const SurfaceLattice rectangular(6.2, 5.2, pi / 2.0);
    const SurfaceLattice square(4.0, 4.0, pi / 2.0);
    const SurfaceLattice hex(4.84, 4.84, pi / 3.0);
    EXPECT_EQ(GridSymmetry(rectangular, PointGroup("C2v"), n, n, 10).n_representatives(), rectangle) << n;
    EXPECT_EQ(GridSymmetry(square, PointGroup("C4v"), n, n, 6).n_representatives(), rectangle) << n;
    EXPECT_EQ(GridSymmetry(square, PointGroup("C4v"), n, n, 8).n_representatives(), square_triangle) << n;
    EXPECT_EQ(GridSymmetry(hex, PointGroup("C6v"), n, n, 12).n_representatives(), hex_triangle_1) << n;
    EXPECT_LE(GridSymmetry(hex, PointGroup("C6v"), n, n, 10).n_representatives(), hex_triangle_a) << n;
  }
  // The values quoted in the review of the plan: n = 8 -> 10 and 24 (21 orbits); n = 32 -> 102.
  const SurfaceLattice hex(4.84, 4.84, pi / 3.0);
  EXPECT_EQ(GridSymmetry(hex, PointGroup("C6v"), 8, 8, 12).n_representatives(), 10);
  EXPECT_EQ(GridSymmetry(hex, PointGroup("C6v"), 8, 8, 10).n_representatives(), 21);
  EXPECT_EQ(GridSymmetry(hex, PointGroup("C6v"), 32, 32, 12).n_representatives(), 102);
}

TEST(GridSymmetry, NoSymmetry) {
  const GridSymmetry s(4, 6, 10);
  EXPECT_TRUE(s.is_trivial());
  EXPECT_EQ(s.point_group(), "none");
  EXPECT_EQ(s.n_representatives(), 24);
  for (int ixy = 0; ixy < 24; ++ixy) {
    EXPECT_EQ(s.representative_xy(s.representative_of(ixy)), ixy);
    EXPECT_EQ(s.operation_of(ixy), 0);
  }
  for (int i = 0; i < 10; ++i) EXPECT_EQ(s.phi_index(0)[i], i);
}

// A potential that is exactly symmetric on the grid, built independently of GridSymmetry: the
// label of a (x, y, phi) grid point is the smallest flat index among its geometric images
// (Cartesian rotation/reflection, rounding back to the grid). Gathering the representative values
// through the maps reproduces it bitwise at every point.
TEST(GridSymmetry, GatherReproducesAnExactlySymmetricFunction) {
  for (const GridCase& c : grid_cases()) {
    if (c.nx > 12 || c.nphi > 14) continue;  // keep the brute force small
    const GridSymmetry s(c.group.lattice, PointGroup(c.group.group), c.nx, c.ny, c.nphi);
    const SurfaceLattice& l = c.group.lattice;
    const double dx = l.lattice_constant_x() / c.nx;
    const double dy = l.lattice_constant_y() / c.ny;
    const double dphi = Constants::two_pi / c.nphi;
    auto flat = [&](int ixy, int iphi) { return iphi + c.nphi * ixy; };
    auto label = [&](int ix, int iy, int iphi) {
      int smallest = flat(iy + c.ny * ix, iphi);
      for (const SurfaceOperation& g : s.operations()) {
        const auto [x, y] = l.to_skewed(apply_to_cartesian(g, l.to_cartesian(ix * dx, iy * dy)));
        const int jx = static_cast<int>(std::lround(x / dx)) % c.nx;
        const int jy = static_cast<int>(std::lround(y / dy)) % c.ny;
        const int jphi = static_cast<int>(std::lround(apply_to_azimuth(g, iphi * dphi) / dphi)) % c.nphi;
        const int j = flat((jy + c.ny) % c.ny + c.ny * ((jx + c.nx) % c.nx), (jphi + c.nphi) % c.nphi);
        smallest = std::min(smallest, j);
      }
      return smallest;
    };
    auto f = [](int label_value) { return std::sin(0.37 * label_value) + 1e-3 * label_value; };

    // Values at the representatives only (phi fastest), then gathered.
    std::vector<double> reduced(static_cast<std::size_t>(s.n_representatives()) * c.nphi);
    for (int irep = 0; irep < s.n_representatives(); ++irep) {
      const int rep = s.representative_xy(irep);
      for (int iphi = 0; iphi < c.nphi; ++iphi) reduced[iphi + c.nphi * irep] = f(label(rep / c.ny, rep % c.ny, iphi));
    }
    for (int ixy = 0; ixy < s.n_xy(); ++ixy) {
      const std::vector<int>& phi_index = s.phi_index(s.operation_of(ixy));
      for (int iphi = 0; iphi < c.nphi; ++iphi) {
        const double gathered = reduced[phi_index[iphi] + c.nphi * s.representative_of(ixy)];
        ASSERT_EQ(gathered, f(label(ixy / c.ny, ixy % c.ny, iphi)))
            << c.group.group << " " << c.nx << "x" << c.ny << " nphi " << c.nphi << " ixy " << ixy << " iphi " << iphi;
      }
    }
  }
}
