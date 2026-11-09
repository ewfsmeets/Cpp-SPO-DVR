// Tests: the atom-surface PES (src/crplib/atom_surface*.cpp, hexagonal_ and square_atom_surface.cpp)
// on synthetic data sets: bitwise equality with the Fortran library, exactness at the data points,
// periodicity and point-group symmetry, finite-difference gradients, the unit-cell folds.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/atom_surface.hpp"
#include "crplib/atom_surface_detail.hpp"
#include "crplib/errors.hpp"
#include "fortran_reference.hpp"
#include "reference_basis.hpp"

using namespace crplib;
using namespace crplib::detail;
using crplib_test::DriverInput;
using crplib_test::DriverOutput;
using crplib_test::reference_file;

namespace {

const double sqrt3 = std::sqrt(3.0);

// Cartesian of lattice coordinates (u, v) (skewed for the hexagonal symmetries).
Vec3 at(Symmetry s, double u, double v, double z, double a) {
  if (s == Symmetry::c4v) return {u * a, v * a, z};
  return {(u + v / 2.0) * a, v * sqrt3 / 2.0 * a, z};
}

// Rotation by `degrees` about the origin (a top site) and the mirror x -> -x.
Vec3 rotated(const Vec3& p, double degrees) {
  const double c = std::cos(degrees * std::numbers::pi / 180.0), s = std::sin(degrees * std::numbers::pi / 180.0);
  return {c * p.x - s * p.y, s * p.x + c * p.y, p.z};
}
Vec3 mirrored(const Vec3& p) { return {-p.x, p.y, p.z}; }

}  // namespace

// ---------------------------------------------------------------------------------------------
// Bitwise equality with Calculate_CRP_AtomSurface_PES (internal units) of the Fortran library.
// ---------------------------------------------------------------------------------------------

class FortranAtomSurface : public testing::TestWithParam<std::string> {};

TEST_P(FortranAtomSurface, BitwiseEqual) {
  const auto output_file = reference_file(GetParam() + "/atom_0_2.out");
  REQUIRE_CRP_REFERENCE(output_file);
  const DataSet data = read_data_set(reference_file(GetParam() + "/dataset"));
  const auto pes = make_atom_surface(data);
  DriverInput in(reference_file(GetParam() + "/atom_0_2.in"));
  DriverOutput out(output_file);
  const int n = in.integer();
  for (int i = 0; i < n; ++i) {
    const Vec3 p{in.real(), in.real(), in.real()};
    const AtomSurfaceValue v = pes->evaluate(p);
    const auto value = out.next("V");
    const auto gradient = out.next_bits("F", 3);
    EXPECT_BITS(v.v, value.bits) << " at (" << p.x << ", " << p.y << ", " << p.z << ")";
    EXPECT_BITS(v.gradient.x, gradient[0]) << " at (" << p.x << ", " << p.y << ", " << p.z << ")";
    EXPECT_BITS(v.gradient.y, gradient[1]) << " at (" << p.x << ", " << p.y << ", " << p.z << ")";
    EXPECT_BITS(v.gradient.z, gradient[2]) << " at (" << p.x << ", " << p.y << ", " << p.z << ")";
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, FortranAtomSurface,
                         testing::Values("c6v_pairwise", "c3v_pairwise", "c4v_pairwise", "c6v_pairwise_angstrom",
                                         "c4v_negative_spacing"));

// ---------------------------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------------------------

class AtomSurfaceProperties : public testing::TestWithParam<Symmetry> {
 protected:
  void SetUp() override {
    data_ = crplib_test::pairwise_data_set(GetParam());
    pes_ = make_atom_surface(data_);
  }
  [[nodiscard]] double a() const { return data_.lattice_constant; }
  DataSet data_;
  std::unique_ptr<const AtomSurface> pes_;
};

// The interpolation reproduces the cut data at the sites and z knots (the Fourier table is the
// inverse of the basis at the sites; C6v only to ~1e-13 because of its truncated DATA literals).
TEST_P(AtomSurfaceProperties, ExactAtTheSitesAndKnots) {
  for (const AtomCut& cut : data_.atom_cuts) {
    for (std::size_t i = 0; i < cut.potential.x.size(); ++i) {
      const Vec3 p{cut.x * a(), cut.y * a(), cut.potential.x[i]};
      EXPECT_NEAR(pes_->evaluate(p).v, cut.potential.f[i], 2e-12 * std::max(1.0, std::abs(cut.potential.f[i])))
          << cut.file_name << " z = " << p.z;
    }
  }
}

TEST_P(AtomSurfaceProperties, PeriodicAndSymmetric) {
  const Symmetry s = GetParam();
  // C6v: 60 degree rotations and the mirror; C3v: 120 degrees and the mirror at 90 degrees (x -> -x);
  // C4v: 90 degree rotations and the mirror.
  const double step = s == Symmetry::c6v ? 60.0 : (s == Symmetry::c3v ? 120.0 : 90.0);
  for (double z : {0.2, 1.3, 3.1, 5.5}) {
    for (const auto& [u, v] : {std::pair{0.13, 0.29}, {0.41, 0.07}, {0.71, 0.66}, {0.05, 0.52}}) {
      const Vec3 p = at(s, u, v, z, a());
      const double v0 = pes_->evaluate(p).v;
      const double tolerance = 1e-12 * std::max(1.0, std::abs(v0));
      EXPECT_NEAR(pes_->evaluate(at(s, u + 1, v, z, a())).v, v0, tolerance);
      EXPECT_NEAR(pes_->evaluate(at(s, u - 2, v + 1, z, a())).v, v0, tolerance);
      for (double angle = step; angle < 360.0; angle += step)
        EXPECT_NEAR(pes_->evaluate(rotated(p, angle)).v, v0, tolerance);
      EXPECT_NEAR(pes_->evaluate(mirrored(p)).v, v0, tolerance);
    }
  }
}

TEST_P(AtomSurfaceProperties, GradientMatchesFiniteDifferences) {
  const double h = 1e-5;
  for (const auto& [u, v, z] : {std::tuple{0.13, 0.29, 0.6}, {0.41, 0.07, 2.2}, {0.71, 0.66, 3.9}, {0.35, 0.12, 1.1}}) {
    const Vec3 p = at(GetParam(), u, v, z, a());
    const AtomSurfaceValue value = pes_->evaluate(p);
    auto f = [&](double dx, double dy, double dz) { return pes_->evaluate({p.x + dx, p.y + dy, p.z + dz}).v; };
    EXPECT_NEAR(value.gradient.x, (f(h, 0, 0) - f(-h, 0, 0)) / (2 * h), 2e-7);
    EXPECT_NEAR(value.gradient.y, (f(0, h, 0) - f(0, -h, 0)) / (2 * h), 2e-7);
    EXPECT_NEAR(value.gradient.z, (f(0, 0, h) - f(0, 0, -h)) / (2 * h), 2e-7);
  }
}

TEST_P(AtomSurfaceProperties, OutsideTheDataRange) {
  const std::vector<double>& z = data_.atom_cuts.front().potential.x;
  EXPECT_THROW((void)pes_->evaluate({0.3, 0.2, z.front() - 0.01}), EvaluationError);
  if (GetParam() == Symmetry::c4v) {
    EXPECT_EQ(pes_->evaluate({0.3, 0.2, z.back()}).v, 0.0);  // Fortran: Z0 >= Z(NZ) -> 0
    EXPECT_EQ(pes_->evaluate({0.3, 0.2, z.back() + 1.0}).v, 0.0);
  } else {
    EXPECT_THROW((void)pes_->evaluate({0.3, 0.2, z.back() + 0.01}), EvaluationError);
  }
}

TEST_P(AtomSurfaceProperties, NoCoverageWarningForTheSyntheticSets) { EXPECT_TRUE(pes_->warnings().empty()); }

INSTANTIATE_TEST_SUITE_P(Symmetries, AtomSurfaceProperties,
                         testing::Values(Symmetry::c3v, Symmetry::c4v, Symmetry::c6v));

// The cluster sums equal the sum of the corrugation function over the infinite surface for points
// in the folded unit cell (the cut-off is inside the coverage radius).
TEST(AtomSurfaceClusters, EqualTheInfiniteLatticeSum) {
  for (Symmetry s : {Symmetry::c3v, Symmetry::c4v, Symmetry::c6v}) {
    const DataSet data = crplib_test::pairwise_data_set(s);
    const CubicSpline1D corrugation = corrugation_function(data);
    const double a = data.lattice_constant;
    const bool square = s == Symmetry::c4v;
    const int layers = s == Symmetry::c3v ? 4 : (square ? 2 : 1);
    auto infinite_sum = [&](double x, double y, double z) {
      double sum = 0.0, layer_z = 0.0;
      for (int k = 0; k < layers; ++k) {
        if (k > 0) layer_z -= data.interlayer_spacings[k - 1];
        const double shift = square ? 0.5 * k : (k % 3) / 3.0;
        for (int i = -12; i <= 12; ++i) {
          for (int j = -12; j <= 12; ++j) {
            const Vec3 atom = at(s, i + shift, j + shift, layer_z, a);
            const double r = std::hypot(x - atom.x, y - atom.y, z - atom.z);
            if (r <= corrugation.x().back()) sum += corrugation.value_and_slope_splint(r).value;
          }
        }
      }
      return sum;
    };
    for (const auto& [u, v, z] : {std::tuple{0.0, 0.0, 0.5}, {0.3, 0.1, 1.2}, {0.9, 0.8, 0.4}, {0.5, 0.5, 2.5}}) {
      double x = at(s, u, v, z, a).x, y = at(s, u, v, z, a).y;
      if (square) {
        // MURO is evaluated in the irreducible triangle
        const auto folded = fold_square(x, y, a);
        x = folded.first;
        y = folded.second;
        (void)fold_into_triangle(x, y, a);
        const double cluster = SquareAtomSurface(data).repulsive_sum(x, y, z).v;
        EXPECT_NEAR(cluster, infinite_sum(x, y, z), 1e-15) << u << ", " << v;
      } else {
        const double cluster = HexagonalAtomSurface(data).repulsive_sum(x, y, z).v;
        EXPECT_NEAR(cluster, infinite_sum(x, y, z), 1e-15) << u << ", " << v;
      }
    }
  }
}

TEST(AtomSurfaceClusters, CoverageRadius) {
  const double c6v = HexagonalAtomSurface(crplib_test::pairwise_data_set(Symmetry::c6v)).coverage_radius() / 4.8;
  const double c3v = HexagonalAtomSurface(crplib_test::pairwise_data_set(Symmetry::c3v)).coverage_radius() / 4.82;
  const double c4v = SquareAtomSurface(crplib_test::pairwise_data_set(Symmetry::c4v)).coverage_radius() / 4.82;
  EXPECT_NEAR(c6v, 2.60, 0.01);
  EXPECT_NEAR(c3v, 2.89, 0.01);
  EXPECT_NEAR(c4v, 2.06, 0.01);
  // A corrugation function reaching beyond it gives a warning.
  DataSet data = crplib_test::pairwise_data_set(Symmetry::c6v);
  for (double& x : data.corrugation.x) x *= 3.0;
  EXPECT_EQ(make_atom_surface(data)->warnings().size(), 1u);
}

TEST(AtomSurfaceFolds, HexagonalCellAndTheNegativeIntegerQuirk) {
  const double a = 4.8;
  const auto [x, y] = fold_hexagonal(2.3 * a, -1.7 * a, a);
  const double u = x / a - y / (a * sqrt3 / 2.0) / 2.0, v = y / (a * sqrt3 / 2.0);
  EXPECT_GE(u, 0.0);
  EXPECT_LT(u, 1.0);
  EXPECT_GE(v, 0.0);
  EXPECT_LT(v, 1.0);
  // u = -1 exactly: INT(-1 - 1) = -2, so the point lands on the far edge u = 1, not on u = 0.
  EXPECT_EQ(fold_hexagonal(-a, 0.0, a).first, a);
  EXPECT_EQ(fold_hexagonal(a, 0.0, a).first, 0.0);
  EXPECT_EQ(fold_square(-a, 0.0, a).first, a);
  EXPECT_THROW((void)fold_square(1e12, 0.0, a), EvaluationError);
}

TEST(AtomSurfaceFolds, RegionsMapIntoTheTriangle) {
  const double a = 1.0;
  int seen[9] = {};
  for (double x = 0.0; x < 1.0; x += 0.0625) {
    for (double y = 0.0; y < 1.0; y += 0.0625) {
      double fx = x, fy = y;
      const int region = fold_into_triangle(fx, fy, a);
      ASSERT_GE(region, 1);
      ASSERT_LE(region, 8);
      ++seen[region];
      EXPECT_GE(fy, -1e-15);
      EXPECT_LE(fy, fx + 1e-15);
      EXPECT_LE(fx, 0.5 + 1e-15);
      // the gradient map undoes the fold: transform a test vector both ways
      double gx = 0.3, gy = -0.7;
      unfold_gradient(region, gx, gy);
      EXPECT_NEAR(std::hypot(gx, gy), std::hypot(0.3, 0.7), 1e-15);
    }
  }
  for (int r = 1; r <= 8; ++r) EXPECT_GT(seen[r], 0) << "region " << r;
  // Boundaries as the Fortran decides them
  double x = 0.25, y = 0.25;
  EXPECT_EQ(fold_into_triangle(x, y, a), 1);  // x >= y
  x = 0.5, y = 0.5;
  EXPECT_EQ(fold_into_triangle(x, y, a), 1);  // y <= a/2, x <= a/2
  x = 0.75, y = 0.25;
  EXPECT_EQ(fold_into_triangle(x, y, a), 4);  // y >= a - x
}

TEST(AtomSurfaceFolds, DistanceToPolygon) {
  const std::vector<std::pair<double, double>> square = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  EXPECT_EQ(distance_to_polygon(0.5, 0.5, square), 0.0);
  EXPECT_NEAR(distance_to_polygon(2.0, 0.5, square), 1.0, 1e-15);
  EXPECT_NEAR(distance_to_polygon(2.0, 2.0, square), std::sqrt(2.0), 1e-15);
}

// A data set in the span of the interpolation: atom cuts = (the library's own repulsive sum at the
// site) + sum_j c_j(z) f_j(site) with cubic c_j and exact end slopes (C4v: c_0(z) only, the
// periodic spline reproduces constants in x, y). The PES must then equal
// sum_j c_j(z) F_j(x, y) + repulsive sum at every point, not only at the data points; the basis
// F_j is the independent one of reference_basis.hpp.
TEST(AtomSurfaceRepresentable, ReproducedOffTheGrid) {
  for (Symmetry s : {Symmetry::c3v, Symmetry::c4v, Symmetry::c6v}) {
    DataSet data = crplib_test::pairwise_data_set(s);
    const double a = data.lattice_constant;
    const bool square = s == Symmetry::c4v;
    const bool c3v = s == Symmetry::c3v;
    const int n_basis = square ? 1 : (c3v ? 10 : 7);
    auto c = [](int j, double z) {
      return 0.002 * (j % 2 ? -1.0 : 1.0) / (j + 1) * (1.0 + 0.3 * z - 0.05 * z * z + 0.004 * z * z * z);
    };
    auto dc = [](int j, double z) { return 0.002 * (j % 2 ? -1.0 : 1.0) / (j + 1) * (0.3 - 0.1 * z + 0.012 * z * z); };
    auto basis_at = [&](double x, double y) {
      if (square) return std::vector<double>{1.0};
      const crplib_test::Site site = crplib_test::skewed(x, y, a);
      return crplib_test::hexagonal_atom_basis(site.u, site.v, c3v);
    };
    // The library's repulsive sum (it only depends on the corrugation function and the slab).
    auto repulsive = [&](double x, double y, double z) {
      if (square) return SquareAtomSurface(data).repulsive_sum(x, y, z);
      return HexagonalAtomSurface(data).repulsive_sum(x, y, z);
    };
    for (AtomCut& cut : data.atom_cuts) {
      const double x = cut.x * a, y = cut.y * a;
      const std::vector<double> f = basis_at(x, y);
      for (std::size_t i = 0; i < cut.potential.x.size(); ++i) {
        const double z = cut.potential.x[i];
        double value = 0.0;
        for (int j = 0; j < n_basis; ++j) value += c(j, z) * f[j];
        cut.potential.f[i] = repulsive(x, y, z).v + value;
      }
      auto slope = [&](double z) {
        double value = 0.0;
        for (int j = 0; j < n_basis; ++j) value += dc(j, z) * f[j];
        return repulsive(x, y, z).dvdz + value;
      };
      cut.slope_first = slope(cut.potential.x.front());
      cut.slope_last = slope(cut.potential.x.back());
    }
    const auto pes = make_atom_surface(data);
    for (const auto& [u, v, z] : {std::tuple{0.137, 0.291, 0.61},
                                  {0.413, 0.072, 2.23},
                                  {0.717, 0.659, 3.91},
                                  {0.355, 0.121, 1.17},
                                  {-1.3, 2.7, 5.3},
                                  {0.5, 0.5, -0.4}}) {
      const Vec3 p = at(s, u, v, z, a);
      double x0 = p.x, y0 = p.y;
      if (square) {
        std::tie(x0, y0) = fold_square(p.x, p.y, a);
        (void)fold_into_triangle(x0, y0, a);
      } else {
        std::tie(x0, y0) = fold_hexagonal(p.x, p.y, a);
      }
      const std::vector<double> f = basis_at(p.x, p.y);
      double expected = repulsive(x0, y0, z).v;
      double smooth = 0.0;
      for (int j = 0; j < n_basis; ++j) smooth += c(j, z) * f[j];
      expected += smooth;
      EXPECT_NEAR(pes->evaluate(p).v, expected, 1e-13)
          << crplib::symmetry_name(s) << " at " << u << ", " << v << ", " << z;
    }
  }
}
