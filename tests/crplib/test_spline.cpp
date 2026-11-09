// Tests: cubic splines and tridiagonal solvers of the CRP library (src/crplib/spline.*).
// Properties (exact reproduction of cubic polynomials, solver accuracy) and bitwise equality with
// the Fortran routines (reference files from tools/fortran_reference/run_crp_reference.sh kernels).
#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "crplib/errors.hpp"
#include "crplib/spline.hpp"
#include "fortran_reference.hpp"

using namespace crplib::detail;
using crplib::DataError;
using crplib_test::DriverInput;
using crplib_test::DriverOutput;
using crplib_test::reference_file;

namespace {

// A fixed sequence of numbers in [-1, 1) (no standard-library distribution: its output depends on
// the implementation).
std::vector<double> fixed_sequence(int n, int seed) {
  std::vector<double> values(n);
  unsigned state = 12345u + static_cast<unsigned>(seed) * 7919u;
  for (double& v : values) {
    state = state * 1103515245u + 12345u;
    v = static_cast<double>((state >> 8) % 20000u) / 10000.0 - 1.0;
  }
  return values;
}

double cubic(double x) { return 0.3 - 1.2 * x + 0.7 * x * x - 0.25 * x * x * x; }
double cubic_slope(double x) { return -1.2 + 1.4 * x - 0.75 * x * x; }
double cubic_curvature(double x) { return 1.4 - 1.5 * x; }

const std::vector<double> knots = {-1.0, -0.4, 0.1, 0.35, 0.9, 1.6, 2.0, 2.8, 3.5};

template <typename Function>
std::vector<double> values_of(Function f, const std::vector<double>& x) {
  std::vector<double> y;
  for (double v : x) y.push_back(f(v));
  return y;
}

SplineEnd end_for_code(int code, double value) {
  if (code == 0) return {EndCondition::not_a_knot, value};
  if (code == 1) return {EndCondition::second_derivative, value};
  if (code == 2) return {EndCondition::first_derivative, value};
  return {EndCondition::natural, value};  // Fortran: any code > 2
}

}  // namespace

TEST(Tridiagonal, SolvesAgainstADenseResidual) {
  const int n = 7;
  const std::vector<double> lower = fixed_sequence(n, 1), upper = fixed_sequence(n, 2), rhs0 = fixed_sequence(n, 3);
  std::vector<double> diagonal0 = fixed_sequence(n, 4);
  for (double& d : diagonal0) d += 4.0;  // diagonally dominant
  std::vector<double> diagonal = diagonal0, rhs = rhs0, x(n);
  solve_tridiagonal(diagonal, lower, upper, rhs, x);
  for (int i = 0; i < n; ++i) {
    double row = diagonal0[i] * x[i];
    if (i > 0) row += lower[i - 1] * x[i - 1];
    if (i + 1 < n) row += upper[i] * x[i + 1];
    EXPECT_NEAR(row, rhs0[i], 1e-14);
  }
}

TEST(Tridiagonal, CyclicSolvesAgainstADenseResidual) {
  const int n = 6;
  const std::vector<double> lower = fixed_sequence(n, 5), upper = fixed_sequence(n, 6), rhs0 = fixed_sequence(n, 7);
  std::vector<double> diagonal0 = fixed_sequence(n, 8);
  for (double& d : diagonal0) d += 4.0;
  const double bottom_left = 0.37, top_right = -0.61;
  std::vector<double> diagonal = diagonal0, rhs = rhs0, x(n);
  solve_cyclic_tridiagonal(diagonal, lower, upper, bottom_left, top_right, rhs, x);
  for (int i = 0; i < n; ++i) {
    double row = diagonal0[i] * x[i];
    if (i > 0) row += lower[i - 1] * x[i - 1];
    if (i + 1 < n) row += upper[i] * x[i + 1];
    if (i == 0) row += top_right * x[n - 1];
    if (i == n - 1) row += bottom_left * x[0];
    EXPECT_NEAR(row, rhs0[i], 1e-14);
  }
}

// Not-a-knot, clamped (exact slopes) and second-derivative (exact curvatures) splines all reproduce
// a cubic polynomial: second derivatives at the knots, values and slopes everywhere inside.
TEST(CubicSpline1D, ReproducesCubicPolynomials) {
  const std::vector<double> y = values_of(cubic, knots);
  const double a = knots.front(), b = knots.back();
  for (const auto& [first, last] : {std::pair{not_a_knot, not_a_knot},
                                    std::pair{SplineEnd{EndCondition::first_derivative, cubic_slope(a)},
                                              SplineEnd{EndCondition::first_derivative, cubic_slope(b)}},
                                    std::pair{SplineEnd{EndCondition::second_derivative, cubic_curvature(a)},
                                              SplineEnd{EndCondition::second_derivative, cubic_curvature(b)}},
                                    std::pair{not_a_knot, SplineEnd{EndCondition::first_derivative, cubic_slope(b)}}}) {
    const CubicSpline1D s(knots, y, first, last);
    for (std::size_t i = 0; i < knots.size(); ++i)
      EXPECT_NEAR(s.second_derivatives()[i], cubic_curvature(knots[i]), 1e-12);
    for (double t = a; t <= b; t += 0.0625) {
      EXPECT_NEAR(s.value_dpl(t), cubic(t), 1e-13) << t;
      EXPECT_NEAR(s.slope_dplp(t), cubic_slope(t), 1e-12) << t;
      const auto vs = s.value_and_slope_splint(t);
      EXPECT_NEAR(vs.value, cubic(t), 1e-13) << t;
      EXPECT_NEAR(vs.slope, cubic_slope(t), 1e-12) << t;
    }
  }
}

TEST(CubicSpline1D, NaturalSplineHasZeroEndCurvature) {
  const CubicSpline1D s(knots, values_of(cubic, knots), {EndCondition::natural, 5.0}, {EndCondition::natural, 5.0});
  EXPECT_EQ(s.second_derivatives().front(), 0.0);
  EXPECT_EQ(s.second_derivatives().back(), 0.0);
  for (std::size_t i = 0; i < knots.size(); ++i) EXPECT_NEAR(s.value_dpl(knots[i]), cubic(knots[i]), 1e-14);
}

// DPL continues linearly outside the knots, SPLINT continues the end cubics.
TEST(CubicSpline1D, ExtrapolationOutsideTheKnots) {
  const CubicSpline1D s(knots, values_of(cubic, knots), not_a_knot, not_a_knot);
  const double below = knots.front() - 0.5, above = knots.back() + 0.5;
  const double slope_low = (s.value_dpl(knots.front() - 0.25) - s.value_dpl(below)) / 0.25;
  EXPECT_NEAR(s.value_dpl(knots.front()) - s.value_dpl(knots.front() - 0.25), slope_low * 0.25, 1e-12);
  EXPECT_NEAR(s.value_and_slope_splint(below).value, cubic(below), 1e-12);  // the cubic itself
  EXPECT_NEAR(s.value_and_slope_splint(above).value, cubic(above), 1e-12);
  EXPECT_GT(std::abs(s.value_dpl(above) - cubic(above)), 1e-3);  // linear, not cubic
}

TEST(CubicSpline1D, RejectsTooFewOrNonMonotonicKnots) {
  EXPECT_THROW(CubicSpline1D({0.0, 1.0, 2.0}, {1.0, 2.0, 0.0}, not_a_knot, not_a_knot), DataError);
  EXPECT_THROW(CubicSpline1D({0.0, 1.0, 0.5, 2.0, 3.0}, {1.0, 2.0, 0.0, 1.0, 2.0}, not_a_knot, not_a_knot), DataError);
  EXPECT_THROW(CubicSpline1D({0.0, 1.0, 1.0, 2.0, 3.0}, {1.0, 2.0, 0.0, 1.0, 2.0}, not_a_knot, not_a_knot), DataError);
}

TEST(LocateInterval, EdgesAndOutside) {
  const std::vector<double> x = {0.5, 1.0, 1.25, 2.0};
  EXPECT_FALSE(locate_interval(x, 0.4999).has_value());
  EXPECT_EQ(locate_interval(x, 0.5), 0);
  EXPECT_EQ(locate_interval(x, 1.0), 1);  // a knot starts its interval
  EXPECT_EQ(locate_interval(x, 1.9), 2);
  EXPECT_EQ(locate_interval(x, 2.0), 2);  // the last knot closes the last interval
  EXPECT_FALSE(locate_interval(x, 2.0000001).has_value());
}

// Not-a-knot in both directions: the product of two cubics is reproduced exactly, with derivatives.
TEST(BicubicSpline, ReproducesProductsOfCubics) {
  const std::vector<double> gx = {0.0, 0.3, 0.8, 1.0, 1.7, 2.1};
  const std::vector<double> gy = {-1.0, -0.2, 0.5, 1.5, 2.0};
  auto q = [](double y) { return 1.0 + 0.5 * y - 0.3 * y * y + 0.1 * y * y * y; };
  auto dq = [](double y) { return 0.5 - 0.6 * y + 0.3 * y * y; };
  std::vector<double> f;
  for (double y : gy)
    for (double x : gx) f.push_back(cubic(x) * q(y));
  const BicubicSpline s(gx, gy, f);
  for (double x = 0.0; x <= 2.1; x += 0.13) {
    for (double y = -1.0; y <= 2.0; y += 0.17) {
      const auto v = s.evaluate(x, y);
      ASSERT_TRUE(v.has_value());
      EXPECT_NEAR(v->value, cubic(x) * q(y), 1e-13);
      EXPECT_NEAR(v->d_dx, cubic_slope(x) * q(y), 1e-12);
      EXPECT_NEAR(v->d_dy, cubic(x) * dq(y), 1e-12);
    }
  }
  EXPECT_FALSE(s.evaluate(-0.01, 0.0).has_value());
  EXPECT_FALSE(s.evaluate(1.0, 2.01).has_value());
}

TEST(BicubicSpline, RejectsBadGrids) {
  const std::vector<double> f(4 * 4, 1.0);
  EXPECT_THROW(BicubicSpline({0, 1, 2}, {0, 1, 2, 3}, std::vector<double>(12, 1.0)), DataError);  // < 4 points
  EXPECT_THROW(BicubicSpline({0, 2, 1, 3}, {0, 1, 2, 3}, f), DataError);                          // not increasing
  EXPECT_THROW(BicubicSpline({3, 2, 1, 0}, {0, 1, 2, 3}, f), DataError);                          // decreasing
  EXPECT_THROW(BicubicSpline({0, 1, 2, 3}, {0, 1, 2, 3}, std::vector<double>(15, 1.0)), DataError);
}

// The periodic spline's second derivatives solve the continuity equations of F' at every knot,
// including the knot where the period closes (x_0 = x_{n-1}).
TEST(PeriodicSpline, SatisfiesThePeriodicContinuityEquations) {
  const std::vector<double> x = {0.0, 0.5, 1.2, 1.6, 2.4, 3.0};
  std::vector<double> y;
  for (double v : x) y.push_back(std::cos(2.0 * M_PI * v / 3.0) + 0.2 * std::sin(4.0 * M_PI * v / 3.0));
  y.back() = y.front();
  const std::vector<double> m = periodic_spline_second_derivatives(x, y);
  EXPECT_EQ(m.front(), m.back());
  const int n = static_cast<int>(x.size());
  const int intervals = n - 1;
  auto h = [&](int i) { return x[(i % intervals) + 1] - x[i % intervals]; };  // interval i (periodic)
  auto slope = [&](int i) { return (y[(i % intervals) + 1] - y[i % intervals]) / h(i); };
  auto curvature = [&](int knot) { return m[knot % intervals]; };
  for (int i = 1; i <= intervals; ++i) {  // knot i between intervals i-1 and i
    const double lhs =
        h(i - 1) / 6.0 * curvature(i - 1) + (h(i - 1) + h(i)) / 3.0 * curvature(i) + h(i) / 6.0 * curvature(i + 1);
    EXPECT_NEAR(lhs, slope(i) - slope(i - 1), 1e-13) << "knot " << i;
  }
}

// ---------------------------------------------------------------------------------------------
// Bitwise equality with the Fortran library (Spline.f90 through the reference driver).
// ---------------------------------------------------------------------------------------------

class FortranSpline1D : public testing::TestWithParam<std::string> {};

TEST_P(FortranSpline1D, BitwiseEqual) {
  const auto input_file = reference_file("kernels/spline1d_" + GetParam() + ".in");
  const auto output_file = reference_file("kernels/spline1d_" + GetParam() + ".out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(input_file);
  const int n = in.integer(), ic1 = in.integer(), icn = in.integer(), m = in.integer();
  const double cm1 = in.real(), cmn = in.real();
  const std::vector<double> x = in.reals(n), y = in.reals(n), t = in.reals(m);

  const CubicSpline1D s(x, y, end_for_code(ic1, cm1), end_for_code(icn, cmn));
  DriverOutput out(output_file);
  const auto cm = out.next_bits("cm", n);
  for (int i = 0; i < n; ++i) EXPECT_BITS(s.second_derivatives()[i], cm[i]) << " i = " << i;
  for (int q = 0; q < m; ++q) {
    const auto fortran = out.next_bits("Q", 4);  // DPL, DPLP, SPLINT value and slope
    const auto vs = s.value_and_slope_splint(t[q]);
    EXPECT_BITS(s.value_dpl(t[q]), fortran[0]) << " t = " << t[q];
    EXPECT_BITS(s.slope_dplp(t[q]), fortran[1]) << " t = " << t[q];
    EXPECT_BITS(vs.value, fortran[2]) << " t = " << t[q];
    EXPECT_BITS(vs.slope, fortran[3]) << " t = " << t[q];
  }
}

INSTANTIATE_TEST_SUITE_P(Kernels, FortranSpline1D,
                         testing::Values("not_a_knot", "clamped", "zdep", "second_derivative", "natural",
                                         "decreasing"));

TEST(FortranSpline2D, BitwiseEqual) {
  const auto output_file = reference_file("kernels/spline2d.out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(reference_file("kernels/spline2d.in"));
  const int nx = in.integer(), ny = in.integer(), m = in.integer();
  const std::vector<double> gx = in.reals(nx), gy = in.reals(ny), f = in.reals(nx * ny);
  const BicubicSpline s(gx, gy, f);
  DriverOutput out(output_file);
  const auto cx = out.next_bits("cx", nx * ny), cy = out.next_bits("cy", nx * ny), cxy = out.next_bits("cxy", nx * ny);
  for (int i = 0; i < nx * ny; ++i) {
    EXPECT_BITS(s.d2f_dx2()[i], cx[i]) << " i = " << i;
    EXPECT_BITS(s.d2f_dy2()[i], cy[i]) << " i = " << i;
    EXPECT_BITS(s.d4f_dx2dy2()[i], cxy[i]) << " i = " << i;
  }
  for (int q = 0; q < m; ++q) {
    const double x = in.real(), y = in.real();
    const auto v = s.evaluate(x, y);
    const auto& first = out.next("Q");
    if (first.text == "Q outside") {
      EXPECT_FALSE(v.has_value()) << x << ", " << y;
      continue;
    }
    ASSERT_TRUE(v.has_value()) << x << ", " << y;
    EXPECT_BITS(v->value, first.bits);
    EXPECT_BITS(v->d_dy, out.next("Q").bits);
    EXPECT_BITS(v->d_dx, out.next("Q").bits);
  }
}

TEST(FortranTrispl, BitwiseEqual) {
  const auto output_file = reference_file("kernels/trispl.out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(reference_file("kernels/trispl.in"));
  const int n = in.integer();
  const std::vector<double> x = in.reals(n), y = in.reals(n);
  const std::vector<double> m = periodic_spline_second_derivatives(x, y);
  DriverOutput out(output_file);
  const auto cm = out.next_bits("cm", n);
  for (int i = 0; i < n; ++i) EXPECT_BITS(m[i], cm[i]) << " i = " << i;
}

TEST(FortranUbicar, SameIntervals) {
  const auto output_file = reference_file("kernels/ubicar.out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(reference_file("kernels/ubicar.in"));
  const int n = in.integer(), m = in.integer();
  const std::vector<double> x = in.reals(n), t = in.reals(m);
  const auto records = crplib_test::read_driver_output(output_file);
  ASSERT_EQ(static_cast<int>(records.size()), m);
  for (int q = 0; q < m; ++q) {
    std::istringstream fields(records[q].text);
    std::string label;
    int index = 0, klo = 0;
    char outside = ' ';
    fields >> label >> index >> klo >> outside;
    const auto interval = locate_interval(x, t[q]);
    if (outside == 'T') {
      EXPECT_FALSE(interval.has_value()) << t[q];
    } else {
      ASSERT_TRUE(interval.has_value()) << t[q];
      EXPECT_EQ(*interval + 1, klo) << t[q];  // Fortran klo is 1-based
    }
  }
}
