#pragma once
// Cubic spline interpolation of the CRP library: tridiagonal solvers, 1D splines with the two
// evaluation formulas the Fortran uses, the periodic 1D spline and the bicubic 2D spline.
// Fortran: Spline.f90 (A. Salin, H. F. Busnengo). Every formula keeps the Fortran's order of
// operations (results bit for bit equal to the Fortran library).

#include <optional>
#include <span>
#include <vector>

namespace crplib::detail {

// Solves the tridiagonal system with diagonal `diagonal`, sub-diagonal `lower` (lower[i - 1] in
// row i) and super-diagonal `upper` (upper[i] in row i). `diagonal` and `rhs` are overwritten.
// Fortran: TRIDIA
void solve_tridiagonal(std::span<double> diagonal, std::span<const double> lower, std::span<const double> upper,
                       std::span<double> rhs, std::span<double> solution);

// As solve_tridiagonal, with the extra corner elements `bottom_left` (row n, column 1) and
// `top_right` (row 1, column n) of a cyclic system. `diagonal` and `rhs` are overwritten.
// Fortran: TRICYCLE
void solve_cyclic_tridiagonal(std::span<double> diagonal, std::span<const double> lower, std::span<const double> upper,
                              double bottom_left, double top_right, std::span<double> rhs, std::span<double> solution);

// Condition of a cubic spline at one end. Fortran: the codes ic1 / icn of DSPLIN.
enum class EndCondition {
  not_a_knot = 0,         // the same cubic in the first (last) two intervals
  second_derivative = 1,  // F'' given
  first_derivative = 2,   // F' given (clamped spline)
  natural = 3,            // F'' = 0 (Fortran: any code > 2)
};
struct SplineEnd {
  EndCondition condition;
  double value;  // the given derivative (ignored for not_a_knot and natural)
};
inline constexpr SplineEnd not_a_knot{EndCondition::not_a_knot, 0.0};

// The second derivatives F''(x_i) of the cubic spline through (x_i, y_i). x must be monotonic
// (increasing or decreasing) with at least 4 points; throws DataError otherwise.
// Fortran: DSPLIN(n, x, y, cm, cm1, ic1, cmn, icn)
std::vector<double> spline_second_derivatives(std::span<const double> x, std::span<const double> y, SplineEnd first,
                                              SplineEnd last);

// The second derivatives of the periodic cubic spline through (x_i, y_i), i = 0..n-1, with period
// x_{n-1} - x_0 (y_{n-1} = y_0 is assumed). Fortran: TRISPL
std::vector<double> periodic_spline_second_derivatives(std::span<const double> x, std::span<const double> y);

// A cubic spline in one variable with its second derivatives at the knots.
class CubicSpline1D {
 public:
  CubicSpline1D() = default;
  // Fortran: Allocate_Spline1D + Compute_Spline1D (DSPLIN)
  CubicSpline1D(std::vector<double> x, std::vector<double> y, SplineEnd first, SplineEnd last);

  // F(t). Outside [x_0, x_{n-1}] linear extrapolation from the first (last) interval.
  // Fortran: DPL (Evaluate_Spline1D)
  [[nodiscard]] double value_dpl(double t) const;
  // F'(t); outside the knots the derivative of the end interval's cubic. Fortran: DPLP
  // (Evaluate_Spline1D_Derivative)
  [[nodiscard]] double slope_dplp(double t) const;

  struct ValueAndSlope {
    double value;
    double slope;
  };
  // F(t) and F'(t) with a different formula (bisection; the end cubics continued outside).
  // Fortran: SPLINT(xa, ya, cya, n, x, y, yp)
  [[nodiscard]] ValueAndSlope value_and_slope_splint(double t) const;

  [[nodiscard]] const std::vector<double>& x() const { return x_; }
  [[nodiscard]] const std::vector<double>& y() const { return y_; }
  [[nodiscard]] const std::vector<double>& second_derivatives() const { return m_; }
  [[nodiscard]] int size() const { return static_cast<int>(x_.size()); }

 private:
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> m_;  // F''(x_i). Fortran: Coefs
};

// The index klo (0-based) of the interval [x_klo, x_klo+1] that contains t, for increasing x;
// nullopt if t < x_0 or t > x_{n-1}. t = x_{n-1} gives the last interval. Fortran: UBICAR
std::optional<int> locate_interval(std::span<const double> x, double t);

// Interpolation weights of one cell in one direction. Value: a, b (of the function at the two
// knots), c, d (of the second derivatives). Fortran: ax0, bx0, cx0, dx0 (and ax0p ... for the
// derivative weights).
struct CellWeights {
  double a, b, c, d;
};

// A point located in a grid cell: the cell's lower-left knot (i, j) and the weights for the value
// and the derivative along x and along y. Fortran: SPLCOFII (TYPE UBIType)
struct GridCell {
  int i, j;
  CellWeights x, y, x_slope, y_slope;
};

// Value and first derivatives of a bicubic spline. Fortran: SPLINII vint(1), vint(2) (d/dy),
// vint(3) (d/dx)
struct BicubicValue {
  double value;
  double d_dy;
  double d_dx;
};

// A bicubic spline on a rectangular grid (both directions not-a-knot). Fortran: Spline2D,
// Compute_Spline2D (SPLIN2D), Evaluate_Spline2D_Alt (UBICAR, SPLCOFII, SPLINII)
class BicubicSpline {
 public:
  BicubicSpline() = default;
  // f(i, j) = f[i + nx * j] at (x_i, y_j). Throws DataError for fewer than 4 or more than 999
  // points in a direction (Fortran: SPLIN2D idn) or non-monotonic grids.
  BicubicSpline(std::vector<double> x, std::vector<double> y, std::vector<double> f);

  // The cell of (x, y); nullopt outside the grid (Fortran: switch = .TRUE.). Fortran: UBICAR, SPLCOFII
  [[nodiscard]] std::optional<GridCell> locate(double x, double y) const;
  // Value and derivatives in a located cell. Fortran: SPLINII
  [[nodiscard]] BicubicValue evaluate_in_cell(const GridCell& cell) const;
  // Both steps. Fortran: Evaluate_Spline2D_Alt
  [[nodiscard]] std::optional<BicubicValue> evaluate(double x, double y) const;

  [[nodiscard]] const std::vector<double>& x() const { return x_; }
  [[nodiscard]] const std::vector<double>& y() const { return y_; }
  [[nodiscard]] const std::vector<double>& f() const { return f_; }
  // The coefficient arrays (same layout as f). Fortran: CoefsX, CoefsY, CoefsXY
  [[nodiscard]] const std::vector<double>& d2f_dx2() const { return d2f_dx2_; }
  [[nodiscard]] const std::vector<double>& d2f_dy2() const { return d2f_dy2_; }
  [[nodiscard]] const std::vector<double>& d4f_dx2dy2() const { return d4f_dx2dy2_; }

 private:
  [[nodiscard]] int index(int i, int j) const { return i + static_cast<int>(x_.size()) * j; }

  std::vector<double> x_, y_, f_;
  std::vector<double> d2f_dx2_, d2f_dy2_, d4f_dx2dy2_;
};

}  // namespace crplib::detail
