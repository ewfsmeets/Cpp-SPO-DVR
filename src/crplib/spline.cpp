// Cubic splines of the CRP library. Fortran: Spline.f90 (TRIDIA, TRICYCLE, DSPLIN, TRISPL, DPL, DPLP,
// SPLINT, UBICAR, SPLIN2D, SPLCOFII, SPLINII). Formulas keep the Fortran order of operations;
// integer powers x**2, x**3 are written as products (gfortran expands them the same way).
#include "crplib/spline.hpp"

#include <string>

#include "crplib/errors.hpp"

namespace crplib::detail {

// Fortran: TRIDIA. Forward elimination, then back substitution.
void solve_tridiagonal(std::span<double> diagonal, std::span<const double> lower, std::span<const double> upper,
                       std::span<double> rhs, std::span<double> solution) {
  const std::size_t n = diagonal.size();
  for (std::size_t i = 1; i < n; ++i) {
    const double rap = lower[i - 1] / diagonal[i - 1];
    diagonal[i] = diagonal[i] - rap * upper[i - 1];
    rhs[i] = rhs[i] - rap * rhs[i - 1];
  }
  solution[n - 1] = rhs[n - 1] / diagonal[n - 1];
  for (std::size_t i = n - 1; i-- > 0;) solution[i] = (rhs[i] - upper[i] * solution[i + 1]) / diagonal[i];
}

// Fortran: TRICYCLE. The cyclic system as a tridiagonal one plus a rank-one correction
// (Sherman-Morrison with the vector (alp, 0, ..., 0, bottom_left)).
void solve_cyclic_tridiagonal(std::span<double> diagonal, std::span<const double> lower, std::span<const double> upper,
                              double bottom_left, double top_right, std::span<double> rhs, std::span<double> solution) {
  const std::size_t n = diagonal.size();
  const double alp = -diagonal[0];
  diagonal[0] = diagonal[0] - alp;
  diagonal[n - 1] = diagonal[n - 1] - bottom_left * top_right / alp;
  std::vector<double> z(diagonal.begin(), diagonal.end());  // Fortran: z = alpha (saved)
  solve_tridiagonal(diagonal, lower, upper, rhs, solution);

  for (std::size_t i = 0; i < n; ++i) diagonal[i] = z[i];
  rhs[0] = alp;
  rhs[n - 1] = bottom_left;
  for (std::size_t i = 1; i + 1 < n; ++i) rhs[i] = 0.0;
  solve_tridiagonal(diagonal, lower, upper, rhs, z);
  const double fact = (solution[0] + top_right * solution[n - 1] / alp) / (1.0 + z[0] + top_right * z[n - 1] / alp);
  for (std::size_t i = 0; i < n; ++i) solution[i] = solution[i] - fact * z[i];
}

namespace {

// Contribution of one end condition to the tridiagonal system of DSPLIN: right-hand side (s),
// diagonal (fac) and off-diagonal (fab) corrections. h_end is the end interval, h_next the one next
// to it, slope_end the secant slope of the end interval. Fortran: s1, fac1, fab1 (sn, facn, fabn)
struct EndTerms {
  double s, fac, fab;
};

// Fortran: DSPLIN "Contour condition at X(1)"
EndTerms first_end_terms(SplineEnd end, std::span<const double> x, std::span<const double> y) {
  switch (end.condition) {
    case EndCondition::not_a_knot:
      return {0.0, (x[1] - x[0]) * (1.0 + (x[1] - x[0]) / (x[2] - x[1])) / 6.0,
              -((x[1] - x[0]) * (x[1] - x[0])) / (6.0 * (x[2] - x[1]))};
    case EndCondition::second_derivative: return {-end.value * (x[1] - x[0]) / 6.0, 0.0, 0.0};
    case EndCondition::first_derivative:
      return {0.5 * (end.value - (y[1] - y[0]) / (x[1] - x[0])), -(x[1] - x[0]) / 12.0, 0.0};
    case EndCondition::natural: return {0.0, 0.0, 0.0};
  }
  return {0.0, 0.0, 0.0};
}

// Fortran: DSPLIN "Contour condition at X(n)"
EndTerms last_end_terms(SplineEnd end, std::span<const double> x, std::span<const double> y) {
  const std::size_t n = x.size();
  switch (end.condition) {
    case EndCondition::not_a_knot:
      return {0.0, (x[n - 1] - x[n - 2]) * (1.0 + (x[n - 1] - x[n - 2]) / (x[n - 2] - x[n - 3])) / 6.0,
              -((x[n - 1] - x[n - 2]) * (x[n - 1] - x[n - 2])) / (6.0 * (x[n - 2] - x[n - 3]))};
    case EndCondition::second_derivative: return {-end.value * (x[n - 1] - x[n - 2]) / 6.0, 0.0, 0.0};
    case EndCondition::first_derivative:
      return {-0.5 * (end.value - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])), -(x[n - 1] - x[n - 2]) / 12.0, 0.0};
    case EndCondition::natural: return {0.0, 0.0, 0.0};
  }
  return {0.0, 0.0, 0.0};
}

bool strictly_increasing(const std::vector<double>& x) {
  for (std::size_t i = 1; i < x.size(); ++i)
    if (!(x[i] > x[i - 1])) return false;
  return true;
}

[[noreturn]] void not_monotonic(double a, double b) {
  throw DataError("spline abscissae are not monotonic: " + std::to_string(a) + " and " + std::to_string(b));
}

}  // namespace

// Fortran: DSPLIN. Unknowns F''(x_1) .. F''(x_{n-2}); the end values follow from the conditions.
std::vector<double> spline_second_derivatives(std::span<const double> x, std::span<const double> y, SplineEnd first,
                                              SplineEnd last) {
  const std::size_t n = x.size();
  if (n < 4) throw DataError("a cubic spline needs at least 4 points, got " + std::to_string(n));
  const EndTerms start = first_end_terms(first, x, y);
  const EndTerms end = last_end_terms(last, x, y);

  const std::size_t m = n - 2;
  std::vector<double> diagonal(m), upper(m), lower(m), rhs(m);  // Fortran: alpha, beta, cm (copy), b
  double c = x[1] - x[0];
  double e = (y[1] - y[0]) / c;
  for (std::size_t i = 2; i < n; ++i) {
    const std::size_t k = i - 2;
    const double a = c;
    c = x[i] - x[i - 1];
    if (a * c <= 0.0) not_monotonic(x[i - 1], x[i]);
    diagonal[k] = (a + c) / 3.0;
    upper[k] = c / 6.0;
    lower[k] = upper[k];
    const double d = e;
    e = (y[i] - y[i - 1]) / c;
    rhs[k] = e - d;
  }
  rhs[0] = rhs[0] + start.s;
  rhs[m - 1] = rhs[m - 1] + end.s;
  diagonal[0] = diagonal[0] + start.fac;
  diagonal[m - 1] = diagonal[m - 1] + end.fac;
  lower[m - 2] = lower[m - 2] + end.fab;
  upper[0] = upper[0] + start.fab;

  // Fortran: CALL TRIDIA(alpha, cm, beta, b, cm(2), n-2). The Fortran passes cm both as the
  // sub-diagonal and (shifted by one) as the solution; the forward sweep reads the sub-diagonal
  // before anything is written, so a separate copy gives the same result.
  std::vector<double> m2(n);
  solve_tridiagonal(diagonal, lower, upper, rhs, std::span<double>(m2).subspan(1, m));

  const double h0 = x[1] - x[0];
  switch (first.condition) {
    case EndCondition::not_a_knot: m2[0] = m2[1] * (1.0 + h0 / (x[2] - x[1])) - m2[2] * h0 / (x[2] - x[1]); break;
    case EndCondition::second_derivative: m2[0] = first.value; break;
    case EndCondition::first_derivative: m2[0] = -6.0 * start.s / h0 - m2[1] / 2.0; break;
    case EndCondition::natural: m2[0] = 0.0; break;
  }
  const double hn = x[n - 1] - x[n - 2];
  switch (last.condition) {
    case EndCondition::not_a_knot:
      m2[n - 1] = m2[n - 2] * (1.0 + hn / (x[n - 2] - x[n - 3])) - m2[n - 3] * hn / (x[n - 2] - x[n - 3]);
      break;
    case EndCondition::second_derivative: m2[n - 1] = last.value; break;
    case EndCondition::first_derivative: m2[n - 1] = -6.0 * end.s / hn - m2[n - 2] / 2.0; break;
    case EndCondition::natural: m2[n - 1] = 0.0; break;
  }
  return m2;
}

// Fortran: TRISPL. The n-1 unknowns F''(x_1) .. F''(x_{n-1}) of a cyclic system; F''(x_0) = F''(x_{n-1}).
std::vector<double> periodic_spline_second_derivatives(std::span<const double> x, std::span<const double> y) {
  const std::size_t n = x.size();
  if (n < 4) throw DataError("a periodic cubic spline needs at least 4 points, got " + std::to_string(n));
  const std::size_t m = n - 1;
  std::vector<double> diagonal(m), lower(m), upper(m), rhs(m);  // Fortran: alpha, beta, gamma, b
  double c = x[1] - x[0];
  double e = (y[1] - y[0]) / c;
  for (std::size_t i = 2; i < n; ++i) {
    const std::size_t k = i - 2;
    const double a = c;
    c = x[i] - x[i - 1];
    if (a * c <= 0.0) not_monotonic(x[i - 1], x[i]);
    diagonal[k] = (a + c) / 3.0;
    lower[k] = c / 6.0;
    upper[k] = lower[k];
    const double d = e;
    e = (y[i] - y[i - 1]) / c;
    rhs[k] = e - d;
  }
  // The row that closes the period: interval n-2 followed by interval 0.
  const double a = c;
  c = x[1] - x[0];
  diagonal[m - 1] = (a + c) / 3.0;
  lower[m - 1] = c / 6.0;
  upper[m - 1] = lower[m - 1];
  const double d = e;
  e = (y[1] - y[0]) / c;
  rhs[m - 1] = e - d;
  const double corner = lower[m - 1];  // Fortran: cbot = ctop = beta(n-1)

  std::vector<double> m2(n);
  solve_cyclic_tridiagonal(diagonal, lower, upper, corner, corner, rhs, std::span<double>(m2).subspan(1, m));
  m2[0] = m2[n - 1];
  return m2;
}

CubicSpline1D::CubicSpline1D(std::vector<double> x, std::vector<double> y, SplineEnd first, SplineEnd last)
    : x_(std::move(x)), y_(std::move(y)) {
  if (x_.size() != y_.size()) throw DataError("spline: different numbers of abscissae and values");
  m_ = spline_second_derivatives(x_, y_, first, last);
}

namespace {

// The interval index k (1-based, as in the Fortran) used by DPL and DPLP: x(k-1) < t <= x(k) for
// increasing x (with the comparisons reversed for decreasing x); k = 1 below and k = n + 1 above the
// knots; k = 2 for t = x(1). The Fortran finds it by hunting from k = 2 and bisection; the result
// is the end of the prefix of knots j with (t > x(j)) == increasing, found here by bisection.
int dpl_interval(const std::vector<double>& x, double t) {
  const int n = static_cast<int>(x.size());
  const bool increasing = x[1] > x[0];
  int lo = 0;  // number of knots in the prefix
  int hi = n;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if ((t > x[mid]) == increasing) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  int k = lo + 1;
  if (t == x[0]) k = 2;
  return k;
}

}  // namespace

// Fortran: DPL. With k the Fortran's 1-based interval index, the knots are i = k - 1 and i - 1.
double CubicSpline1D::value_dpl(double t) const {
  const int n = size();
  const int k = dpl_interval(x_, t);
  if (k <= 1) {
    const double e = x_[1] - x_[0];
    return ((y_[1] - y_[0]) / e - m_[1] * e / 6.0) * (t - x_[0]) + y_[0];
  }
  if (k > n) {
    const double e = x_[n - 1] - x_[n - 2];
    return ((y_[n - 1] - y_[n - 2]) / e + m_[n - 2] * e / 6.0) * (t - x_[n - 1]) + y_[n - 1];
  }
  const int i = k - 1;
  const double f = x_[i] - t;
  const double g = t - x_[i - 1];
  const double e = f + g;
  return (-g * f * (m_[i - 1] * (f + e) + m_[i] * (g + e)) + 6.0 * (g * y_[i] + f * y_[i - 1])) / (6.0 * e);
}

// Fortran: DPLP
double CubicSpline1D::slope_dplp(double t) const {
  const int n = size();
  const int k = dpl_interval(x_, t);
  if (k <= 1) {
    const double e = x_[1] - x_[0];
    const double g = t - x_[0];
    const double f = x_[1] - t;
    return (y_[1] - y_[0]) / e + (m_[1] * (2.0 * g - f) + m_[0] * (g - 2.0 * f)) / 6.0;
  }
  if (k > n) {
    const double e = x_[n - 1] - x_[n - 2];
    const double g = t - x_[n - 2];
    const double f = x_[n - 1] - t;
    return (y_[n - 1] - y_[n - 2]) / e + (m_[n - 1] * (2.0 * g - f) + m_[n - 2] * (g - 2.0 * f)) / 6.0;
  }
  const int i = k - 1;
  const double f = x_[i] - t;
  const double g = t - x_[i - 1];
  const double e = f + g;
  const double e2 = e * e;
  return (-m_[i - 1] * (3.0 * f * f - e2) + m_[i] * (3.0 * g * g - e2) + 6.0 * (y_[i] - y_[i - 1])) / (6.0 * e);
}

// Fortran: SPLINT
CubicSpline1D::ValueAndSlope CubicSpline1D::value_and_slope_splint(double t) const {
  int klo = 0;
  int khi = size() - 1;
  while (khi - klo > 1) {
    const int k = (khi + klo) / 2;
    if (x_[k] > t) {
      khi = k;
    } else {
      klo = k;
    }
  }
  const double h = x_[khi] - x_[klo];
  const double a = (x_[khi] - t) / h;
  const double b = (t - x_[klo]) / h;
  const double value =
      a * y_[klo] + b * y_[khi] + ((a * a * a - a) * m_[klo] + (b * b * b - b) * m_[khi]) * (h * h) / 6.0;
  const double slope =
      (y_[khi] - y_[klo]) / h + h * (m_[khi] * (3.0 * b * b - 1.0) - m_[klo] * (3.0 * a * a - 1.0)) / 6.0;
  return {value, slope};
}

// Fortran: UBICAR (klo returned 0-based)
std::optional<int> locate_interval(std::span<const double> x, double t) {
  const int n = static_cast<int>(x.size());
  if (t < x[0] || t > x[n - 1]) return std::nullopt;
  int klo = 0;
  int khi = n - 1;
  while (khi - klo > 1) {
    const int k = (khi + klo) / 2;
    if (x[k] > t) {
      khi = k;
    } else {
      klo = k;
    }
  }
  return klo;
}

// Fortran: SPLIN2D. Not-a-knot splines along x for every y, along y for every x, and along x of the
// y second derivatives.
BicubicSpline::BicubicSpline(std::vector<double> x, std::vector<double> y, std::vector<double> f)
    : x_(std::move(x)), y_(std::move(y)), f_(std::move(f)) {
  const std::size_t nx = x_.size();
  const std::size_t ny = y_.size();
  if (nx > 999 || ny > 999)
    throw DataError("a 2D spline grid may have at most 999 points per direction (Fortran SPLIN2D)");
  if (f_.size() != nx * ny) throw DataError("2D spline: the number of values does not match the grid");
  // UBICAR (locate) assumes increasing knots.
  if (!strictly_increasing(x_) || !strictly_increasing(y_))
    throw DataError("2D spline: the grid is not strictly increasing");
  d2f_dx2_.resize(nx * ny);
  d2f_dy2_.resize(nx * ny);
  d4f_dx2dy2_.resize(nx * ny);
  std::vector<double> column(ny);
  for (std::size_t j = 0; j < ny; ++j) {
    const std::vector<double> m =
        spline_second_derivatives(x_, std::span<const double>(f_).subspan(nx * j, nx), not_a_knot, not_a_knot);
    for (std::size_t i = 0; i < nx; ++i) d2f_dx2_[i + nx * j] = m[i];
  }
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) column[j] = f_[i + nx * j];
    const std::vector<double> m = spline_second_derivatives(y_, column, not_a_knot, not_a_knot);
    for (std::size_t j = 0; j < ny; ++j) d2f_dy2_[i + nx * j] = m[j];
  }
  for (std::size_t j = 0; j < ny; ++j) {
    const std::vector<double> m =
        spline_second_derivatives(x_, std::span<const double>(d2f_dy2_).subspan(nx * j, nx), not_a_knot, not_a_knot);
    for (std::size_t i = 0; i < nx; ++i) d4f_dx2dy2_[i + nx * j] = m[i];
  }
}

namespace {

// Fortran: SPLCOFII, one direction: weights of the value and of the derivative at t in [x_i, x_{i+1}].
void cell_weights(std::span<const double> x, int i, double t, CellWeights& value, CellWeights& slope) {
  const double h2s6 = (x[i + 1] - x[i]) * (x[i + 1] - x[i]) / 6.0;  // Fortran: dxia2s6
  value.a = (x[i + 1] - t) / (x[i + 1] - x[i]);
  value.b = 1.0 - value.a;
  value.c = (value.a * value.a * value.a - value.a) * h2s6;
  value.d = (value.b * value.b * value.b - value.b) * h2s6;
  slope.a = -1.0 / (x[i + 1] - x[i]);
  slope.b = -slope.a;
  slope.c = (3.0 * (value.a * value.a) * slope.a - slope.a) * h2s6;
  slope.d = (3.0 * (value.b * value.b) * slope.b - slope.b) * h2s6;
}

}  // namespace

std::optional<GridCell> BicubicSpline::locate(double x, double y) const {
  const std::optional<int> i = locate_interval(x_, x);
  if (!i) return std::nullopt;
  const std::optional<int> j = locate_interval(y_, y);
  if (!j) return std::nullopt;
  GridCell cell{};
  cell.i = *i;
  cell.j = *j;
  cell_weights(x_, cell.i, x, cell.x, cell.x_slope);
  cell_weights(y_, cell.j, y, cell.y, cell.y_slope);
  return cell;
}

// Fortran: SPLINII. vint(1) and vint(3) sum over the x weights outside, vint(2) over the y weights.
BicubicValue BicubicSpline::evaluate_in_cell(const GridCell& cell) const {
  const int i = cell.i;
  const int j = cell.j;
  const std::vector<double>& v = f_;
  const std::vector<double>& fxx = d2f_dx2_;
  const std::vector<double>& fyy = d2f_dy2_;
  const std::vector<double>& fxxyy = d4f_dx2dy2_;
  const int i00 = index(i, j), i10 = index(i + 1, j), i01 = index(i, j + 1), i11 = index(i + 1, j + 1);

  auto x_outer = [&](const CellWeights& wx, const CellWeights& wy) {
    return wx.a * (wy.a * v[i00] + wy.b * v[i01] + wy.c * fyy[i00] + wy.d * fyy[i01]) +
           wx.b * (wy.a * v[i10] + wy.b * v[i11] + wy.c * fyy[i10] + wy.d * fyy[i11]) +
           wx.c * (wy.a * fxx[i00] + wy.b * fxx[i01] + wy.c * fxxyy[i00] + wy.d * fxxyy[i01]) +
           wx.d * (wy.a * fxx[i10] + wy.b * fxx[i11] + wy.c * fxxyy[i10] + wy.d * fxxyy[i11]);
  };
  auto y_outer = [&](const CellWeights& wx, const CellWeights& wy) {
    return wy.a * (wx.a * v[i00] + wx.b * v[i10] + wx.c * fxx[i00] + wx.d * fxx[i10]) +
           wy.b * (wx.a * v[i01] + wx.b * v[i11] + wx.c * fxx[i01] + wx.d * fxx[i11]) +
           wy.c * (wx.a * fyy[i00] + wx.b * fyy[i10] + wx.c * fxxyy[i00] + wx.d * fxxyy[i10]) +
           wy.d * (wx.a * fyy[i01] + wx.b * fyy[i11] + wx.c * fxxyy[i01] + wx.d * fxxyy[i11]);
  };
  return {x_outer(cell.x, cell.y), y_outer(cell.x, cell.y_slope), x_outer(cell.x_slope, cell.y)};
}

std::optional<BicubicValue> BicubicSpline::evaluate(double x, double y) const {
  const std::optional<GridCell> cell = locate(x, y);
  if (!cell) return std::nullopt;
  return evaluate_in_cell(*cell);
}

}  // namespace crplib::detail
