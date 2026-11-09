// 3D spline periodic in x, y and clamped in z (C4v atom surface). Fortran: Spline.f90 SPLINE3D,
// int_PES6CuSPL2DPER, SPLINT3D. The Fortran's x**3.0 is a call of pow (6 sites in SPLINT3D) and is
// kept as std::pow; x**2.0 is a product.
#include "crplib/periodic_spline_3d.hpp"

#include <cmath>
#include <span>
#include <string>

#include "crplib/errors.hpp"
#include "crplib/spline.hpp"

namespace crplib::detail {

namespace {

// Second derivatives of a periodic 2D function g(i, j) = g[i + nx * j]: along y, along x, and along
// x of the y derivatives. Fortran: int_PES6CuSPL2DPER(X, L, Y, M, F, DFDX, DFDY, DFDXDY)
void periodic_2d(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& g,
                 std::vector<double>& gxx, std::vector<double>& gyy, std::vector<double>& gxxyy) {
  const std::size_t nx = x.size();
  const std::size_t ny = y.size();
  gxx.resize(nx * ny);
  gyy.resize(nx * ny);
  gxxyy.resize(nx * ny);
  std::vector<double> line(ny);
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) line[j] = g[i + nx * j];
    const std::vector<double> m = periodic_spline_second_derivatives(y, line);
    for (std::size_t j = 0; j < ny; ++j) gyy[i + nx * j] = m[j];
  }
  for (std::size_t j = 0; j < ny; ++j) {
    const std::vector<double> m = periodic_spline_second_derivatives(x, std::span<const double>(g).subspan(nx * j, nx));
    for (std::size_t i = 0; i < nx; ++i) gxx[i + nx * j] = m[i];
  }
  for (std::size_t j = 0; j < ny; ++j) {
    const std::vector<double> m =
        periodic_spline_second_derivatives(x, std::span<const double>(gyy).subspan(nx * j, nx));
    for (std::size_t i = 0; i < nx; ++i) gxxyy[i + nx * j] = m[i];
  }
}

// Fortran SPLINT3D, one direction: value weights (AX0, BX0, CX0, DX0) and derivative weights
// (AX0P ...) at t in [x_i, x_{i+1}].
struct DirectionWeights {
  CellWeights value;
  CellWeights slope;
};

DirectionWeights weights(const std::vector<double>& x, int i, double t) {
  const double h = x[i + 1] - x[i];  // Fortran: DEX
  DirectionWeights w{};
  w.value.a = (x[i + 1] - t) / h;
  w.value.b = 1.0 - w.value.a;
  w.value.c = (std::pow(w.value.a, 3.0) - w.value.a) / 6.0 * (h * h);
  w.value.d = (std::pow(w.value.b, 3.0) - w.value.b) / 6.0 * (h * h);
  w.slope.a = -1.0 / h;
  w.slope.b = -w.slope.a;
  w.slope.c = (3.0 * (w.value.a * w.value.a) * w.slope.a - w.slope.a) / 6.0 * (h * h);
  w.slope.d = (3.0 * (w.value.b * w.value.b) * w.slope.b - w.slope.b) / 6.0 * (h * h);
  return w;
}

}  // namespace

// Fortran: SPLINE3D. Clamped splines along z for every (x, y); then the periodic 2D second
// derivatives of F and of F_zz in every z plane.
PeriodicSpline3D::PeriodicSpline3D(std::vector<double> x, std::vector<double> y, std::vector<double> z,
                                   std::vector<double> f, const std::vector<double>& dfdz_first,
                                   const std::vector<double>& dfdz_last)
    : x_(std::move(x)), y_(std::move(y)), z_(std::move(z)), f_(std::move(f)) {
  const int nx = static_cast<int>(x_.size());
  const int ny = static_cast<int>(y_.size());
  const int nz = static_cast<int>(z_.size());
  const std::size_t n = f_.size();
  if (n != x_.size() * y_.size() * z_.size() || dfdz_first.size() != x_.size() * y_.size() ||
      dfdz_last.size() != dfdz_first.size())
    throw DataError("3D spline: the number of values does not match the grid");
  c_.zz.resize(n);
  std::vector<double> line(nz);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) line[k] = f_[index(i, j, k)];
      const SplineEnd first{EndCondition::first_derivative, dfdz_first[i + nx * j]};
      const SplineEnd last{EndCondition::first_derivative, dfdz_last[i + nx * j]};
      const std::vector<double> m = spline_second_derivatives(z_, line, first, last);
      for (int k = 0; k < nz; ++k) c_.zz[index(i, j, k)] = m[k];
    }
  }
  c_.xx.resize(n);
  c_.yy.resize(n);
  c_.xxyy.resize(n);
  c_.zzxx.resize(n);
  c_.zzyy.resize(n);
  c_.zzxxyy.resize(n);
  std::vector<double> plane(nx * ny), plane_zz(nx * ny);
  std::vector<double> xx, yy, xxyy, zzxx, zzyy, zzxxyy;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        plane[i + nx * j] = f_[index(i, j, k)];
        plane_zz[i + nx * j] = c_.zz[index(i, j, k)];
      }
    }
    periodic_2d(x_, y_, plane, xx, yy, xxyy);
    periodic_2d(x_, y_, plane_zz, zzxx, zzyy, zzxxyy);
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        const int at = index(i, j, k);
        c_.xx[at] = xx[i + nx * j];
        c_.yy[at] = yy[i + nx * j];
        c_.xxyy[at] = xxyy[i + nx * j];
        c_.zzxx[at] = zzxx[i + nx * j];
        c_.zzyy[at] = zzyy[i + nx * j];
        c_.zzxxyy[at] = zzxxyy[i + nx * j];
      }
    }
  }
}

// Fortran: SPLINT3D
PeriodicSpline3D::Value PeriodicSpline3D::evaluate(double x, double y, double z) const {
  const std::optional<int> ia = locate_interval(x_, x);
  const std::optional<int> ja = locate_interval(y_, y);
  const std::optional<int> ka = locate_interval(z_, z);
  if (!ia || !ja || !ka) {
    throw EvaluationError("C4v atom-surface spline: point (" + std::to_string(x) + ", " + std::to_string(y) + ", " +
                          std::to_string(z) + ") outside the interpolation grid (Fortran SPLINT3D: out of range)");
  }
  const DirectionWeights wx = weights(x_, *ia, x);
  const DirectionWeights wy = weights(y_, *ja, y);
  const DirectionWeights wz = weights(z_, *ka, z);

  // One z plane k: Σ over the y weights outside, the x weights inside (Fortran ATERM, BTERM, ...),
  // for the function array g with its x, y and xy second derivatives.
  auto plane_sum = [this](const CellWeights& ex, const CellWeights& ey, const std::vector<double>& g,
                          const std::vector<double>& gxx, const std::vector<double>& gyy,
                          const std::vector<double>& gxxyy, int i, int j, int k) {
    const int i00 = index(i, j, k), i10 = index(i + 1, j, k), i01 = index(i, j + 1, k), i11 = index(i + 1, j + 1, k);
    return ey.a * (ex.a * g[i00] + ex.b * g[i10] + ex.c * gxx[i00] + ex.d * gxx[i10]) +
           ey.b * (ex.a * g[i01] + ex.b * g[i11] + ex.c * gxx[i01] + ex.d * gxx[i11]) +
           ey.c * (ex.a * gyy[i00] + ex.b * gyy[i10] + ex.c * gxxyy[i00] + ex.d * gxxyy[i10]) +
           ey.d * (ex.a * gyy[i01] + ex.b * gyy[i11] + ex.c * gxxyy[i01] + ex.d * gxxyy[i11]);
  };
  // The four plane terms of F (planes ka, ka+1) and of F_zz (planes ka, ka+1).
  struct Terms {
    double a, b, c, d;
  };
  auto terms = [&](const CellWeights& ex, const CellWeights& ey) {
    return Terms{plane_sum(ex, ey, f_, c_.xx, c_.yy, c_.xxyy, *ia, *ja, *ka),
                 plane_sum(ex, ey, f_, c_.xx, c_.yy, c_.xxyy, *ia, *ja, *ka + 1),
                 plane_sum(ex, ey, c_.zz, c_.zzxx, c_.zzyy, c_.zzxxyy, *ia, *ja, *ka),
                 plane_sum(ex, ey, c_.zz, c_.zzxx, c_.zzyy, c_.zzxxyy, *ia, *ja, *ka + 1)};
  };
  auto combine = [](const CellWeights& ez, const Terms& t) {
    return ez.a * t.a + ez.b * t.b + ez.c * t.c + ez.d * t.d;
  };

  const Terms base = terms(wx.value, wy.value);  // Fortran: ATERM, BTERM, CTERM, DTERM
  Value result{};
  result.value = combine(wz.value, base);                          // Fortran: F0
  result.d_dx = combine(wz.value, terms(wx.slope, wy.value));      // Fortran: DX (ATERMX ...)
  result.d_dy = combine(wz.value, terms(wx.value, wy.slope));      // Fortran: DY (ATERMY ...)
  result.d_dz = combine(wz.slope, base);                           // Fortran: DZ
  return result;
}

}  // namespace crplib::detail
