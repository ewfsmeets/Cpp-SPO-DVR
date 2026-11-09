// The interpolation function of the C3v and C6v CRP PES at one molecular geometry: angular
// interpolation at the 12 Fourier sites (each site's own cuts, rotated for its symmetry images),
// then the Fourier interpolation over the molecule's centre. Fortran: the site part of RRINTII,
// SISII and FOUR6II (CRP_C6v_MoleculeSurface.f90, CRP_C3v_MoleculeSurface.f90); the site and cut
// tables are Appendix A of docs/crp_library.md.
#include <array>
#include <cmath>

#include "crplib/constants.hpp"
#include "crplib/fourier_tables.hpp"
#include "crplib/molecule_surface_detail.hpp"
#include "crplib/site_interpolation.hpp"

namespace crplib::detail {

// Fortran: FOUR6II. Basis in crystal coordinates (2 pi u, 2 pi v); the y derivative is transformed
// to Cartesian at the end.
void hexagonal_molecule_basis(double x, double y, double lattice_constant, std::span<double> f, std::span<double> dfdx,
                              std::span<double> dfdy, Want want) {
  const double fac = 2.00 * pi / lattice_constant;
  const double bigx = x - y * usq3;
  const double bigy = 2.00 * y * usq3;
  const double fx = fac * bigx;  // Fortran: fBigX
  const double fy = fac * bigy;  // Fortran: fBigY

  const double c2x = std::cos(2.00 * fx), c2y = std::cos(2.00 * fy), c2xy = std::cos(2.00 * (fx + fy));
  const double s2x = std::sin(2.00 * fx), s2y = std::sin(2.00 * fy), s2xy = std::sin(2.00 * (fx + fy));

  f[0] = 1.00;
  f[1] = std::cos(fx);
  f[2] = std::cos(fy);
  f[3] = std::cos(fx + fy);
  f[4] = std::sin(fx);
  f[5] = std::sin(fy);
  f[6] = std::sin(fx + fy);
  f[7] = std::cos(2.00 * fx + fy);
  f[8] = std::cos(fx + 2.00 * fy);
  f[9] = std::cos(fx - fy);
  f[10] = c2x + c2y + c2xy;
  f[11] = -s2x - s2y + s2xy;
  if (want == Want::value) return;

  // Derivatives with respect to the crystal coordinates
  dfdx[0] = 0.00;
  dfdx[1] = -fac * f[4];
  dfdx[2] = 0.00;
  dfdx[3] = -fac * f[6];
  dfdx[4] = fac * f[1];
  dfdx[5] = 0.00;
  dfdx[6] = fac * f[3];
  dfdx[7] = -2.00 * fac * std::sin(2.00 * fx + fy);
  dfdx[8] = -fac * std::sin(fx + 2.00 * fy);
  dfdx[9] = -fac * std::sin(fx - fy);
  dfdx[10] = -2.0 * fac * (s2x + s2xy);
  dfdx[11] = 2.0 * fac * (-c2x + c2xy);

  dfdy[0] = 0.00;
  dfdy[1] = 0.00;
  dfdy[2] = -fac * f[5];
  dfdy[3] = -fac * f[6];
  dfdy[4] = 0.00;
  dfdy[5] = fac * f[2];
  dfdy[6] = fac * f[3];
  dfdy[7] = 0.500 * dfdx[7];
  dfdy[8] = 2.00 * dfdx[8];
  dfdy[9] = -dfdx[9];
  dfdy[10] = -2.0 * fac * (s2y + s2xy);
  dfdy[11] = 2.0 * fac * (-c2y + c2xy);

  // Derivatives with respect to Cartesian coordinates
  for (int i = 0; i < 12; ++i) dfdy[i] = (-dfdx[i] + 2.00 * dfdy[i]) * usq3;
}

namespace {

// Component i of a site value (Fortran: Psite(i + 1, ·)).
double component(const SiteValue& s, int i) {
  switch (i) {
    case 0: return s.v;
    case 1: return s.d_z;
    case 2: return s.d_r;
    case 3: return s.d_phi;
    default: return s.d_theta;
  }
}

// Fortran: SISII. a(j) = sum over the sites of fouri6D(ksit, j) Psite(i, ksit), then sum over j.
// Want::value: Vinter(1) only.
Interaction fourier_interpolation(const std::array<SiteValue, 12>& sites, double x, double y, double lattice_constant,
                                  Want want) {
  std::array<double, 12> f{}, dfdx{}, dfdy{};
  hexagonal_molecule_basis(x, y, lattice_constant, f, dfdx, dfdy, want);
  const std::array<double, 144>& table = fourier_table_molecule_hexagonal();
  Interaction v{};
  v[0] = 0.00;
  v[5] = 0.00;
  v[6] = 0.00;
  for (int j = 0; j < 12; ++j) {
    double a = 0.00;
    for (int ksit = 0; ksit < 12; ++ksit) a = a + table[ksit + 12 * j] * sites[ksit].v;
    v[0] = v[0] + a * f[j];
    v[5] = v[5] + a * dfdx[j];
    v[6] = v[6] + a * dfdy[j];
  }
  if (want == Want::value) return v;
  // First derivatives with respect to Z, r, phi, theta
  for (int i = 1; i < 5; ++i) {
    v[i] = 0.00;
    for (int j = 0; j < 12; ++j) {
      double a = 0.00;
      for (int ksit = 0; ksit < 12; ++ksit) a = a + table[ksit + 12 * j] * component(sites[ksit], i);
      v[i] = v[i] + a * f[j];
    }
  }
  return v;
}

}  // namespace

// Fortran: RRINTII, the part between the cut evaluation and SISII. Fourier site k (Psite index
// k + 1) takes the values of a cut site at an azimuth rotated to that site's frame.
Interaction hexagonal_interaction(Symmetry symmetry, std::span<const BicubicValue> cuts, double theta, double phi,
                                  double x, double y, double lattice_constant, Want want) {
  constexpr double pis3 = pi / 3.00;
  constexpr double pi2s3 = 2.00 * pis3;
  const ThetaTrig t = theta_trig(theta, want);
  std::array<SiteValue, 12> s{};

  s[0] = hexagonal_top(cuts.subspan(0, 3), t, phi, want);  // top
  const auto bridge = cuts.subspan(3, 4);
  s[1] = hexagonal_bridge(bridge, t, phi, want);
  s[4] = hexagonal_bridge(bridge, t, phi + pis3, want);
  s[5] = hexagonal_bridge(bridge, t, phi - pis3, want);
  const auto hollow = cuts.subspan(7, 5);
  if (symmetry == Symmetry::c6v) {
    s[2] = hexagonal_hollow(hollow, t, phi, want);
    s[8] = hexagonal_hollow(hollow, t, phi + pi, want);
    const auto t2h = cuts.subspan(12, 6);
    s[3] = c6v_t2h(t2h, t, phi, want);
    s[6] = c6v_t2h(t2h, t, phi - pi2s3, want);
    s[7] = c6v_t2h(t2h, t, phi + pi2s3, want);
    s[9] = c6v_t2h(t2h, t, phi + pi, want);
    s[10] = c6v_t2h(t2h, t, phi + pis3, want);
    s[11] = c6v_t2h(t2h, t, phi - pis3, want);
  } else {
    s[2] = hexagonal_hollow(hollow, t, phi, want);  // hcp
    const auto t2h = cuts.subspan(12, 6);
    s[3] = c3v_t2h(t2h, t, phi, want);
    s[6] = c3v_t2h(t2h, t, phi - pi2s3, want);
    s[7] = c3v_t2h(t2h, t, phi + pi2s3, want);
    const auto t2f = cuts.subspan(18, 6);
    s[10] = c3v_t2f(t2f, t, phi, want);
    s[9] = c3v_t2f(t2f, t, phi + pi2s3, want);
    s[11] = c3v_t2f(t2f, t, phi - pi2s3, want);
    s[8] = hexagonal_hollow(cuts.subspan(24, 5), t, phi, want);  // fcc
  }
  return fourier_interpolation(s, x, y, lattice_constant, want);
}

}  // namespace crplib::detail
