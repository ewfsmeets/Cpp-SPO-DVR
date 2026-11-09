// The interpolation function of the C4v CRP PES at one molecular geometry: angular interpolation at
// the 8 Fourier sites (top, bridges, hollow, t2h and their rotated images) and the Fourier
// interpolation over the molecule's centre. Fortran: the site part of RRINT and SIS8x8
// (CRP_C4v_MoleculeSurface.f90).
#include <array>
#include <cmath>

#include "crplib/constants.hpp"
#include "crplib/molecule_surface_detail.hpp"
#include "crplib/site_interpolation.hpp"

namespace crplib::detail {

namespace {

double component(const SiteValue& s, int i) {
  switch (i) {
    case 0: return s.v;
    case 1: return s.d_z;
    case 2: return s.d_r;
    case 3: return s.d_phi;
    default: return s.d_theta;
  }
}

// The coefficients A1..A8 of SIS8x8 for component i of the sites. Psite index k + 1 = s[k]; the
// coefficients are the inverse of the basis at the sites top (0,0), bridge (1/2,0), hollow (1/2,1/2),
// t2h (1/4,1/4), bridge (0,1/2) and t2h (3/4,1/4), (3/4,3/4), (1/4,3/4) (Cartesian, lattice units).
std::array<double, 8> sis8x8_coefficients(const std::array<SiteValue, 8>& s, int i) {
  const double p1 = component(s[0], i), p2 = component(s[1], i), p3 = component(s[2], i), p4 = component(s[3], i);
  const double p5 = component(s[4], i), p6 = component(s[5], i), p7 = component(s[6], i), p8 = component(s[7], i);
  return {1.0 / 8.0 * (p1 + p2 + p5 + p3 + p4 + p6 + p7 + p8),
          1.0 / 4.0 * (p1 - p2 + p5 - p3),
          1.0 / 4.0 * (p1 + p2 - p5 - p3),
          1.0 / 4.0 * (p4 - p6 - p7 + p8),
          1.0 / 4.0 * (p4 + p6 - p7 - p8),
          1.0 / 8.0 * (p1 - p2 - p5 + p3 - p4 + p6 - p7 + p8),
          1.0 / 8.0 * (p1 - p2 - p5 + p3 + p4 - p6 + p7 - p8),
          1.0 / 16.0 * (p1 + p2 + p5 + p3 - p4 - p6 - p7 - p8)};
}

// Fortran: SIS8x8. Want::value: Vinter(1) only.
Interaction sis8x8(const std::array<SiteValue, 8>& s, double xcm, double ycm, double lattice_constant, Want want) {
  const double fac = 2.0 * pi / lattice_constant;
  const double f1 = 1.0;
  const double f2 = std::cos(fac * xcm);
  const double f3 = std::cos(fac * ycm);
  const double f4 = std::sin(fac * xcm);
  const double f5 = std::sin(fac * ycm);
  const double f6 = std::cos(fac * (xcm + ycm));
  const double f7 = std::cos(fac * (xcm - ycm));
  const double f8 = std::cos(fac * (2.0 * xcm)) + std::cos(fac * (2.0 * ycm));

  Interaction v{};
  {
    const auto [a1, a2, a3, a4, a5, a6, a7, a8] = sis8x8_coefficients(s, 0);
    v[0] = a1 * f1 + a2 * f2 + a3 * f3 + a4 * f4 + a5 * f5 + a6 * f6 + a7 * f7 + a8 * f8;
  }
  if (want == Want::value) return v;

  const double df1dx = 0.0, df1dy = 0.0;
  const double df2dx = -fac * std::sin(fac * xcm), df2dy = 0.0;
  const double df3dx = 0.0, df3dy = -fac * std::sin(fac * ycm);
  const double df4dx = fac * std::cos(fac * xcm), df4dy = 0.0;
  const double df5dx = 0.0, df5dy = fac * std::cos(fac * ycm);
  const double df6dx = -fac * std::sin(fac * (xcm + ycm)), df6dy = -fac * std::sin(fac * (xcm + ycm));
  const double df7dx = -fac * std::sin(fac * (xcm - ycm)), df7dy = fac * std::sin(fac * (xcm - ycm));
  const double df8dx = -2.0 * fac * std::sin(fac * (2.0 * xcm)), df8dy = -2.0 * fac * std::sin(fac * (2.0 * ycm));
  {
    const auto [a1, a2, a3, a4, a5, a6, a7, a8] = sis8x8_coefficients(s, 0);
    v[5] = a1 * df1dx + a2 * df2dx + a3 * df3dx + a4 * df4dx + a5 * df5dx + a6 * df6dx + a7 * df7dx + a8 * df8dx;
    v[6] = a1 * df1dy + a2 * df2dy + a3 * df3dy + a4 * df4dy + a5 * df5dy + a6 * df6dy + a7 * df7dy + a8 * df8dy;
  }
  for (int i = 1; i < 5; ++i) {
    const auto [a1, a2, a3, a4, a5, a6, a7, a8] = sis8x8_coefficients(s, i);
    v[i] = a1 * f1 + a2 * f2 + a3 * f3 + a4 * f4 + a5 * f5 + a6 * f6 + a7 * f7 + a8 * f8;
  }
  return v;
}

}  // namespace

// Fortran: RRINT, the part between the cut evaluation and SIS8x8.
Interaction square_interaction(std::span<const BicubicValue> cuts, double theta, double phi, double x, double y,
                               double lattice_constant, Want want) {
  const ThetaTrig t = theta_trig(theta, want);
  std::array<SiteValue, 8> s{};
  s[0] = square_top(cuts.subspan(0, 3), t, phi, want);
  const auto bridge = cuts.subspan(3, 4);
  s[1] = square_bridge(bridge, t, phi, want);
  s[4] = square_bridge(bridge, t, phi + pi / 2.0, want);
  s[2] = square_hollow(cuts.subspan(7, 3), t, phi, want);
  const auto t2h = cuts.subspan(10, 6);
  s[3] = square_t2h(t2h, t, phi, want);
  s[5] = square_t2h(t2h, t, phi - pi / 2.0, want);
  s[6] = square_t2h(t2h, t, phi - pi, want);
  s[7] = square_t2h(t2h, t, phi + pi / 2.0, want);
  return sis8x8(s, x, y, lattice_constant, want);
}

}  // namespace crplib::detail
