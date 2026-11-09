// Angular interpolation at one surface site. Fortran: the InterpolateOn* routines of
// CRP_C6v_MoleculeSurface.f90, CRP_C3v_MoleculeSurface.f90 and CRP_C4v_MoleculeSurface.f90. Each
// site function evaluates its formula for the cut values, their Z derivatives and their r
// derivatives (Fortran: the loop i = 1, 3 over Vint(i, k)), and the angular derivatives for the
// values.
#include "crplib/site_interpolation.hpp"

#include <cmath>

#include "crplib/constants.hpp"

namespace crplib::detail {

namespace {

// Component i of cut k (0-based): 0 value, 1 d/dZ, 2 d/dr. Fortran: Vint(i + 1, k + 1)
double vint(std::span<const BicubicValue> cuts, int i, int k) {
  const BicubicValue& c = cuts[k];
  return i == 0 ? c.value : (i == 1 ? c.d_dy : c.d_dx);
}

// Fortran: Psite(i + 1) = x for i = 0, 1, 2
void set(SiteValue& site, int i, double x) {
  if (i == 0) site.v = x;
  if (i == 1) site.d_z = x;
  if (i == 2) site.d_r = x;
}

// The value from the orientations theta = 0, 45, 90, 135 degrees (hexagonal hollow and t2 sites).
double four_orientations(double v0, double v45, double v90, double v135, const ThetaTrig& t) {
  return 0.2500 * ((v0 + v45 + v90 + v135) + (v0 + v90 - v45 - v135) * t.cos4) +
         0.500 * ((v0 - v90) * t.cos2 + (v45 - v135) * t.sin2);
}
double four_orientations_theta_slope(double v0, double v45, double v90, double v135, const ThetaTrig& t) {
  return -(v0 - v90) * t.sin2 + (v45 - v135) * t.cos2 - (v0 + v90 - v45 - v135) * t.sin4;
}
// The phi derivative of the t2 sites (V0 does not depend on phi).
double four_orientations_phi_slope(double dv45, double dv90, double dv135, const ThetaTrig& t) {
  return 0.2500 * ((dv45 + dv90 + dv135) + (dv90 - dv45 - dv135) * t.cos4) +
         0.500 * ((-dv90) * t.cos2 + (dv45 - dv135) * t.sin2);
}

// cos/sin of phi' and 2 phi', with phi' the azimuth shifted to the t2 site's reference direction
// (Fortran: CosPhi, SinPhi, Cos2Phi, Sin2Phi of Phip6). The three hexagonal t2 site functions differ
// only in which cuts play which role; each is written out for the exact Fortran operation order.
struct PhiTrig {
  double cos1, sin1, cos2, sin2;  // of phi' and 2 phi'
};
PhiTrig phi_trig(double phi_shifted) {
  return {std::cos(phi_shifted), std::sin(phi_shifted), std::cos(2.00 * phi_shifted), std::sin(2.00 * phi_shifted)};
}

}  // namespace

ThetaTrig theta_trig(double theta) {
  return {std::cos(2.00 * theta), std::sin(2.00 * theta), std::cos(4.00 * theta), std::sin(4.00 * theta)};
}

// Fortran: InterpolateOnTop. Cuts: theta = 0; 90/0; 90/30.
SiteValue hexagonal_top(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  const double cos6phi = std::cos(6.00 * phi), sin6phi = std::sin(6.00 * phi);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0);
    const double v90 =
        0.500 * ((vint(cuts, i, 1) + vint(cuts, i, 2)) + (vint(cuts, i, 1) - vint(cuts, i, 2)) * cos6phi);
    set(p, i, 0.500 * ((v0 + v90) + (v0 - v90) * t.cos2));
    if (i == 0) {
      const double dv90dphi = -3.00 * (vint(cuts, 0, 1) - vint(cuts, 0, 2)) * sin6phi;
      p.d_phi = 0.500 * dv90dphi * (1.00 - t.cos2);
      p.d_theta = -(v0 - v90) * t.sin2;
    }
  }
  return p;
}

// Fortran: InterpolateOnBridge. Cuts: theta = 0; 90/0; 90/60; 90/90.
SiteValue hexagonal_bridge(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  constexpr double u3 = 1.00 / 3.00;
  const double cos2phi = std::cos(2.00 * phi), sin2phi = std::sin(2.00 * phi);
  const double cos4phi = std::cos(4.00 * phi), sin4phi = std::sin(4.00 * phi);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2), c3 = vint(cuts, i, 3);
    const double v90 =
        (c1 + 2.00 * c2) * u3 - 0.500 * (c3 - c1) * cos2phi + (0.500 * (c3 - c1) - 2.00 * (c2 - c1) * u3) * cos4phi;
    set(p, i, 0.500 * (v0 * (1.00 + t.cos2) + v90 * (1.00 - t.cos2)));
    if (i == 0) {
      const double dv90dphi = (c3 - c1) * sin2phi - (2.00 * (c3 - c1) - 8.00 * (c2 - c1) * u3) * sin4phi;
      p.d_phi = 0.500 * dv90dphi * (1.00 - t.cos2);
      p.d_theta = (v90 - v0) * t.sin2;
    }
  }
  return p;
}

// Fortran: InterpolateOnFCC (C6v), InterpolateOnHCP (C3v). Cuts: theta = 0; 90/0; 90/30; 45/30; 45/210.
SiteValue hexagonal_hollow(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  const double cos3phi = std::cos(3.00 * phi), sin3phi = std::sin(3.00 * phi);
  const double cos6phi = std::cos(6.00 * phi), sin6phi = std::sin(6.00 * phi);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double c3 = vint(cuts, i, 3), c4 = vint(cuts, i, 4);
    const double v45 = 0.500 * ((c3 + c4) + (c3 - c4) * sin3phi);
    const double v135 = 0.500 * ((c3 + c4) - (c3 - c4) * sin3phi);
    const double v90 = 0.500 * ((c1 + c2) + (c1 - c2) * cos6phi);
    set(p, i, four_orientations(v0, v45, v90, v135, t));
    if (i == 0) {
      const double dv45dphi = 1.500 * (c3 - c4) * cos3phi;
      const double dv90dphi = -3.00 * (c1 - c2) * sin6phi;
      p.d_phi = dv90dphi * (0.2500 * (1.00 + t.cos4) - 0.500 * t.cos2) + dv45dphi * t.sin2;
      p.d_theta = four_orientations_theta_slope(v0, v45, v90, v135, t);
    }
  }
  return p;
}

// Fortran: C6v InterpolateOnT2F (the t2h site of C6v). Cuts: theta = 0; 90/30; 90/120; 45/30; 45/210; 45/120.
SiteValue c6v_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  constexpr double pis6 = pi / 6.00;
  const PhiTrig f = phi_trig(phi - pis6);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double c3 = vint(cuts, i, 3), c4 = vint(cuts, i, 4), c5 = vint(cuts, i, 5);
    const double v45 = 0.2500 * ((c3 + 2.00 * c5 + c4) + (c3 - 2.00 * c5 + c4) * f.cos2) + 0.500 * (c3 - c4) * f.cos1;
    const double v135 = 0.2500 * ((c3 + 2.00 * c5 + c4) + (c3 - 2.00 * c5 + c4) * f.cos2) - 0.500 * (c3 - c4) * f.cos1;
    const double v90 = 0.500 * ((c1 + c2) + (c1 - c2) * f.cos2);
    set(p, i, four_orientations(v0, v45, v90, v135, t));
    if (i == 0) {
      const double dv45dphi = -0.500 * ((c3 - c4) * f.sin1 + (c3 - 2.00 * c5 + c4) * f.sin2);
      const double dv135dphi = 0.500 * ((c3 - c4) * f.sin1 - (c3 - 2.00 * c5 + c4) * f.sin2);
      const double dv90dphi = -(c1 - c2) * f.sin2;
      p.d_phi = four_orientations_phi_slope(dv45dphi, dv90dphi, dv135dphi, t);
      p.d_theta = four_orientations_theta_slope(v0, v45, v90, v135, t);
    }
  }
  return p;
}

// Fortran: C3v InterpolateOnT2H. Cuts: theta = 0; 90/30; 90/120; 45/30; 45/120; 45/210.
SiteValue c3v_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  constexpr double pis6 = pi / 6.00;
  const PhiTrig f = phi_trig(phi - pis6);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double c3 = vint(cuts, i, 3), c4 = vint(cuts, i, 4), c5 = vint(cuts, i, 5);
    const double v45 = 0.2500 * ((c3 + 2.00 * c4 + c5) + (c3 - 2.00 * c4 + c5) * f.cos2) + 0.500 * (c3 - c5) * f.cos1;
    const double v135 = 0.2500 * ((c3 + 2.00 * c4 + c5) + (c3 - 2.00 * c4 + c5) * f.cos2) - 0.500 * (c3 - c5) * f.cos1;
    const double v90 = 0.500 * ((c1 + c2) + (c1 - c2) * f.cos2);
    set(p, i, four_orientations(v0, v45, v90, v135, t));
    if (i == 0) {
      const double dv45dphi = -0.500 * ((c3 - c5) * f.sin1 + (c3 - 2.00 * c4 + c5) * f.sin2);
      const double dv135dphi = 0.500 * ((c3 - c5) * f.sin1 - (c3 - 2.00 * c4 + c5) * f.sin2);
      const double dv90dphi = -(c1 - c2) * f.sin2;
      p.d_phi = four_orientations_phi_slope(dv45dphi, dv90dphi, dv135dphi, t);
      p.d_theta = four_orientations_theta_slope(v0, v45, v90, v135, t);
    }
  }
  return p;
}

// Fortran: C3v InterpolateOnT2F. Cuts: theta = 0; 90/240; 90/330; 45/150; 45/240; 45/330.
SiteValue c3v_t2f(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  constexpr double pis6 = pi / 6.00;
  const PhiTrig f = phi_trig(phi - 11.00 * pis6);
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double c3 = vint(cuts, i, 3), c4 = vint(cuts, i, 4), c5 = vint(cuts, i, 5);
    const double v45 = 0.2500 * ((c3 + 2.00 * c4 + c5) + (c3 - 2.00 * c4 + c5) * f.cos2) + 0.500 * (-c3 + c5) * f.cos1;
    const double v135 = 0.2500 * ((c3 + 2.00 * c4 + c5) + (c3 - 2.00 * c4 + c5) * f.cos2) - 0.500 * (-c3 + c5) * f.cos1;
    const double v90 = 0.500 * ((c1 + c2) + (-c1 + c2) * f.cos2);
    set(p, i, four_orientations(v0, v45, v90, v135, t));
    if (i == 0) {
      const double dv45dphi = -0.500 * ((-c3 + c5) * f.sin1 + (c3 - 2.00 * c4 + c5) * f.sin2);
      const double dv135dphi = 0.500 * ((-c3 + c5) * f.sin1 - (c3 - 2.00 * c4 + c5) * f.sin2);
      const double dv90dphi = -(-c1 + c2) * f.sin2;
      p.d_phi = four_orientations_phi_slope(dv45dphi, dv90dphi, dv135dphi, t);
      p.d_theta = four_orientations_theta_slope(v0, v45, v90, v135, t);
    }
  }
  return p;
}

// Fortran: CRP_C4v_MoleculeSurface.f90 InterpolateOnTop (and InterpolateOnHollow, identical).
// Cuts: theta = 0; 90/0; 90/45.
SiteValue square_top(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double v90 = 0.5 * (c1 + c2) + 0.5 * (c1 - c2) * std::cos(4.0 * phi);
    set(p, i, 0.5 * (v0 + v90) + 0.5 * (v0 - v90) * t.cos2);
    if (i == 0) {
      const double dv0dphi = 0.0;
      const double dv90dphi = -4.0 * 0.5 * (c1 - c2) * std::sin(4.0 * phi);
      p.d_phi = 0.5 * (dv0dphi + dv90dphi) + 0.5 * (dv0dphi - dv90dphi) * t.cos2;
      p.d_theta = -2.0 * 0.5 * (v0 - v90) * t.sin2;
    }
  }
  return p;
}

SiteValue square_hollow(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  return square_top(cuts, t, phi);
}

// Fortran: CRP_C4v_MoleculeSurface.f90 InterpolateOnBridge. Cuts: theta = 0; 90/0; 90/45; 90/90.
SiteValue square_bridge(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2), c3 = vint(cuts, i, 3);
    const double v90 = (0.25 * c1 + 0.5 * c2 + 0.25 * c3) + 0.5 * (c1 - c3) * std::cos(2.0 * phi) +
                       (0.25 * c1 - 0.5 * c2 + 0.25 * c3) * std::cos(4.0 * phi);
    set(p, i, 0.5 * (v0 + v90) + 0.5 * (v0 - v90) * t.cos2);
    if (i == 0) {
      const double dv0dphi = 0.0;
      const double dv90dphi =
          -2.0 * 0.5 * (c1 - c3) * std::sin(2.0 * phi) - 4.0 * (0.25 * c1 - 0.5 * c2 + 0.25 * c3) * std::sin(4.0 * phi);
      p.d_phi = 0.5 * (dv0dphi + dv90dphi) + 0.5 * (dv0dphi - dv90dphi) * t.cos2;
      p.d_theta = -2.0 * 0.5 * (v0 - v90) * t.sin2;
    }
  }
  return p;
}

// Fortran: CRP_C4v_MoleculeSurface.f90 InterpolateOnT2H.
// Cuts: theta = 0; 45/225; 45/45; 90/45; 90/135; 45/135.
SiteValue square_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi) {
  const double phi45 = phi - pi / 4.0;          // Fortran: Phi - CRP_Constants_PI / 4.0
  const double phi135 = phi + pi - pi / 4.0;    // Fortran: Phi + CRP_Constants_PI - CRP_Constants_PI / 4.0
  SiteValue p{};
  for (int i = 0; i < 3; ++i) {
    const double v0 = vint(cuts, i, 0), c1 = vint(cuts, i, 1), c2 = vint(cuts, i, 2);
    const double c3 = vint(cuts, i, 3), c4 = vint(cuts, i, 4), c5 = vint(cuts, i, 5);
    // V(theta = 45, phi) = a0 + b0 cos(phi - 45) + c0 cos(2 (phi - 45)) through the three 45 cuts
    const double cc0 = 0.25 * (c1 + c2 - 2.0 * c5);
    const double b0 = -c5 + c2 - 2.0 * cc0;
    const double a0 = c2 - b0 - cc0;
    const double v45 = a0 + b0 * std::cos(phi45) + cc0 * std::cos(2.0 * phi45);
    const double v135 = a0 + b0 * std::cos(phi135) + cc0 * std::cos(2.0 * phi135);
    const double v90 = 0.5 * (c3 + c4) + 0.5 * (c3 - c4) * std::cos(2.0 * phi45);
    set(p, i,
        0.25 * (v0 + v45 + v90 + v135) + 0.5 * (v0 - v90) * t.cos2 + 0.5 * (v45 - v135) * t.sin2 +
            0.25 * (v0 + v90 - v45 - v135) * t.cos4);
    if (i == 0) {
      const double dv0dphi = 0.0;
      const double dv45dphi = -b0 * std::sin(phi45) - 2.0 * cc0 * std::sin(2.0 * phi45);
      const double dv135dphi = -b0 * std::sin(phi135) - 2.0 * cc0 * std::sin(2.0 * phi135);
      const double dv90dphi = -2.0 * 0.5 * (c3 - c4) * std::sin(2.0 * phi45);
      p.d_phi = 0.25 * (dv0dphi + dv45dphi + dv90dphi + dv135dphi) + 0.5 * (dv0dphi - dv90dphi) * t.cos2 +
                0.5 * (dv45dphi - dv135dphi) * t.sin2 + 0.25 * (dv0dphi + dv90dphi - dv45dphi - dv135dphi) * t.cos4;
      p.d_theta = -2.0 * 0.5 * (v0 - v90) * t.sin2 + 2.0 * 0.5 * (v45 - v135) * t.cos2 -
                  4.0 * 0.25 * (v0 + v90 - v45 - v135) * t.sin4;
    }
  }
  return p;
}

}  // namespace crplib::detail
