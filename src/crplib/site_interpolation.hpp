#pragma once
// Interpolation over the molecular orientation (theta, phi) at one surface site: from the values of
// the site's cuts (each a fixed orientation, see docs/crp_library.md Appendix A) at the current
// (r, Z) to the value at (theta, phi), with its derivatives. Fortran: the InterpolateOn* routines
// of CRP_C3v/C4v/C6v_MoleculeSurface.f90. Formulas in the Fortran order of operations.

#include <span>

#include "crplib/spline.hpp"

namespace crplib::detail {

// A site's interpolated value and derivatives. Fortran: Psite(1:5) = value, d/dZ, d/dr, d/dphi,
// d/dtheta (the cut derivatives come from the bicubic spline: d_dy = d/dZ, d_dx = d/dr).
struct SiteValue {
  double v, d_z, d_r, d_phi, d_theta;
};

// The cos/sin of 2 theta and 4 theta, shared by every site of one evaluation (speed-up S5; the
// same expressions as the Fortran, so the values are identical). sin 4 theta enters only the
// derivatives (0 for Want::value, S3).
struct ThetaTrig {
  double cos2, sin2, cos4, sin4;
};
ThetaTrig theta_trig(double theta, Want want = Want::value_and_gradient);

// Every site function: cuts are the site's cut values in the cut order of Appendix A; with
// Want::value only SiteValue::v is computed (the rest is 0; speed-up S3).

// Hexagonal (C3v, C6v). Fortran: InterpolateOnTop (3 cuts), InterpolateOnBridge (4),
// InterpolateOnFCC (C6v) = InterpolateOnHCP (C3v) (5), C6v InterpolateOnT2F (6), C3v
// InterpolateOnT2H (6), C3v InterpolateOnT2F (6)
SiteValue hexagonal_top(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                        Want want = Want::value_and_gradient);
SiteValue hexagonal_bridge(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                           Want want = Want::value_and_gradient);
SiteValue hexagonal_hollow(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                           Want want = Want::value_and_gradient);
SiteValue c6v_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                  Want want = Want::value_and_gradient);
SiteValue c3v_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                  Want want = Want::value_and_gradient);
SiteValue c3v_t2f(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                  Want want = Want::value_and_gradient);

// Square (C4v). Fortran: CRP_C4v_MoleculeSurface.f90 InterpolateOnTop (3 cuts), InterpolateOnBridge
// (4), InterpolateOnHollow (3), InterpolateOnT2H (6)
SiteValue square_top(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                     Want want = Want::value_and_gradient);
SiteValue square_bridge(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                        Want want = Want::value_and_gradient);
SiteValue square_hollow(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                        Want want = Want::value_and_gradient);
SiteValue square_t2h(std::span<const BicubicValue> cuts, const ThetaTrig& t, double phi,
                     Want want = Want::value_and_gradient);

}  // namespace crplib::detail
