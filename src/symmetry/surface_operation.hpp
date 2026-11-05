#pragma once

#include <string>

namespace spodvr {

// A symmetry operation of a surface point group: a rotation about the surface normal or a
// reflection in a plane containing the normal. With beta = 2 pi angle_twelfths / 12 it maps
//   the azimuth          phi -> beta + phi   (rotation)   or   beta - phi   (mirror),
//   Cartesian (x, y)     by R(beta)          (rotation)   or   R(beta) diag(1, -1) (mirror),
// i.e. a mirror reflects in the line at angle beta/2. Z, r and theta are unchanged. Multiples of
// 30 deg cover all 2D crystallographic point groups (C1..C6, Cs, C2v, C3v, C4v, C6v).
struct SurfaceOperation {
  bool mirror = false;
  int angle_twelfths = 0;  // 0 .. 11

  friend bool operator==(const SurfaceOperation&, const SurfaceOperation&) = default;
};

// The identity operation.
[[nodiscard]] constexpr SurfaceOperation identity_operation() { return {false, 0}; }

// a o b: first b, then a. phi -> s_a (s_b phi + beta_b) + beta_a.
[[nodiscard]] SurfaceOperation compose(const SurfaceOperation& a, const SurfaceOperation& b);

// The inverse: a rotation by -beta; a mirror is its own inverse.
[[nodiscard]] SurfaceOperation inverse(const SurfaceOperation& g);

// The image of the azimuth phi (rad), not reduced to [0, 2 pi).
[[nodiscard]] double apply_to_azimuth(const SurfaceOperation& g, double phi);

// The image of the Cartesian vector (x, y).
struct CartesianVector {
  double x;
  double y;
};
[[nodiscard]] CartesianVector apply_to_cartesian(const SurfaceOperation& g, const CartesianVector& v);

// Short description for reports, e.g. "rotation 120 deg", "mirror at 90 deg".
[[nodiscard]] std::string describe(const SurfaceOperation& g);

}  // namespace spodvr
