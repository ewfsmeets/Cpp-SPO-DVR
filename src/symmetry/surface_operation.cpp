// Rotations and reflections of a surface point group (see surface_operation.hpp).
#include "symmetry/surface_operation.hpp"

#include <cmath>

#include "core/constants.hpp"

namespace spodvr {

namespace {
int reduce_twelfths(int k) { return ((k % 12) + 12) % 12; }
double twelfths_to_rad(int k) { return Constants::two_pi * k / 12.0; }
}  // namespace

SurfaceOperation compose(const SurfaceOperation& a, const SurfaceOperation& b) {
  const int sign_a = a.mirror ? -1 : 1;
  return {a.mirror != b.mirror, reduce_twelfths(a.angle_twelfths + sign_a * b.angle_twelfths)};
}

SurfaceOperation inverse(const SurfaceOperation& g) {
  if (g.mirror) return g;
  return {false, reduce_twelfths(-g.angle_twelfths)};
}

double apply_to_azimuth(const SurfaceOperation& g, double phi) {
  const double beta = twelfths_to_rad(g.angle_twelfths);
  return g.mirror ? beta - phi : beta + phi;
}

CartesianVector apply_to_cartesian(const SurfaceOperation& g, const CartesianVector& v) {
  const double beta = twelfths_to_rad(g.angle_twelfths);
  const double c = std::cos(beta);
  const double s = std::sin(beta);
  const double y = g.mirror ? -v.y : v.y;
  return {c * v.x - s * y, s * v.x + c * y};
}

std::string describe(const SurfaceOperation& g) {
  if (!g.mirror) {
    if (g.angle_twelfths == 0) return "identity";
    return "rotation " + std::to_string(30 * g.angle_twelfths) + " deg";
  }
  // The mirror line is at beta/2 = 15 * angle_twelfths degrees (0 <= angle < 180).
  return "mirror at " + std::to_string(15 * g.angle_twelfths) + " deg";
}

}  // namespace spodvr
