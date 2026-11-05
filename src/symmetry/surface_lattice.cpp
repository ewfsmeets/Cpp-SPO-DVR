// The surface lattice: reciprocal metric and the diffraction order of a channel.
// Fortran: RotDifBasis.f90 getDiffractionOrder (generalised from per point group to per lattice)
#include "symmetry/surface_lattice.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace spodvr {

namespace {
// Relative tolerance for "not longer than": the hexagonal case |b1 + b2| = |b1| is an exact tie.
constexpr double length_tolerance = 1e-12;
}  // namespace

SurfaceLattice::SurfaceLattice(double lattice_constant_x, double lattice_constant_y, double skewing_angle)
    : a_(lattice_constant_x), b_(lattice_constant_y), skew_(skewing_angle) {
  const double longest_axis = std::max(reciprocal_length_squared(1, 0), reciprocal_length_squared(0, 1));
  const double limit = longest_axis * (1.0 + length_tolerance);
  step_sum_ = reciprocal_length_squared(1, 1) <= limit;
  step_difference_ = reciprocal_length_squared(1, -1) <= limit;
}

SurfaceLattice::SurfaceLattice(const SystemParameters& system)
    : SurfaceLattice(system.lattice_constant_x, system.lattice_constant_y, system.skewing_angle) {}

CartesianVector SurfaceLattice::to_cartesian(double x, double y) const {
  return {x + y * std::cos(skew_), y * std::sin(skew_)};
}

std::pair<double, double> SurfaceLattice::to_skewed(const CartesianVector& v) const {
  const double y = v.y / std::sin(skew_);
  return {v.x - y * std::cos(skew_), y};
}

// Images of the lattice vectors a1, a2 under g, expressed in the basis (a1, a2): the columns of A.
std::optional<IntegerMatrix2> SurfaceLattice::lattice_matrix(const SurfaceOperation& g) const {
  constexpr double integer_tolerance = 1e-10;
  const CartesianVector image_a1 = apply_to_cartesian(g, to_cartesian(a_, 0.0));
  const CartesianVector image_a2 = apply_to_cartesian(g, to_cartesian(0.0, b_));
  const auto [x1, y1] = to_skewed(image_a1);
  const auto [x2, y2] = to_skewed(image_a2);
  const double entries[4] = {x1 / a_, x2 / a_, y1 / b_, y2 / b_};  // a11, a12, a21, a22
  int rounded[4] = {0, 0, 0, 0};
  for (int i = 0; i < 4; ++i) {
    rounded[i] = static_cast<int>(std::lround(entries[i]));
    if (std::abs(entries[i] - rounded[i]) > integer_tolerance) return std::nullopt;
  }
  return IntegerMatrix2{rounded[0], rounded[1], rounded[2], rounded[3]};
}

// |n b1 + m b2|^2 with a_i . b_j = delta_ij: (n^2/a^2 + m^2/b^2 - 2 n m cos/(a b)) / sin^2.
double SurfaceLattice::reciprocal_length_squared(int n, int m) const {
  const double sin_skew = std::sin(skew_);
  const double nx = n / a_;
  const double my = m / b_;
  return (nx * nx + my * my - 2.0 * std::cos(skew_) * nx * my) / (sin_skew * sin_skew);
}

// Steps along a diagonal that is short cover |n| and |m| at the same time, so where n and m have
// the sign of that diagonal the distance is max(|n|, |m|); elsewhere only axis steps help.
// (At most one diagonal can be short in a cell with nonzero area.)
int SurfaceLattice::diffraction_order(int n, int m) const {
  const int abs_n = std::abs(n);
  const int abs_m = std::abs(m);
  const bool along_sum = step_sum_ && n * m > 0;
  const bool along_difference = step_difference_ && n * m < 0;
  if (along_sum || along_difference) return std::max(abs_n, abs_m);
  return abs_n + abs_m;
}

}  // namespace spodvr
