#pragma once

#include <optional>
#include <utility>

#include "potential/potential.hpp"
#include "symmetry/surface_operation.hpp"

namespace spodvr {

// An integer 2x2 matrix acting on lattice (fractional) coordinates (u, w): (u', w') = A (u, w).
struct IntegerMatrix2 {
  int a11, a12, a21, a22;
  friend bool operator==(const IntegerMatrix2&, const IntegerMatrix2&) = default;
};

// The 2D surface lattice: lattice vectors a1 = alat (1, 0) and a2 = blat (cos skew, sin skew) in
// Cartesian coordinates. Grid coordinates (x, y) are skewed: the Cartesian position of (x, y) is
// x (1, 0) + y (cos skew, sin skew), with x in [0, alat) and y in [0, blat).
// Fortran: alat, blat, skewingAngle (H2SurfaceSystem.f90 GET_PARAMETERS)
class SurfaceLattice {
 public:
  SurfaceLattice(double lattice_constant_x, double lattice_constant_y, double skewing_angle);
  explicit SurfaceLattice(const SystemParameters& system);

  [[nodiscard]] double lattice_constant_x() const { return a_; }
  [[nodiscard]] double lattice_constant_y() const { return b_; }
  [[nodiscard]] double skewing_angle() const { return skew_; }

  // Skewed grid coordinates (x, y) (bohr along a1 and a2) <-> Cartesian coordinates (bohr).
  [[nodiscard]] CartesianVector to_cartesian(double x, double y) const;
  [[nodiscard]] std::pair<double, double> to_skewed(const CartesianVector& v) const;

  // The matrix of operation g in lattice coordinates, A = B^-1 M_g B with B = [a1 a2] (the
  // fractional coordinates are u = x/alat, w = y/blat). Empty if g is not a symmetry of the
  // lattice (A not integer within 1e-10), e.g. a 90 deg rotation of a rectangular cell.
  [[nodiscard]] std::optional<IntegerMatrix2> lattice_matrix(const SurfaceOperation& g) const;

  // Squared length of the reciprocal lattice vector n b1 + m b2 in units of (2 pi)^2 / bohr^2.
  // Same metric as the parallel kinetic energy (kx^2 + ky^2 - 2 cos(skew) kx ky).
  [[nodiscard]] double reciprocal_length_squared(int n, int m) const;

  // Diffraction order of channel (n, m): the number of steps from (0, 0) to (n, m) on the
  // reciprocal lattice, each step one of the shortest reciprocal vectors (+-b1, +-b2, and
  // +-(b1 + b2) or +-(b1 - b2) when that is not longer than the longer of b1 and b2).
  // Rectangular and square cells: |n| + |m|; hexagonal cells (skew 60 deg, a = b, b1 + b2 short):
  // max(|n|, |m|) if n m > 0, else |n| + |m|. This reproduces the Fortran table, which is per
  // point group (C2v, C4v; C3, C3v, C6v), and also works for any other cell.
  // Fortran: RotDifBasis.f90 getDiffractionOrder
  [[nodiscard]] int diffraction_order(int n, int m) const;

  // Which diagonal reciprocal vectors are among the shortest steps.
  [[nodiscard]] bool has_step_sum() const { return step_sum_; }            // b1 + b2
  [[nodiscard]] bool has_step_difference() const { return step_difference_; }  // b1 - b2

 private:
  double a_;
  double b_;
  double skew_;
  bool step_sum_ = false;
  bool step_difference_ = false;
};

}  // namespace spodvr
