#pragma once
// A 3D cubic spline, periodic in x and y and clamped (given first derivatives) in z: the C4v
// atom-surface interpolation. Fortran: Spline.f90 Spline3DPeriodic, SPLINE3D, int_PES6CuSPL2DPER,
// SPLINT3D (used through Compute_/Evaluate_Spline3D_Periodic).

#include <vector>

#include "crplib/spline.hpp"

namespace crplib::detail {

class PeriodicSpline3D {
 public:
  PeriodicSpline3D() = default;
  // f(i, j, k) = f[i + nx * (j + ny * k)] at (x_i, y_j, z_k); the periods are x_{nx-1} - x_0 and
  // y_{ny-1} - y_0. dfdz_first(i, j) = dfdz_first[i + nx * j] is dF/dz at z_0, dfdz_last at
  // z_{nz-1}. Fortran: Compute_Spline3D_Periodic (SPLINE3D with IC1 = ICN = 2)
  PeriodicSpline3D(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<double> f,
                   const std::vector<double>& dfdz_first, const std::vector<double>& dfdz_last);

  struct Value {
    double value;
    double d_dx;
    double d_dy;
    double d_dz;
  };
  // Value and gradient at (x, y, z) inside the grid (Want::value: the gradient is 0); throws
  // EvaluationError outside. Fortran: Evaluate_Spline3D_Periodic (SPLINT3D, which STOPs outside)
  [[nodiscard]] Value evaluate(double x, double y, double z, Want want = Want::value_and_gradient) const;

  [[nodiscard]] const std::vector<double>& z() const { return z_; }

  // The coefficient arrays (layout of f), for comparisons with the Fortran. Fortran: DFDX, DFDY,
  // DFDXDY, DFDZ, DFDZDX, DFDZDY, DFDZDXDY (second derivatives, despite the names)
  struct Coefficients {
    std::vector<double> xx, yy, xxyy, zz, zzxx, zzyy, zzxxyy;
  };
  [[nodiscard]] const Coefficients& coefficients() const { return c_; }

 private:
  [[nodiscard]] int index(int i, int j, int k) const {
    return i + static_cast<int>(x_.size()) * (j + static_cast<int>(y_.size()) * k);
  }

  std::vector<double> x_, y_, z_, f_;
  Coefficients c_;
};

}  // namespace crplib::detail
