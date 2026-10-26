#pragma once

#include <span>
#include <vector>

namespace spodvr {

// The (theta, phi) DVR grid: Gauss-Legendre points in theta (half a grid, 0 < theta <= pi/2,
// for a homonuclear molecule) and equidistant points in phi with an FFT-compatible size.
// Fortran: setupdvr.f90 SetupDVRProj + angulargrid.f90 InitAngularGrid (SymmetryInPhi = false)
class AngularGrid {
 public:
  // j_max: maximum J of the rotational basis; mj_max_basis: maximum |mJ| (already with the
  // 'MJMax < 1 -> JMax' rule applied); num_phis_multiple: nphi must be a multiple of this
  // (1 unless surface symmetry is exploited). Throws InputError for invalid combinations.
  AngularGrid(int j_max, int mj_max_basis, bool homonuclear, int j_in, int mj_in,
              int num_phis_multiple = 1);

  [[nodiscard]] int n_theta() const { return static_cast<int>(thetas_.size()); }
  [[nodiscard]] int n_phi() const { return static_cast<int>(phis_.size()); }
  [[nodiscard]] double d_phi() const { return d_phi_; }
  [[nodiscard]] int j_max() const { return j_max_; }  // after parity adjustment
  [[nodiscard]] int j_min() const { return j_min_; }
  [[nodiscard]] int delta_j() const { return delta_j_; }
  [[nodiscard]] int mj_max() const { return mj_max_; }
  [[nodiscard]] bool homonuclear() const { return homonuclear_; }
  [[nodiscard]] std::span<const double> thetas() const { return thetas_; }
  [[nodiscard]] std::span<const double> theta_weights() const { return weights_; }
  [[nodiscard]] std::span<const double> phis() const { return phis_; }

 private:
  int j_max_;
  int j_min_;
  int delta_j_;
  int mj_max_;
  bool homonuclear_;
  double d_phi_;
  std::vector<double> thetas_;
  std::vector<double> weights_;
  std::vector<double> phis_;
};

}  // namespace spodvr
