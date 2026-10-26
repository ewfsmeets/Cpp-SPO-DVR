// (theta, phi) DVR grid. Fortran: setupdvr.f90 SetupDVRProj, angulargrid.f90 InitAngularGrid
#include "grids/angular_grid.hpp"

#include <cmath>
#include <cstdlib>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "fft/fft_sizes.hpp"
#include "grids/gauss_legendre.hpp"

namespace spodvr {

AngularGrid::AngularGrid(int j_max, int mj_max_basis, bool homonuclear, int j_in, int mj_in,
                         int num_phis_multiple)
    : j_max_(j_max), j_min_(0), delta_j_(1), mj_max_(mj_max_basis), homonuclear_(homonuclear) {
  // Fortran: setupdvr.f90 SetupDVRProj input checks (also covered by ValidateInput).
  if (j_in < 0) throw InputError("JIn", "JIn cannot be negative.");
  if (j_max < 0) throw InputError("JMax", "JMax cannot be negative.");
  if (std::abs(mj_in) > j_in) throw InputError("MJIn", "The absolute value of MjIn is bigger than JIn.");
  if (homonuclear && (j_in % 2) != (j_max % 2))
    throw InputError("JMax", "JMax and JIn must have the same parity for a homonuclear molecule.");
  if (mj_max_basis > j_max)
    throw InputError("MJMax", "The MjMax parameter was set higher than the JMax parameter.");
  if (j_max < j_in) throw InputError("JIn", "JIn > JMax.");  // Fortran: silent RETURN in InitAngularGrid

  // Number of phi points. Fortran: setupdvr.f90 Setup_FFTSize, Setup_dPhi
  const int n_phi = nearest_fft_size(2 * mj_max_basis + 1, num_phis_multiple);
  d_phi_ = Constants::two_pi / n_phi;

  // Fortran: transform.f90 InitGenComplexTransform silently returns in these cases; here they
  // are hard errors.
  if (n_phi != nearest_fft_size(n_phi, num_phis_multiple))
    throw InputError("MJMax", "the number of phi points is not a valid FFT size");
  if (n_phi < 2 * mj_max_basis + 1)
    throw InputError("MJMax", "the number of phi points must be at least 2*MJMax+1");

  // Fortran: angulargrid.f90 InitAngularGrid (SymInPhi = .FALSE. => JMin = 0)
  int n_theta = 0;
  if (homonuclear) {
    if (j_in % 2 == 0) {
      if (j_max_ % 2 != 0) j_max_ -= 1;
      if (j_min_ % 2 != 0) j_min_ += 1;
    } else {
      if (j_max_ % 2 == 0) j_max_ -= 1;
      if (j_min_ % 2 == 0) j_min_ += 1;
    }
    n_theta = (j_max_ - j_min_) / 2 + 1;
    delta_j_ = 2;
  } else {
    n_theta = j_max_ - j_min_ + 1;
    delta_j_ = 1;
  }

  phis_.resize(n_phi);
  for (int i = 0; i < n_phi; ++i) phis_[i] = i * d_phi_;

  const QuadratureRule rule = gauss_associated_legendre(1.0, -1.0, j_max_ + 1, 0);
  thetas_.resize(n_theta);
  weights_.assign(rule.weights.begin(), rule.weights.begin() + n_theta);

  double norm = -1.0;
  if (homonuclear) {
    norm = -2.0;
    // With an odd number of roots the middle one (theta = pi/2) is not doubled.
    if (j_max_ % 2 == 0) weights_.back() = 0.5 * weights_.back();
  }
  for (std::size_t t = 0; t < weights_.size(); ++t) {
    weights_[t] = weights_[t] * norm;
    thetas_[t] = std::acos(rule.roots[t]);
  }
}

}  // namespace spodvr
