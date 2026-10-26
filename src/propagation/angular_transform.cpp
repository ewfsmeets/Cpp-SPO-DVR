// FBR <-> DVR transformation in (theta, phi): Legendre sums per mJ block plus an FFT over phi.
// Fortran: transform.f90
#include "propagation/angular_transform.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "grids/associated_legendre.hpp"

namespace spodvr {

// Fortran: transform.f90 InitGenComplexTransform
AngularTransform::AngularTransform(const AngularGrid& grid, const RotationalBasis& basis, const Fft& fft)
    : n_theta_(grid.n_theta()), n_phi_(grid.n_phi()), n_rot_(basis.size()), mj_max_(basis.mj_max()) {
  if (n_phi_ < 2 * mj_max_ + 1) throw std::invalid_argument("AngularTransform: nphi < 2*MjMax+1");

  for (int m = -mj_max_; m <= mj_max_; ++m) {
    const int phi_slot = m >= 0 ? m : n_phi_ + m;
    blocks_.push_back({basis.mj_block(m), phi_slot});
  }

  for (int m = -mj_max_; m <= mj_max_; ++m) {
    for (int t = 0; t < n_theta_; ++t) {
      const std::vector<double> p = associated_legendre(grid.j_max(), m, grid.thetas()[t]);
      const double factor = std::sqrt(grid.theta_weights()[t] / n_phi_);
      for (int j = grid.j_min(); j <= grid.j_max(); j += grid.delta_j()) {
        if (std::abs(m) > j) continue;
        matrix_.push_back(p[j - std::abs(m)] * factor);
      }
    }
  }

  // One FFT of length nphi per theta, the theta columns nphi + 1 elements apart.
  const FftLayout layout{n_phi_, n_theta_, 1, n_phi_ + 1};
  to_dvr_fft_ = fft.plan_complex(layout, FftSign::Plus, 1.0);
  to_fbr_fft_ = fft.plan_complex(layout, FftSign::Minus, 1.0);
}

// Fortran: transform.f90 FBRToDVRGenComplexTransform
void AngularTransform::fbr_to_dvr(std::span<const std::complex<double>> fbr,
                                  std::span<std::complex<double>> dvr) const {
  const int ld = leading_dimension();
  int index = 0;
  for (const MjBlock& block : blocks_) {
    for (int t = 0; t < n_theta_; ++t) {
      std::complex<double> sum(0.0, 0.0);
      for (int irot = block.fbr.first; irot <= block.fbr.last; ++irot) {
        sum = sum + matrix_[index] * fbr[irot];
        ++index;
      }
      dvr[block.phi_slot + ld * t] = sum;
    }
  }

  // The phi slots of |mJ| > MjMax carry no basis function.
  for (int t = 0; t < n_theta_; ++t)
    for (int k = mj_max_ + 1; k <= n_phi_ - mj_max_ - 1; ++k) dvr[k + ld * t] = 0.0;

  to_dvr_fft_->execute(dvr);
}

// Fortran: transform.f90 DVRToFBRGenComplexTransform
void AngularTransform::dvr_to_fbr(std::span<std::complex<double>> dvr,
                                  std::span<std::complex<double>> fbr) const {
  to_fbr_fft_->execute(dvr);

  const int ld = leading_dimension();
  std::fill_n(fbr.begin(), n_rot_, std::complex<double>(0.0, 0.0));
  int index = 0;
  for (const MjBlock& block : blocks_) {
    for (int t = 0; t < n_theta_; ++t) {
      const std::complex<double> value = dvr[block.phi_slot + ld * t];
      for (int irot = block.fbr.first; irot <= block.fbr.last; ++irot) {
        fbr[irot] = fbr[irot] + matrix_[index] * value;
        ++index;
      }
    }
  }
}

}  // namespace spodvr
