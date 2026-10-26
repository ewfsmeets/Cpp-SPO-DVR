#pragma once

#include <complex>
#include <memory>
#include <span>
#include <vector>

#include "fft/fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"

namespace spodvr {

// Transformation of one rotational wave-function column between the FBR basis |J mJ>
// (RotationalBasis order) and the (theta, phi) DVR grid. The DVR block has the Fortran layout
// (iphi, itheta) with leading dimension nphi + 1 (one padding element per theta).
//   FBR -> DVR: psi(theta_t, phi_k) = sum_{J,m} sqrt(w_t / nphi) P_J^m(theta_t) exp(i m phi_k) c_{Jm}
//   DVR -> FBR: the adjoint, so DVR -> FBR after FBR -> DVR is the identity.
// Both directions may be called concurrently on different arrays.
// Fortran: transform.f90 (FBRToDVRGenComplexTransform, DVRToFBRGenComplexTransform)
class AngularTransform {
 public:
  AngularTransform(const AngularGrid& grid, const RotationalBasis& basis, const Fft& fft);

  [[nodiscard]] int n_theta() const { return n_theta_; }
  [[nodiscard]] int n_phi() const { return n_phi_; }
  [[nodiscard]] int leading_dimension() const { return n_phi_ + 1; }
  [[nodiscard]] int dvr_block_size() const { return (n_phi_ + 1) * n_theta_; }
  [[nodiscard]] int fbr_size() const { return n_rot_; }

  // fbr: n_rot values; dvr: dvr_block_size() values (overwritten, padding left untouched).
  void fbr_to_dvr(std::span<const std::complex<double>> fbr, std::span<std::complex<double>> dvr) const;
  // dvr is used as work space (overwritten); fbr: n_rot values (overwritten).
  void dvr_to_fbr(std::span<std::complex<double>> dvr, std::span<std::complex<double>> fbr) const;

 private:
  // The FBR functions with one mJ value and the phi slot (FFT index) of that mJ.
  struct MjBlock {
    IndexRange fbr;
    int phi_slot;
  };

  int n_theta_;
  int n_phi_;
  int n_rot_;
  int mj_max_;
  std::vector<MjBlock> blocks_;
  // sqrt(w_t / nphi) P_J^m(theta_t), ordered (mJ block, theta, J). Fortran: FBRToDVRMatrices
  std::vector<double> matrix_;
  std::unique_ptr<ComplexFftPlan> to_dvr_fft_;  // Fortran: FFTForFBRToDVRInPhi (sign +1)
  std::unique_ptr<ComplexFftPlan> to_fbr_fft_;  // Fortran: FFTForDVRToFBRInPhi (sign -1)
};

}  // namespace spodvr
