#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <span>
#include <vector>

#include "fft/fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"

namespace spodvr {

// Transformation of rotational wave-function columns between the FBR basis |J mJ>
// (RotationalBasis order) and the (theta, phi) DVR grid. A DVR block has the Fortran layout
// (iphi, itheta) with leading dimension nphi + 1 (one padding element per theta). n_columns
// columns are transformed together, so the phi FFTs of all of them are one batched call (the CPU
// potential step: columns_per_pass of the (x, y) points of one (Z, r) point at a time; the frames:
// all of them): column c is at fbr[n_rot * c ...] and dvr[dvr_block_size() * c ...].
//   FBR -> DVR: psi(theta_t, phi_k) = sum_{J,m} sqrt(w_t / nphi) P_J^m(theta_t) exp(i m phi_k) c_{Jm}
//   DVR -> FBR: the adjoint, so DVR -> FBR after FBR -> DVR is the identity.
// Both directions may be called concurrently on different arrays.
// Fortran: transform.f90 (FBRToDVRGenComplexTransform, DVRToFBRGenComplexTransform)
class AngularTransform {
 public:
  AngularTransform(const AngularGrid& grid, const RotationalBasis& basis, const Fft& fft, int n_columns = 1);

  [[nodiscard]] int n_theta() const { return n_theta_; }
  [[nodiscard]] int n_phi() const { return n_phi_; }
  [[nodiscard]] int leading_dimension() const { return n_phi_ + 1; }
  [[nodiscard]] int dvr_block_size() const { return dvr_block_size(n_theta_, n_phi_); }
  // Size of one DVR block (with padding) for a grid of n_theta x n_phi points.
  [[nodiscard]] static int dvr_block_size(int n_theta, int n_phi) { return (n_phi + 1) * n_theta; }
  [[nodiscard]] int fbr_size() const { return n_rot_; }
  [[nodiscard]] int n_columns() const { return n_columns_; }

  // The columns the CPU transforms per call: the largest divisor of n_xy whose DVR blocks fit
  // cpu_pass_bytes, at least 1. A constant, so results never depend on the machine; the sums of a
  // column keep their order whatever the count. Small passes keep one thread's DVR work in its L2
  // instead of megabytes per thread (that this is why they are faster is inferred, not measured).
  // The limit is four DZ blocks (4 x 14 640 B), rounded up to a KiB: on the DZ grid (144 columns,
  // 8 threads, Ryzen 7 5800X, 2026-11-22) the potential step took 1.45, 1.36, 1.34, 1.31, 1.32, 1.29
  // and 1.62 s with 1, 2, 3, 4, 6, 12 and 144 columns per pass (mean of two runs; 3 to 12 are equal
  // within the run-to-run noise of about 6 %), and 4 is the smallest count within 2 % of the best.
  [[nodiscard]] static int columns_per_pass(int n_xy, int n_theta, int n_phi);
  static constexpr std::size_t cpu_pass_bytes = std::size_t(58) << 10;

  // fbr: n_rot * n_columns values; dvr: dvr_block_size() * n_columns values (overwritten,
  // padding left untouched).
  void fbr_to_dvr(std::span<const std::complex<double>> fbr, std::span<std::complex<double>> dvr) const;
  // dvr is used as work space (overwritten); fbr: n_rot * n_columns values (overwritten).
  void dvr_to_fbr(std::span<std::complex<double>> dvr, std::span<std::complex<double>> fbr) const;

  // The FBR functions with one mJ value and the phi slot (FFT index) of that mJ.
  struct MjBlock {
    IndexRange fbr;
    int phi_slot;
  };
  // The tables of the transform (for other implementations of it, e.g. on a GPU): the mJ blocks
  // in the order mJ = -MjMax .. MjMax, and sqrt(w_t / nphi) P_J^m(theta_t) ordered (mJ block,
  // theta, J).
  [[nodiscard]] const std::vector<MjBlock>& mj_blocks() const { return blocks_; }
  [[nodiscard]] const std::vector<double>& matrix() const { return matrix_; }
  [[nodiscard]] int mj_max() const { return mj_max_; }

 private:
  int n_theta_;
  int n_phi_;
  int n_rot_;
  int mj_max_;
  int n_columns_;
  std::vector<MjBlock> blocks_;
  // sqrt(w_t / nphi) P_J^m(theta_t), ordered (mJ block, theta, J). Fortran: FBRToDVRMatrices
  std::vector<double> matrix_;
  std::unique_ptr<ComplexFftPlan> to_dvr_fft_;  // Fortran: FFTForFBRToDVRInPhi (sign +1)
  std::unique_ptr<ComplexFftPlan> to_fbr_fft_;  // Fortran: FFTForDVRToFBRInPhi (sign -1)
};

}  // namespace spodvr
