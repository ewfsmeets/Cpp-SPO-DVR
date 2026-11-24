#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <span>

#include "fft/fft.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

enum class TimeStep { Full, Half };

// Shape of the scattering wave function: (irot, iy + ny*ix, ir + nr*iz).
struct ScatteringShape {
  int n_rot;
  int nx;
  int ny;
  int nr;
  int nz;
};

// Number of adjacent (irot, iy) values gathered per (r, Z) work unit: up to 4 (one 64-byte
// cache line), and a divisor of n_roty = n_rot * ny so that all units have the same shape.
[[nodiscard]] int zr_chunk_size(int n_roty);
// Number of r points per (x, y) group: the largest divisor of nr whose blocks fit in max_bytes
// (1 MiB; at least 1), so all groups have the same shape and the work array stays cache-sized.
[[nodiscard]] int xy_group_points(const ScatteringShape& shape, std::size_t max_bytes = std::size_t(1) << 20);

// Whether T_xy can transform a group in the wave function itself (see
// ScatteringKineticOperator::transforms_in_place): one point per group, no momentum shift, and a block
// of n_rot * nx * ny complex values that is a multiple of 64 bytes.
[[nodiscard]] bool xy_in_place_possible(const ScatteringShape& shape, int group_points, bool shift);

// Whether T_xy shifts to the frame of the initial parallel momentum: some factor of shift_forward is not exactly 1
// (with EParallel = 0 every phase is +-0).
[[nodiscard]] bool xy_momentum_shift_needed(std::span<const std::complex<double>> shift_forward);

// How T_xy runs. The defaults are the program's; only tests change them, to reach both paths on small
// grids (a group byte limit other than xy_group_points' default gives groups other than those of the
// first touch and the potential step: a difference in scheduling only).
struct XyKineticOptions {
  std::size_t group_byte_limit = std::size_t(1) << 20;  // xy_group_points
  bool allow_in_place = true;                           // see ScatteringKineticOperator::transforms_in_place
};

// exp(-i T dt) on the scattering grid: T_Z + T_r by a 2D FFT in (r, Z), and T_xy (with the
// momentum shift of the parallel motion) by a 2D FFT in (x, y), done for a group of
// consecutive r points of one Z value at once (group size: see group_points()).
// Fortran: KineticOp.f90 PerformEXPKinOpZScat, PerformEXPKinOpRScat (there: separate 1D FFTs in
// Z and then r; the 2D form is the same operator, equal up to rounding); SplitOperator.f90
// PerformEXPKineticOp (MomentumShiftWF, TransToMomenta, MultByEXPKinEnergy, TransToCoords)
class ScatteringKineticOperator {
 public:
  ScatteringKineticOperator(const ScatteringShape& shape, const PropagationTables& tables, const Fft& fft,
                            const XyKineticOptions& options = {});

  // exp(-i (T_Z + T_r) dt). Fortran: KineticOpAlongZrOnScatGrid
  void apply_zr(ScatteringWave& psi, TimeStep step) const;
  // exp(-i T_xy dt) at every (Z, r) point. Fortran: PerformEXPKineticOp
  void apply_xy(ScatteringWave& psi, TimeStep step) const;
  // exp(-i T_xy dt) for group_points() consecutive (Z, r) points of one Z value; group = their
  // (irot, ixy) blocks, work: group_work_size() elements.
  void apply_xy_group(std::span<std::complex<double>> group, TimeStep step, std::span<std::complex<double>> work) const;
  // Number of r points per (x, y) group (xy_group_points).
  [[nodiscard]] int group_points() const { return group_points_; }
  // Elements of the work array apply_xy_group needs: none when T_xy runs in place.
  [[nodiscard]] std::size_t group_work_size() const;
  // Whether T_xy transforms a group in the wave function itself: with one point per group and no
  // momentum shift, the gathered work array would be a plain copy of the block (same layout), so the
  // same FFT plans run on the block directly (same results, bitwise). The block must also be 64-byte
  // aligned like the work array, so that the FFT library chooses the same (aligned) algorithm: n_rot ·
  // nx · ny · 16 bytes a multiple of 64, which every valid input meets (all FFT sizes are even).
  [[nodiscard]] bool transforms_in_place() const { return in_place_; }

 private:
  ScatteringShape shape_;
  const PropagationTables& tables_;
  int chunk_;   // number of adjacent (irot, iy) values per (r, Z) work unit
  bool shift_;  // momentum shift needed (non-zero parallel momentum)
  int group_points_;
  bool in_place_;
  // Coordinate -> momentum ("backward", sign -1, scale d/sqrt(2 pi) resp. sqrt(d/n) per
  // dimension) and back (sign +1, scale 1/(n * scale)). Fortran: BackwardFFT* / ForwardFFT*
  std::unique_ptr<ComplexFftPlan> zr_to_momentum_, zr_to_coordinate_;
  std::unique_ptr<ComplexFftPlan> xy_to_momentum_, xy_to_coordinate_;
};

}  // namespace spodvr
