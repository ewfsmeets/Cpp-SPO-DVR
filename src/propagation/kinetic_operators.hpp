#pragma once

#include <memory>

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

// exp(-i T dt) on the scattering grid: T_Z + T_r by a 2D FFT in (r, Z), and T_xy (with the
// momentum shift of the parallel motion) by a 2D FFT in (x, y), done for all r points of one Z
// value (a "group") at once.
// Fortran: KineticOp.f90 PerformEXPKinOpZScat, PerformEXPKinOpRScat (there: separate 1D FFTs in
// Z and then r; the 2D form is the same operator, equal up to rounding); SplitOperator.f90
// PerformEXPKineticOp (MomentumShiftWF, TransToMomenta, MultByEXPKinEnergy, TransToCoords)
class ScatteringKineticOperator {
 public:
  ScatteringKineticOperator(const ScatteringShape& shape, const PropagationTables& tables, const Fft& fft);

  // exp(-i (T_Z + T_r) dt). Fortran: KineticOpAlongZrOnScatGrid
  void apply_zr(ScatteringWave& psi, TimeStep step) const;
  // exp(-i T_xy dt) at every (Z, r) point. Fortran: PerformEXPKineticOp
  void apply_xy(ScatteringWave& psi, TimeStep step) const;
  // exp(-i T_xy dt) for the nr points of one Z value; group = their consecutive (irot, ixy)
  // blocks, work: group_work_size() elements.
  void apply_xy_group(std::span<std::complex<double>> group, TimeStep step, std::span<std::complex<double>> work) const;
  [[nodiscard]] std::size_t group_work_size() const;

 private:
  ScatteringShape shape_;
  const PropagationTables& tables_;
  int chunk_;   // number of adjacent (irot, iy) values per (r, Z) work unit
  bool shift_;  // momentum shift needed (non-zero parallel momentum)
  // Coordinate -> momentum ("backward", sign -1, scale d/sqrt(2 pi) resp. sqrt(d/n) per
  // dimension) and back (sign +1, scale 1/(n * scale)). Fortran: BackwardFFT* / ForwardFFT*
  std::unique_ptr<ComplexFftPlan> zr_to_momentum_, zr_to_coordinate_;
  std::unique_ptr<ComplexFftPlan> xy_to_momentum_, xy_to_coordinate_;
};

}  // namespace spodvr
