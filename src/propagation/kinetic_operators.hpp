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

// exp(-i T dt) on the scattering grid: T_Z and T_r by FFT over the whole wave function, and
// T_xy (with the momentum shift of the parallel motion) per (Z, r) point.
// Fortran: KineticOp.f90 PerformEXPKinOpZScat, PerformEXPKinOpRScat; SplitOperator.f90
// PerformEXPKineticOp (MomentumShiftWF, TransToMomenta, MultByEXPKinEnergy, TransToCoords)
class ScatteringKineticOperator {
 public:
  ScatteringKineticOperator(const ScatteringShape& shape, const PropagationTables& tables, const Fft& fft);

  // exp(-i T_Z dt) then exp(-i T_r dt). Fortran: KineticOpAlongZrOnScatGrid
  void apply_zr(ScatteringWave& psi, TimeStep step) const;
  // exp(-i T_xy dt) at every (Z, r) point. Fortran: PerformEXPKineticOp
  void apply_xy(ScatteringWave& psi, TimeStep step) const;

 private:
  void apply_xy_block(std::span<std::complex<double>> block, const std::vector<std::complex<double>>& factors) const;

  ScatteringShape shape_;
  const PropagationTables& tables_;
  // Coordinate -> momentum ("backward", sign -1, scale d/sqrt(2 pi) or sqrt(d/n)) and back
  // (sign +1, scale 1/(n * scale)). Fortran: BackwardFFT* / ForwardFFT*
  std::unique_ptr<ComplexFftPlan> z_to_momentum_, z_to_coordinate_;
  std::unique_ptr<ComplexFftPlan> r_to_momentum_, r_to_coordinate_;
  std::unique_ptr<ComplexFftPlan> x_to_momentum_, x_to_coordinate_;
  std::unique_ptr<ComplexFftPlan> y_to_momentum_, y_to_coordinate_;
};

}  // namespace spodvr
