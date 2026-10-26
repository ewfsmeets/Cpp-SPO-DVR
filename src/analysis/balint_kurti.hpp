#pragma once

#include <complex>
#include <memory>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Time-dependent Balint-Kurti coefficients c(is, v; t): the overlap of the wave packet at the
// analysis surface Z = Zana with the asymptotic channel states (rovibrational eigenfunction
// times diffraction plane wave). Layout (is, v): is + nstata * v. Fortran: balcf(nstata, nviba)
struct BkCoefficients {
  int n_states = 0;
  int n_v = 0;
  std::vector<std::complex<double>> values;

  BkCoefficients(int states, int vibrations)
      : n_states(states), n_v(vibrations), values(states * vibrations, std::complex<double>(0.0, 0.0)) {}
  [[nodiscard]] std::complex<double>& operator()(int is, int v) { return values[is + n_states * v]; }
  [[nodiscard]] const std::complex<double>& operator()(int is, int v) const { return values[is + n_states * v]; }
};

// Computes the BK coefficients of the current wave packet, once per analysis step.
// Scattering grid, per analysed (v, J, mJ): project psi(Zana) onto phi_vJ(r), shift the
// parallel momentum, FFT over (x, y) and pick the diffraction channels. With ProjOp = 2 and a
// fixed specular analysis surface, the coefficient of the initial channel is *replaced* by the
// projection of the specular wave at that surface (Fortran semantics).
// Fortran: Analyzor.f90 BalintKurtiCoeffs (the part before the S-matrix accumulation)
class BalintKurtiAnalyzer {
 public:
  BalintKurtiAnalyzer(const Config& config, const SpatialGrids& grids, const AnalysisBasis& basis, const GasPhase& gas,
                      const PropagationTables& tables, const Fft& fft);

  [[nodiscard]] BkCoefficients coefficients(const WavePacket& psi, const AnalysisSurface& surface) const;

 private:
  const AnalysisBasis& basis_;
  const GasPhase& gas_;
  const PropagationTables& tables_;
  int nx_;
  int ny_;
  int nr_;
  int iz_analysis_;   // Fortran: izanly - 1
  int j_in_;
  double dr_;         // Fortran: delr
  double xy_scale_;   // sqrt(dx dy sin(skew) / (nx ny))
  std::unique_ptr<ComplexFftPlan> xy_to_momentum_;  // scale 1. Fortran: yOverlap, xOverlap
};

}  // namespace spodvr
