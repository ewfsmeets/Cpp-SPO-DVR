#pragma once

#include <functional>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/s_matrix.hpp"
#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/averages.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/specular_operator.hpp"
#include "propagation/split_operator.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Callbacks of the time loop (unset ones are skipped): the averages at tStart and every nave
// analysis steps; the Z of the dynamic specular analysis surface when it becomes fixed; the BK
// coefficients of every analysis step with the time at its end.
struct PropagationObserver {
  std::function<void(const Averages&)> on_averages;
  std::function<void(double z)> on_surface_fixed;
  std::function<void(double time, const BkCoefficients&)> on_analysis_step;
};

// A new wave-packet propagation (ProjOp = 2): sets up the initial wave packet and all
// operators, then runs the time loop of Main.f90 (nread = 0).
// Fortran: Main.f90 (SetupRDBasis ... SetupSplitOperator, DO propStep loop), Averages.f90,
// Analyzor.f90 BalintKurtiCoeffs
class WavePacketPropagation {
 public:
  // The potential grid is taken over and premultiplied by -dt. config, grids and gas must
  // outlive the propagation (they are referenced).
  WavePacketPropagation(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                        PotentialGrid potential, const GasPhase& gas, const Fft& fft);

  [[nodiscard]] const RotationalBasis& rotational_basis() const { return basis_; }
  [[nodiscard]] const InitialWavePacketInfo& initial_wave_packet() const { return initial_.info; }
  [[nodiscard]] const AnalysisSurface& analysis_surface() const { return surface_; }
  [[nodiscard]] const PropagationTables& tables() const { return tables_; }
  [[nodiscard]] const WavePacket& wave_packet() const { return psi_; }
  [[nodiscard]] const AnalysisBasis& analysis_basis() const { return analysis_basis_; }
  [[nodiscard]] const AnalysisEnergies& analysis_energies() const { return s_matrix_.energies(); }

  // Propagates from tStart to tEnd, analysing the wave packet after every analysis step.
  void run(const PropagationObserver& observer);

  // S-matrix and probabilities from the analysis steps done so far.
  [[nodiscard]] ScatteringResults scattering_results() const { return s_matrix_.results(); }

 private:
  const Config& config_;
  const SpatialGrids& grids_;
  RotationalBasis basis_;
  AngularTransform transform_;
  PropagationTables tables_;
  AlignedVector<double> premultiplied_potential_;  // -V dt. Fortran: vdvr
  ScatteringKineticOperator kinetic_;
  PotentialOperator potential_;
  SpecularOperator specular_;
  SplitOperatorPropagator propagator_;
  InitialWavePacket initial_;
  WavePacket psi_;
  AnalysisSurface surface_;
  AnalysisBasis analysis_basis_;
  BalintKurtiAnalyzer balint_kurti_;
  SMatrixAccumulator s_matrix_;
};

}  // namespace spodvr
