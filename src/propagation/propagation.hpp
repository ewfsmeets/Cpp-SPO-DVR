#pragma once

#include <functional>
#include <memory>
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
#include "propagation/averages.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/propagation_backend.hpp"
#include "propagation/propagation_tables.hpp"
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

// A new wave-packet propagation (ProjOp = 2): sets up the initial wave packet, the tables and
// the backend that holds and propagates the wave packet (CPU or GPU), then runs the time loop of
// Main.f90 (nread = 0). The analysis (surface, BK coefficients, S-matrix, averages) runs on the
// host, on the pieces of the wave packet the backend hands out after every analysis step.
// Fortran: Main.f90 (SetupRDBasis ... SetupSplitOperator, DO propStep loop), Averages.f90,
// Analyzor.f90 BalintKurtiCoeffs
class WavePacketPropagation {
 public:
  // The potential grid is taken over by the backend (and premultiplied by -dt there). config,
  // grids, angular, gas and fft must outlive the propagation (they are referenced).
  WavePacketPropagation(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                        PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                        const BackendFactory& make_backend = make_cpu_backend);

  [[nodiscard]] const RotationalBasis& rotational_basis() const { return basis_; }
  [[nodiscard]] const InitialWavePacketInfo& initial_wave_packet() const { return initial_.info; }
  [[nodiscard]] const AnalysisSurface& analysis_surface() const { return surface_; }
  [[nodiscard]] const PropagationTables& tables() const { return tables_; }
  // The complete wave packet in host memory (checkpoint, dump; expensive for a GPU backend).
  [[nodiscard]] const WavePacket& wave_packet_on_host() { return backend_->wave_packet_on_host(); }
  [[nodiscard]] const AnalysisBasis& analysis_basis() const { return analysis_basis_; }
  [[nodiscard]] const AnalysisEnergies& analysis_energies() const { return s_matrix_.energies(); }

  // Propagates from tStart to tEnd, analysing the wave packet after every analysis step.
  void run(const PropagationObserver& observer);

  // Restart: replace the wave packet (`read` fills it, see PropagationBackend::read_wave_packet),
  // continue with the saved analysis surface, and add the BK coefficients of the earlier runs to
  // the S-matrix.
  void read_wave_packet(const std::function<void(WavePacket&)>& read) { backend_->read_wave_packet(read); }
  void restore_analysis_surface(const AnalysisSurfaceState& state) { surface_.restore(state); }
  void add_earlier_analysis_step(double time, const BkCoefficients& bk) { s_matrix_.add(bk, time); }

  // S-matrix and probabilities from the analysis steps done so far.
  [[nodiscard]] ScatteringResults scattering_results() const { return s_matrix_.results(); }

 private:
  [[nodiscard]] Averages averages(double time) const;

  // Declaration order = construction order (the backend references basis_ and tables_).
  const Config& config_;
  const SpatialGrids& grids_;
  RotationalBasis basis_;
  PropagationTables tables_;
  InitialWavePacket initial_;
  std::unique_ptr<PropagationBackend> backend_;
  AnalysisSurface surface_;
  AnalysisBasis analysis_basis_;
  BalintKurtiAnalyzer balint_kurti_;
  SMatrixAccumulator s_matrix_;
};

}  // namespace spodvr
