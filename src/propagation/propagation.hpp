#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/flux_analysis.hpp"
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
// analysis steps; the time and Z at which the dynamic specular analysis surface becomes fixed; the BK
// coefficients of every analysis step with the time at its end; with the flux analysis, psi and
// d psi / dr at its surface r = r_c in every analysis step, with the time at its end.
struct PropagationObserver {
  std::function<void(const Averages&)> on_averages;
  std::function<void(double time, double z)> on_surface_fixed;
  std::function<void(double time, const BkCoefficients&)> on_analysis_step;
  std::function<void(double time, const FluxSurfaceWave&)> on_flux_surface;
};

// A new wave-packet propagation (ProjOp = 2): sets up the initial wave packet, the tables and
// the backend that holds and propagates the wave packet (CPU or GPU), then runs the time loop of
// Main.f90 (nread = 0). The analysis (surface, BK coefficients, S-matrix, averages, flux analysis)
// runs on the host, on the pieces of the wave packet the backend hands out after every analysis step.
// Fortran: Main.f90 (SetupRDBasis ... SetupSplitOperator, DO propStep loop), Averages.f90,
// Analyzor.f90 BalintKurtiCoeffs, FluxAnalyzor.f90 storeFluxCoeffs
class WavePacketPropagation {
 public:
  // The potential grid is taken over by the backend (and premultiplied by -dt there). config,
  // grids, angular, gas and fft must outlive the propagation (they are referenced). With the flux
  // analysis (doFluxAnalysis), flux_slot_count (>= 1) analysis steps are collected per pass over its
  // coefficients: MemoryEstimate::flux_slot_count. The results do not depend on it, the time does
  // (DZ: the flush takes about 5 times longer per step with 1 slot than with 8, docs/cuda.md).
  // Ignored without the flux analysis.
  WavePacketPropagation(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                        PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                        const BackendFactory& make_backend, int flux_slot_count);

  [[nodiscard]] const RotationalBasis& rotational_basis() const { return basis_; }
  [[nodiscard]] const InitialWavePacketInfo& initial_wave_packet() const { return initial_.info; }
  [[nodiscard]] const AnalysisSurface& analysis_surface() const { return surface_; }
  [[nodiscard]] const PropagationTables& tables() const { return tables_; }
  // The complete wave packet in host memory (checkpoint, dump; expensive for a GPU backend).
  [[nodiscard]] const WavePacket& wave_packet_on_host() { return backend_->wave_packet_on_host(); }
  [[nodiscard]] const AnalysisBasis& analysis_basis() const { return analysis_basis_; }
  [[nodiscard]] const AnalysisEnergies& analysis_energies() const { return s_matrix_.energies(); }

  // Propagates from tStart to tEnd, analysing the wave packet after every analysis step. Every step
  // is in the flux coefficients when it returns: nothing is pending outside run().
  void run(const PropagationObserver& observer);

  // Restart: replace the wave packet (`read` fills it, see PropagationBackend::read_wave_packet),
  // continue with the saved analysis surface, and add the BK coefficients of the earlier runs to
  // the S-matrix.
  void read_wave_packet(const std::function<void(WavePacket&)>& read) { backend_->read_wave_packet(read); }
  void restore_analysis_surface(const AnalysisSurfaceState& state) { surface_.restore(state); }
  void add_earlier_analysis_step(double time, const BkCoefficients& bk) { s_matrix_.add(bk, time); }

  // S-matrix and probabilities from the analysis steps done so far.
  [[nodiscard]] ScatteringResults scattering_results() const { return s_matrix_.results(); }

  // The flux analysis (energies, slot count, coefficients of the steps done so far, including those
  // read on a restart); empty without doFluxAnalysis. Its coefficients() throw std::logic_error
  // while steps wait in an unflushed slot: possible only inside run() (e.g. in an observer).
  [[nodiscard]] const std::optional<FluxAccumulator>& flux() const { return flux_; }
  // The reaction probability from the flux coefficients. Throws std::logic_error without
  // doFluxAnalysis, and while steps are pending (inside run(), as flux()).
  // Fortran: FluxAnalyzor.f90 computeFluxProbability
  [[nodiscard]] FluxResults flux_results() const;
  // Restart: `read` fills the flux coefficients of the earlier run in place (FluxAccumulator::
  // read_coefficients); call it before run(). Throws std::logic_error without doFluxAnalysis.
  void read_flux_coefficients(const std::function<void(FluxCoefficients&)>& read);

 private:
  [[nodiscard]] Averages averages(double time) const;
  // psi and d psi / dr at the flux surface after the analysis step that ends at `time`, into the
  // flux coefficients. Fortran: FluxAnalyzor.f90 storeFluxCoeffs
  void add_flux_step(double time, const PropagationObserver& observer);
  void require_flux_analysis() const;  // throws std::logic_error without doFluxAnalysis

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
  // After backend_: allocated once the backend has freed the potential (not in the set-up peak).
  std::optional<FluxAccumulator> flux_;
};

}  // namespace spodvr
