#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/energy_resolved.hpp"
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
// analysis steps; with them the data they are computed from (writeDensities: the block densities of the
// scattering wave in izr order, the specular wave and the analysis-surface index; no extra work, the
// densities are computed once for both); the time and Z at which the dynamic specular analysis surface
// becomes fixed; the BK coefficients of every analysis step with the time at its end; with dtFrames, the backend
// at tStart and every analyses_per_frame analysis steps (the frames read it through z_plane and specular, after the
// analysis of that step); with NormBudget, the block densities and the specular wave at tStart and after every
// analysis step (the averages of those times then use the same densities: one reduction, not two); with the
// energy-resolved sums, those sums once, at the first analysis step at or after 0.8 tEnd (the partial sums that show
// the effect of stopping the time sums early).
struct PropagationObserver {
  std::function<void(const Averages&)> on_averages;
  std::function<void(const Averages&, std::span<const double> block_densities, const SpecularWave& specular,
                     int surface_index)>
      on_average_densities;
  std::function<void(double time, double z)> on_surface_fixed;
  std::function<void(double time, const BkCoefficients&)> on_analysis_step;
  std::function<void(double time, const PropagationBackend& backend, int surface_index)> on_frame;
  std::function<void(double time, std::span<const double> block_densities, const SpecularWave& specular)>
      on_analysis_densities;
  std::function<void(const EnergyResolvedSums&)> on_energy_resolved_partial;
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

  // NormBudget: the wall-clock time of the block densities of every analysis step and the observer's work on them
  // (the averages rows' own reductions included), for the cost line at the end of the run.
  [[nodiscard]] double analysis_densities_seconds() const { return analysis_densities_seconds_; }

  // Propagates from tStart to tEnd, analysing the wave packet after every analysis step.
  void run(const PropagationObserver& observer);

  // Restart: replace the wave packet (`read` fills it, see PropagationBackend::read_wave_packet),
  // continue with the saved analysis surface, and add the BK coefficients of the earlier runs to
  // the S-matrix.
  void read_wave_packet(const std::function<void(WavePacket&)>& read) { backend_->read_wave_packet(read); }
  void restore_analysis_surface(const AnalysisSurfaceState& state) { surface_.restore(state); }
  void add_earlier_analysis_step(double time, const BkCoefficients& bk) { s_matrix_.add(bk, time); }

  // EnergyResolvedSites/Energies: sum every analysis step from `start` on (empty_energy_resolved_sums for a new
  // run, the earlier segment's sums for a restart). Without this call nothing is summed and site_column never runs.
  void start_energy_resolved(EnergyResolvedSums start) { energy_resolved_.emplace(std::move(start)); }
  [[nodiscard]] bool sums_energy_resolved() const { return energy_resolved_.has_value(); }
  // The sums so far (only after start_energy_resolved).
  [[nodiscard]] const EnergyResolvedSums& energy_resolved_sums() const { return energy_resolved_.value().sums(); }
  // The wall-clock time of the sums (column copies and accumulation), for the cost line at the end of the run.
  [[nodiscard]] double energy_resolved_seconds() const { return energy_resolved_seconds_; }

  // S-matrix and probabilities from the analysis steps done so far.
  [[nodiscard]] ScatteringResults scattering_results() const { return s_matrix_.results(); }

 private:
  // The averages at `time` for the observer (and the densities they come from, when it asks for them), from
  // `densities` when they were computed already at this time (NormBudget), else from a new reduction.
  void report_averages(double time, const PropagationObserver& observer,
                       const std::optional<std::vector<double>>& densities) const;

  // NormBudget: the block densities of this time handed to the observer (which keeps nothing of them) and returned
  // for the averages of the same time; nothing without NormBudget (then no reduction runs).
  [[nodiscard]] std::optional<std::vector<double>> report_analysis_densities(double time,
                                                                             const PropagationObserver& observer);

  // Adds the analysis step ending at `time` to the energy-resolved sums, and hands the partial sums to the observer
  // at the first step at or after partial_sums_fraction tEnd (partial_written: done in this segment).
  void sum_energy_resolved(double time, const PropagationObserver& observer, bool& partial_written);

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
  std::optional<EnergyResolvedAccumulator> energy_resolved_;  // only with EnergyResolvedSites/Energies
  double analysis_densities_seconds_ = 0.0;  // NormBudget's time (analysis_densities_seconds)
  double energy_resolved_seconds_ = 0.0;     // the energy-resolved sums' time (energy_resolved_seconds)
};

}  // namespace spodvr
