// Set-up and time loop of a new wave-packet propagation. Fortran: Main.f90
#include "propagation/propagation.hpp"

#include <utility>

namespace spodvr {

WavePacketPropagation::WavePacketPropagation(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas,
                                             const Fft& fft, const BackendFactory& make_backend)
    : config_(config),
      grids_(grids),
      basis_(angular, config.input.j_in, config.input.mj_in),
      tables_(make_propagation_tables(config, grids, basis_, potential.asymptotic)),
      initial_(make_initial_wave_packet(config, grids, gas, fft)),
      backend_(make_backend(BackendInput{config, grids, angular, basis_, tables_, fft}, std::move(potential),
                            initial_.specular)),
      surface_(config, grids, initial_.specular),
      analysis_basis_(config, angular, basis_, gas, tables_.xy_kinetic_energies),
      balint_kurti_(config, grids, analysis_basis_, gas, tables_, fft),
      s_matrix_(config, analysis_basis_, gas) {}

Averages WavePacketPropagation::averages(double time) const {
  return compute_averages(time, backend_->block_densities(), backend_->specular(), config_, grids_, surface_.index());
}

// Fortran: Main.f90 DO propStep = nanspo, nstep, nanspo
void WavePacketPropagation::run(const PropagationObserver& observer) {
  const InputParameters& p = config_.input;
  const int n_per_analysis = config_.derived.steps_per_analysis;
  const int n_per_average = config_.derived.analyses_per_average;

  auto report_averages = [&](double time) {
    if (observer.on_averages) observer.on_averages(averages(time));
  };
  report_averages(p.t_end - p.n_steps * p.time_step);

  int analysis_count = 0;
  for (int step = n_per_analysis; step <= p.n_steps; step += n_per_analysis) {
    const double time = p.t_end + (step - p.n_steps) * p.time_step;
    backend_->propagate_analysis_step();

    // Analysis: move the specular analysis surface, project onto the asymptotic channels and
    // accumulate the S-matrix. Fortran: BalintKurtiCoeffs
    if (const auto z_fixed = surface_.update(backend_->specular()); z_fixed && observer.on_surface_fixed)
      observer.on_surface_fixed(*z_fixed);
    const BkCoefficients bk =
        balint_kurti_.coefficients(backend_->z_plane(grids_.iz_analysis), backend_->specular(), surface_);
    s_matrix_.add(bk, time);
    if (observer.on_analysis_step) observer.on_analysis_step(time, bk);

    analysis_count = analysis_count + 1;
    if (analysis_count % n_per_average == 0) report_averages(time);
  }
}

}  // namespace spodvr
