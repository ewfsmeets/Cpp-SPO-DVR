// Set-up and time loop of a new wave-packet propagation. Fortran: Main.f90
#include "propagation/propagation.hpp"

#include <stdexcept>
#include <utility>

namespace spodvr {

namespace {

// The flux analysis of the input (doFluxAnalysis), else nothing. Its energies use the parallel
// energy (Fortran edifin) and the initial state's energy (eviroi); c_n weighs each step with the
// effective analysis step nanspo dt, not the input dtAna (flux Q4).
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor
std::optional<FluxAccumulator> make_flux_accumulator(const Config& config, const SpatialGrids& grids,
                                                     const RotationalBasis& basis, const PropagationTables& tables,
                                                     const GasPhase& gas, int slot_count) {
  if (!config.input.do_flux_analysis) return std::nullopt;
  const int n_blocks = basis.size() * config.input.nx * config.input.ny;
  return std::optional<FluxAccumulator>(std::in_place, n_blocks, grids.nz_potential,
                                        make_flux_energies(config, tables.parallel_energy, gas.initial.energy),
                                        config.derived.analysis_time_step, slot_count);
}

}  // namespace

WavePacketPropagation::WavePacketPropagation(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas,
                                             const Fft& fft, const BackendFactory& make_backend, int flux_slot_count)
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
      s_matrix_(config, analysis_basis_, gas),
      flux_(make_flux_accumulator(config, grids, basis_, tables_, gas, flux_slot_count)) {}

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
      observer.on_surface_fixed(time, *z_fixed);
    const BkCoefficients bk =
        balint_kurti_.coefficients(backend_->z_plane(grids_.iz_analysis), backend_->specular(), surface_);
    s_matrix_.add(bk, time);
    // After the BK analysis, as Main.f90 calls storeFluxCoeffs after BalintKurtiCoeffs.
    if (flux_) add_flux_step(time, observer);
    if (observer.on_analysis_step) observer.on_analysis_step(time, bk);

    analysis_count = analysis_count + 1;
    if (analysis_count % n_per_average == 0) report_averages(time);
  }
  // The steps still in the slots: from here on the coefficients are complete.
  if (flux_) flux_->flush();
}

// Fortran: FluxAnalyzor.f90 storeFluxCoeffs
void WavePacketPropagation::add_flux_step(double time, const PropagationObserver& observer) {
  flux_->add_step(time, [&](FluxSurfaceWave& slot) {
    backend_->flux_surface(slot);
    if (observer.on_flux_surface) observer.on_flux_surface(time, slot);
  });
}

void WavePacketPropagation::require_flux_analysis() const {
  if (!flux_) throw std::logic_error("WavePacketPropagation: no flux analysis in this run (doFluxAnalysis)");
}

FluxResults WavePacketPropagation::flux_results() const {
  require_flux_analysis();
  return spodvr::flux_results(flux_->coefficients(), flux_->energies(), make_flux_grid(config_, grids_, basis_.size()));
}

void WavePacketPropagation::read_flux_coefficients(const std::function<void(FluxCoefficients&)>& read) {
  require_flux_analysis();
  flux_->read_coefficients(read);
}

}  // namespace spodvr
