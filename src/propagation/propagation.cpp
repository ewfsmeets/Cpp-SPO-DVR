// Set-up and time loop of a new wave-packet propagation. Fortran: Main.f90
#include "propagation/propagation.hpp"

#include <utility>

namespace spodvr {

namespace {

ScatteringShape scattering_shape(const Config& config, const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  return {basis.size(), p.nx, p.ny, p.nr, p.nz};
}

std::vector<double> premultiplied(PotentialGrid potential, double time_step) {
  premultiply_with_time_step(potential.values, time_step);
  return std::move(potential.values);
}

}  // namespace

WavePacketPropagation::WavePacketPropagation(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, PotentialGrid potential,
                                             const GasPhase& gas, const Fft& fft)
    : config_(config),
      grids_(grids),
      basis_(angular, config.input.j_in, config.input.mj_in),
      transform_(angular, basis_, fft),
      tables_(make_propagation_tables(config, grids, basis_, potential.asymptotic)),
      premultiplied_potential_(premultiplied(std::move(potential), config.input.time_step)),
      kinetic_(scattering_shape(config, basis_), tables_, fft),
      potential_(scattering_shape(config, basis_), grids.nz_potential, basis_.initial_index(), tables_, transform_,
                 premultiplied_potential_),
      specular_(config.input.nr, grids.z_specular.size(), grids.nz_potential, basis_.initial_index(),
                *grids.z_specular_optical, grids.r_optical, tables_, fft),
      propagator_(config.derived.steps_per_analysis, kinetic_, potential_, specular_),
      initial_(make_initial_wave_packet(config, grids, gas, fft)),
      // The initial wave packet lives on the specular grid; the scattering grid starts empty.
      psi_{ScatteringWave(basis_.size(), config.input.nx * config.input.ny, config.input.nr * config.input.nz),
           initial_.specular},
      surface_(config, grids, psi_.specular) {
  // With precomputed exp(-i V dt) the real potential is no longer needed.
  if (potential_.uses_precomputed_exponentials()) {
    premultiplied_potential_.clear();
    premultiplied_potential_.shrink_to_fit();
  }
}

// Fortran: Main.f90 DO propStep = nanspo, nstep, nanspo
void WavePacketPropagation::run(const std::function<void(const Averages&)>& on_averages,
                                const std::function<void(double)>& on_surface_fixed) {
  const InputParameters& p = config_.input;
  const int n_per_analysis = config_.derived.steps_per_analysis;
  const int n_per_average = config_.derived.analyses_per_average;

  on_averages(compute_averages(p.t_end - p.n_steps * p.time_step, psi_, config_, grids_, surface_.index()));

  int analysis_count = 0;
  for (int step = n_per_analysis; step <= p.n_steps; step += n_per_analysis) {
    const double time = p.t_end + (step - p.n_steps) * p.time_step;
    propagator_.propagate_analysis_step(psi_);

    // Fortran: BalintKurtiCoeffs (only its FindNewAnalysisPoint part so far)
    if (const auto z_fixed = surface_.update(psi_.specular)) on_surface_fixed(*z_fixed);

    analysis_count = analysis_count + 1;
    if (analysis_count % n_per_average == 0)
      on_averages(compute_averages(time, psi_, config_, grids_, surface_.index()));
  }
}

}  // namespace spodvr
