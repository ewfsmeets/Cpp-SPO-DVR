// The energy-resolved sums at chosen sites and energies (EnergyResolvedSites, EnergyResolvedEnergies; C++ only).
// Their time-to-energy transform is the one of the on-the-fly S matrix (Fortran: Analyzor.f90 BalintKurtiCoeffs,
// "direct accumulation for some scattering energies"), applied to the wave packet instead of the BK coefficients.
#include "analysis/energy_resolved.hpp"

#include <cmath>
#include <utility>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "input/validate_outputs.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

int site_ixy_of(const InputParameters& p, const std::string& name) {
  for (const PlotSite& site : p.plot_sites)
    if (site.name == name) return plot_site_ixy(site, p.nx, p.ny);
  throw InputError("EnergyResolvedSites", "'" + name + "' is not a site of PlotSites");  // checked before
}

}  // namespace

std::size_t energy_resolved_memory_values(const InputParameters& p, int n_rot, int n_z_specular) {
  const std::size_t column = std::size_t(p.nz) * p.nr * n_rot;
  const std::size_t pairs = p.energy_resolved_sites.size() * p.energy_resolved_energies.size();
  return pairs * column + p.energy_resolved_energies.size() * std::size_t(n_z_specular) * p.nr + column;
}

EnergyResolvedSums empty_energy_resolved_sums(const Config& config, const GasPhase& gas, int n_rot, int n_z_specular,
                                              double t_start) {
  const InputParameters& p = config.input;
  const auto& units = units::legacy;
  EnergyResolvedSums s;
  s.site_names = p.energy_resolved_sites;
  for (const std::string& name : s.site_names) s.site_ixy.push_back(site_ixy_of(p, name));
  // E_tot as make_analysis_energies (s_matrix.cpp): the same operations, so a ρ_E energy equal to an analysis
  // energy has the same phase factors.
  for (const double collision : p.energy_resolved_energies) {
    const double total = collision + units.energy_factor() * gas.initial.energy;
    s.collision_ev.push_back(collision);
    s.total_hartree.push_back(total / units.energy_factor());
  }
  s.n_z = p.nz;
  s.n_r = p.nr;
  s.n_rot = n_rot;
  s.n_z_specular = n_z_specular;
  s.analysis_time_step = config.derived.steps_per_analysis * p.time_step;  // as SMatrixAccumulator's
  s.last_time = t_start;
  s.sums.assign(s.site_names.size() * s.collision_ev.size() * s.site_block(), cplx(0.0, 0.0));
  s.specular_sums.assign(s.collision_ev.size() * s.specular_block(), cplx(0.0, 0.0));
  return s;
}

EnergyResolvedAccumulator::EnergyResolvedAccumulator(EnergyResolvedSums start)
    : sums_(std::move(start)), column_(sums_.site_block()) {}

void EnergyResolvedAccumulator::add(double time, const ColumnReader& read_column, const SpecularWave& specular) {
  const std::size_t n_energies = sums_.collision_ev.size();
  const std::size_t n_zr = std::size_t(sums_.n_z) * sums_.n_r;
  const std::size_t n_rot = sums_.n_rot;
  std::vector<cplx> phase_factors(n_energies);
  for (std::size_t ien = 0; ien < n_energies; ++ien)
    phase_factors[ien] = std::exp(Constants::imaginary_unit * sums_.total_hartree[ien] * time);

  for (std::size_t isite = 0; isite < sums_.site_ixy.size(); ++isite) {
    read_column(sums_.site_ixy[isite], column_);
    const std::span<cplx> site = std::span(sums_.sums).subspan(isite * n_energies * sums_.site_block());
    // One owner per element: the result does not depend on the thread count (one parallel region per site).
#pragma omp parallel for schedule(static)
    for (std::size_t izr = 0; izr < n_zr; ++izr)
      for (std::size_t ien = 0; ien < n_energies; ++ien) {
        const std::span<cplx> block = site.subspan(ien * sums_.site_block(), sums_.site_block());
        for (std::size_t irot = 0; irot < n_rot; ++irot)
          block[irot + n_rot * izr] =
              block[irot + n_rot * izr] + column_[izr + n_zr * irot] * phase_factors[ien] * sums_.analysis_time_step;
      }
  }
  for (std::size_t ien = 0; ien < n_energies; ++ien) {
    const std::span<cplx> block =
        std::span(sums_.specular_sums).subspan(ien * sums_.specular_block(), sums_.specular_block());
    for (std::size_t k = 0; k < block.size(); ++k)
      block[k] = block[k] + specular.values[k] * phase_factors[ien] * sums_.analysis_time_step;
  }
  sums_.last_time = time;
  sums_.steps_summed += 1;
}

}  // namespace spodvr
