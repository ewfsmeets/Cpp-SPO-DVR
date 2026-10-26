#include "potential/potential_grid.hpp"

#include <cmath>

#include "core/constants.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// V(0, 0, Z, r, pi/2, 0) in hartree for every r; iz = nzpot marks the gas phase for grid-based
// PESs. Fortran: Potential.f90 setupTheAsymptoticPotential / asypot in switchPotentialToGasPhase
std::vector<double> gas_phase_potential(const Potential& potential, double z, const SpatialGrids& grids) {
  const auto& units = units::legacy;
  std::vector<double> v(static_cast<std::size_t>(grids.r.size()));
  for (int ir = 0; ir < grids.r.size(); ++ir)
    v[static_cast<std::size_t>(ir)] = units.ev_to_hartree(
        potential.value(0.0, 0.0, z, grids.r[ir], Constants::pi / 2.0, 0.0, grids.nz_potential, ir));
  return v;
}

// Fortran: Potential.f90 switchPotentialToGasPhase
void switch_to_gas_phase(PotentialGrid& grid, const Potential& potential, const Config& config,
                         const SpatialGrids& grids) {
  using fortran_compare::rge;
  using fortran_compare::rle;
  const double zsw_min = config.input.z_min_switch;
  const double zsw_max = config.input.z_max_switch;
  const std::vector<double> asymptotic =
      gas_phase_potential(potential, config.derived.z_asymptotic_unswitched, grids);

  const std::size_t nr = static_cast<std::size_t>(grids.r.size());
  const std::size_t block = grid.index.n0 * grid.index.n1;
  for (std::size_t iz = 0; iz < static_cast<std::size_t>(grids.nz_potential); ++iz) {
    const double z = grids.z_scattering[static_cast<int>(iz)];
    if (rle(z, zsw_min)) continue;
    const double arg = Constants::pi * (z - zsw_min) / (zsw_max - zsw_min);
    const double f = 0.5 + 0.5 * std::cos(arg);
    for (std::size_t ir = 0; ir < nr; ++ir) {
      const double va = asymptotic[ir];
      const auto block_values = std::span<double>(grid.values).subspan((ir + nr * iz) * block, block);
      for (double& v : block_values) v = rge(z, zsw_max) ? va : f * (v - va) + va;
    }
  }
}

// Fortran: Potential.f90 ImposeLimitOnPotential / ImposeLimitOnPotentialAsy
void impose_limits(std::span<double> v, const Config& config) {
  if (config.input.impose_v_min)
    for (auto& x : v)
      if (x < config.derived.v_min) x = config.derived.v_min;
  if (config.input.impose_v_max)
    for (auto& x : v)
      if (x > config.derived.v_max) x = config.derived.v_max;
}

}  // namespace

PotentialGrid compute_potential_grid(const Potential& potential, const Config& config,
                                     const SpatialGrids& grids, const AngularGrid& angular) {
  const auto& units = units::legacy;
  const double ev2au = units.ev_to_hartree(1.0);  // Fortran multiplies by 1/chaev here

  const auto n_theta = static_cast<std::size_t>(angular.n_theta());
  const auto n_phi = static_cast<std::size_t>(angular.n_phi());
  const auto nx = static_cast<std::size_t>(grids.x.size());
  const auto ny = static_cast<std::size_t>(grids.y.size());
  const auto nr = static_cast<std::size_t>(grids.r.size());
  const auto nz_pot = static_cast<std::size_t>(grids.nz_potential);

  PotentialGrid grid{Index3D{n_theta * n_phi, nx * ny, nz_pot * nr}, {}, {}};
  grid.values.resize(grid.index.size());

  const auto n_zr = static_cast<long>(nz_pot * nr);
#pragma omp parallel for schedule(static)
  for (long izr = 0; izr < n_zr; ++izr) {
    const auto uzr = static_cast<std::size_t>(izr);
    const auto iz = static_cast<int>(uzr / nr);
    const auto ir = static_cast<int>(uzr % nr);
    const double z = grids.z_scattering[iz];
    const double r = grids.r[ir];
    for (std::size_t ix = 0; ix < nx; ++ix) {
      for (std::size_t iy = 0; iy < ny; ++iy) {
        const double x = grids.x[static_cast<int>(ix)];
        const double y = grids.y[static_cast<int>(iy)];
        const std::size_t ixy = iy + ny * ix;
        for (std::size_t itheta = 0; itheta < n_theta; ++itheta) {
          for (std::size_t iphi = 0; iphi < n_phi; ++iphi) {
            const double theta = angular.thetas()[itheta];
            const double phi = angular.phis()[iphi];
            grid.values[grid.index(iphi + n_phi * itheta, ixy, uzr)] =
                potential.value(x, y, z, r, theta, phi, iz, ir) * ev2au;
          }
        }
      }
    }
  }

  if (config.input.use_switching) switch_to_gas_phase(grid, potential, config, grids);
  impose_limits(grid.values, config);

  grid.asymptotic = gas_phase_potential(potential, config.derived.z_asymptotic, grids);
  impose_limits(grid.asymptotic, config);
  return grid;
}

void premultiply_with_time_step(std::span<double> values, double time_step) {
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < values.size(); ++i) values[i] = -values[i] * time_step;
}

}  // namespace spodvr
