// The potential on the DVR grid of the interaction region and the gas-phase potential V(r).
// Fortran: Potential.f90
#include "potential/potential_grid.hpp"

#include <cmath>

#include "core/constants.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// Gas-phase potential V(0, 0, Z, r, pi/2, 0) in hartree for every r; iz = nzpot marks the gas phase for grid-based
// PESs. Fortran: Potential.f90 setupTheAsymptoticPotential / asypot in switchPotentialToGasPhase
std::vector<double> gas_phase_potential(const Potential& potential, double z, const SpatialGrids& grids) {
  std::vector<double> v(grids.r.size());
  for (int ir = 0; ir < grids.r.size(); ++ir)
    v[ir] = units::legacy.ev_to_hartree(
        potential.value(0.0, 0.0, z, grids.r[ir], Constants::pi / 2.0, 0.0, grids.nz_potential, ir));
  return v;
}

// Blends the potential into the gas-phase potential V(r) at the unswitched zasymp:
// V -> f(Z) (V - V_gas) + V_gas with a cosine switch f from 1 at ZminSwitch to 0 at ZmaxSwitch.
// Fortran: Potential.f90 switchPotentialToGasPhase
void switch_to_gas_phase(PotentialGrid& grid, const Potential& potential, const Config& config,
                         const SpatialGrids& grids) {
  using fortran_compare::rge;
  using fortran_compare::rle;
  const double zsw_min = config.input.z_min_switch;
  const double zsw_max = config.input.z_max_switch;
  const std::vector<double> v_gas = gas_phase_potential(potential, config.derived.z_asymptotic_unswitched, grids);

  const int nr = grids.r.size();
  const std::size_t block_size = grid.index.n0 * grid.index.n1;  // all (theta, phi, x, y) values
  for (int iz = 0; iz < grids.nz_potential; ++iz) {
    const double z = grids.z_scattering[iz];
    if (rle(z, zsw_min)) continue;
    const double f = 0.5 + 0.5 * std::cos(Constants::pi * (z - zsw_min) / (zsw_max - zsw_min));
    for (int ir = 0; ir < nr; ++ir) {
      const double va = v_gas[ir];
      const auto block = std::span<double>(grid.values).subspan((ir + nr * iz) * block_size, block_size);
      for (double& v : block) v = rge(z, zsw_max) ? va : f * (v - va) + va;
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
  const double ev2au = units::legacy.ev_to_hartree(1.0);  // Fortran multiplies by 1/chaev here

  const int n_theta = angular.n_theta();
  const int n_phi = angular.n_phi();
  const int nx = grids.x.size();
  const int ny = grids.y.size();
  const int nr = grids.r.size();
  const int n_zr = grids.nz_potential * nr;

  PotentialGrid grid{Index3D(n_theta * n_phi, nx * ny, n_zr), {}, {}};
  grid.values.resize(grid.index.size());

  // Every (Z, r) block is filled independently, so the result does not depend on the threads.
#pragma omp parallel for schedule(static)
  for (int izr = 0; izr < n_zr; ++izr) {
    const int iz = izr / nr;
    const int ir = izr % nr;
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
        for (int itheta = 0; itheta < n_theta; ++itheta) {
          for (int iphi = 0; iphi < n_phi; ++iphi) {
            const double v = potential.value(grids.x[ix], grids.y[iy], grids.z_scattering[iz], grids.r[ir],
                                             angular.thetas()[itheta], angular.phis()[iphi], iz, ir);
            grid.values[grid.index(iphi + n_phi * itheta, iy + ny * ix, izr)] = v * ev2au;
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
