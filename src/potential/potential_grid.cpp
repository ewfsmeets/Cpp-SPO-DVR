// The potential on the DVR grid of the interaction region and the gas-phase potential V(r).
// Fortran: Potential.f90
#include "potential/potential_grid.hpp"

#include <cmath>
#include <utility>

#include "core/constants.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "potential/potential_symmetry.hpp"

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
  const std::size_t block_size = grid.index.n0 * grid.index.n1;  // all (theta, phi) x representatives values
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

PotentialGrid compute_potential_grid(const Potential& potential, const Config& config, const SpatialGrids& grids,
                                     const AngularGrid& angular, GridSymmetry symmetry, int group_points) {
  const double ev2au = units::legacy.ev_to_hartree(1.0);  // Fortran multiplies by 1/chaev here

  const int n_theta = angular.n_theta();
  const int n_phi = angular.n_phi();
  const int ny = grids.y.size();
  const int nr = grids.r.size();
  const int n_zr = grids.nz_potential * nr;
  const int n_rep = symmetry.n_representatives();

  PotentialGrid grid{Index3D(n_theta * n_phi, n_rep, n_zr), {}, {}, std::move(symmetry)};
  grid.values.resize(grid.index.size());  // uninitialised; every element is written below

  // Every (Z, r) block is filled independently, so the result does not depend on the threads.
  // The (x, y) points are the representatives in memory order (all points without symmetry).
  // Written per (Z, r) group: first touch with the threads that use the group in the
  // propagation (schedule(runtime), core/loop_schedule.hpp).
#pragma omp parallel for schedule(runtime)
  for (int g = 0; g < n_zr / group_points; ++g)
    for (int izr = g * group_points; izr < (g + 1) * group_points; ++izr) {
      const int iz = izr / nr;
      const int ir = izr % nr;
      for (int irep = 0; irep < n_rep; ++irep) {
        const int ixy = grid.symmetry.representative_xy(irep);
        const int ix = ixy / ny;
        const int iy = ixy % ny;
        for (int itheta = 0; itheta < n_theta; ++itheta) {
          for (int iphi = 0; iphi < n_phi; ++iphi) {
            const double v = potential.value(grids.x[ix], grids.y[iy], grids.z_scattering[iz], grids.r[ir],
                                             angular.thetas()[itheta], angular.phis()[iphi], iz, ir);
            grid.values[grid.index(iphi + n_phi * itheta, irep, izr)] = v * ev2au;
          }
        }
      }
    }

  if (config.input.use_switching) switch_to_gas_phase(grid, potential, config, grids);
  impose_limits(grid.values, config);

  grid.asymptotic = compute_asymptotic_potential(potential, config, grids);
  return grid;
}

PotentialGrid compute_potential_grid(const Potential& potential, const Config& config, const SpatialGrids& grids,
                                     const AngularGrid& angular) {
  GridSymmetry symmetry = make_grid_symmetry(config, angular);
  if (config.input.use_symmetry)
    require_symmetric_potential(check_potential_symmetry(potential, symmetry, config, grids), config, symmetry);
  return compute_potential_grid(potential, config, grids, angular, std::move(symmetry));
}

std::vector<double> compute_asymptotic_potential(const Potential& potential, const Config& config,
                                                 const SpatialGrids& grids) {
  std::vector<double> v = gas_phase_potential(potential, config.derived.z_asymptotic, grids);
  impose_limits(v, config);
  return v;
}

void premultiply_with_time_step(std::span<double> values, double time_step) {
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < values.size(); ++i) values[i] = -values[i] * time_step;
}

}  // namespace spodvr
