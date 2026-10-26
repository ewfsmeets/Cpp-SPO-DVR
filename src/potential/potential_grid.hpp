#pragma once

#include <span>
#include <vector>

#include "core/array_index.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential.hpp"

namespace spodvr {

// The potential on the DVR grid of the interaction region (Z < zasymp) in hartree, layout
// (itp, ixy, izr) with itp = iphi + nphi*itheta, ixy = iy + ny*ix, izr = ir + nr*iz
// (Fortran memory order, first index fastest), plus the gas-phase potential V(r).
// Fortran: Potential.f90 vdvr (before PremultiplyWithTimestep) and GlobalData vasy
struct PotentialGrid {
  Index3D index;                // extents (ntheta*nphi, nx*ny, nzpot*nr)
  std::vector<double> values;   // Fortran: vdvr
  std::vector<double> asymptotic;  // V(r) at Z = zasymp, theta = pi/2 (hartree). Fortran: vasy

  [[nodiscard]] double operator()(std::size_t itp, std::size_t ixy, std::size_t izr) const {
    return values[index(itp, ixy, izr)];
  }
};

// Evaluates the PES on the grid, applies the optional switching to the gas phase and the
// optional minimum/maximum. Parallelised with OpenMP over (Z, r); every value is computed
// independently, so the result does not depend on the number of threads.
// Fortran: Potential.f90 setupThePES/computeThePES (up to, not including, the -dt
// premultiplication), switchPotentialToGasPhase, ImposeLimitOnPotential(Asy),
// setupTheAsymptoticPotential
[[nodiscard]] PotentialGrid compute_potential_grid(const Potential& potential, const Config& config,
                                                   const SpatialGrids& grids,
                                                   const AngularGrid& angular);

// vdvr = -V * dt, the form used by the propagator. Fortran: Potential.f90 premultiplyWithTimestep
void premultiply_with_time_step(std::span<double> values, double time_step);

}  // namespace spodvr
