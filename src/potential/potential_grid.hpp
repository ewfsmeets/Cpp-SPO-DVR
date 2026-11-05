#pragma once

#include <span>
#include <vector>

#include "core/aligned_vector.hpp"
#include "core/array_index.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential.hpp"
#include "symmetry/grid_symmetry.hpp"

namespace spodvr {

// The potential on the DVR grid of the interaction region (Z < zasymp) in hartree, layout
// (itp, irep, izr) with itp = iphi + nphi*itheta, izr = ir + nr*iz (Fortran memory order, first
// index fastest), plus the gas-phase potential V(r). irep runs over the (x, y) representatives of
// the surface symmetry; without symmetry irep = ixy = iy + ny*ix.
// Fortran: Potential.f90 vdvr (before PremultiplyWithTimestep) and GlobalData vasy
struct PotentialGrid {
  Index3D index;                   // extents (ntheta*nphi, n_representatives, nzpot*nr)
  AlignedVector<double> values;    // Fortran: vdvr (first touch in parallel, see AlignedAllocator)
  std::vector<double> asymptotic;  // V(r) at Z = zasymp, theta = pi/2 (hartree). Fortran: vasy
  GridSymmetry symmetry;           // the (x, y) representatives and maps (identity without symmetry)

  // V at the full-grid point (itp, ixy, izr), through the symmetry maps.
  // Fortran: vdvr(symmetryIndex(...)) in PerformEXPPotentialOp
  [[nodiscard]] double operator()(std::size_t itp, std::size_t ixy, std::size_t izr) const {
    const int n_phi = symmetry.nphi();
    const int iphi = static_cast<int>(itp) % n_phi;
    const int itheta = static_cast<int>(itp) / n_phi;
    const int ixy_index = static_cast<int>(ixy);
    const int iphi_image = symmetry.phi_index(symmetry.operation_of(ixy_index))[iphi];
    return values[index(iphi_image + n_phi * itheta, symmetry.representative_of(ixy_index), izr)];
  }
};

// Evaluates the PES on the grid (at the (x, y) representatives of the symmetry), applies the
// optional switching to the gas phase and the optional minimum/maximum. Parallelised with OpenMP
// over (Z, r); every value is computed independently, so the result does not depend on the
// number of threads.
// Fortran: Potential.f90 setupThePES/computeThePES (up to, not including, the -dt
// premultiplication), switchPotentialToGasPhase, ImposeLimitOnPotential(Asy),
// setupTheAsymptoticPotential
[[nodiscard]] PotentialGrid compute_potential_grid(const Potential& potential, const Config& config,
                                                   const SpatialGrids& grids, const AngularGrid& angular,
                                                   GridSymmetry symmetry);

// The same with the symmetry of the run (make_grid_symmetry); with UseSymmetry the PES is first
// checked (check_potential_symmetry, require_symmetric_potential: throws InputError).
[[nodiscard]] PotentialGrid compute_potential_grid(const Potential& potential, const Config& config,
                                                   const SpatialGrids& grids, const AngularGrid& angular);

// The gas-phase potential V(r) alone (PotentialGrid::asymptotic, hartree), for runs that do not
// need the 6D grid (fine analysis). Fortran: Potential.f90 setupTheAsymptoticPotential
[[nodiscard]] std::vector<double> compute_asymptotic_potential(const Potential& potential, const Config& config,
                                                               const SpatialGrids& grids);

// vdvr = -V * dt, the form used by the propagator. Fortran: Potential.f90 premultiplyWithTimestep
void premultiply_with_time_step(std::span<double> values, double time_step);

}  // namespace spodvr
