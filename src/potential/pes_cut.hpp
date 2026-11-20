#pragma once
// A 2D cut of the PES for write2DPES: V(r, Z) at one fixed (x, y, theta, phi).
// Fortran: Potential.f90 Write2DPESToFile, H2SurfaceSystem.f90 write2DPESH2SurfaceSystem

#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential.hpp"

namespace spodvr {

// V of the PES itself (no switching to the gas phase, no VMin/VMax limits) at every point of the r grid and of
// the scattering Z grid (all nZ points, also those at Z >= zasymp), at the input's x2DPES, y2DPES (bohr, along
// the lattice vectors), theta2DPES and phi2DPES (degrees).
struct PesCut {
  double x;          // bohr. Fortran: x2DPES
  double y;          // bohr. Fortran: y2DPES
  double theta_deg;  // Fortran: theta2DPES
  double phi_deg;    // Fortran: phi2DPES
  double theta;      // rad: theta_deg * deg_to_rad. Fortran: theta_fixed
  double phi;        // rad. Fortran: phi_fixed
  std::vector<double> r;     // bohr, the r grid (nR)
  std::vector<double> z;     // bohr, the scattering Z grid (nZ)
  std::vector<double> v_ev;  // V in eV, (iz, ir) at iz + nZ * ir: Z fastest, the Fortran's loop order

  [[nodiscard]] double v(int ir, int iz) const { return v_ev[iz + z.size() * ir]; }
};

// The cut of `potential` at the write2DPES keys of `input` on the grids. Fortran: Potential.f90 Write2DPESToFile
[[nodiscard]] PesCut compute_pes_cut(const Potential& potential, const InputParameters& input,
                                     const SpatialGrids& grids);

}  // namespace spodvr
