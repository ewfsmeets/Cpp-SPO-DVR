// V(r, Z) of the PES at one fixed (x, y, theta, phi) for write2DPES.
// Fortran: Potential.f90 Write2DPESToFile, H2SurfaceSystem.f90 write2DPESH2SurfaceSystem
#include "potential/pes_cut.hpp"

#include "core/unit_conversion.hpp"

namespace spodvr {

PesCut compute_pes_cut(const Potential& potential, const InputParameters& input, const SpatialGrids& grids) {
  PesCut cut;
  cut.x = input.x_2d_pes;
  cut.y = input.y_2d_pes;
  cut.theta_deg = input.theta_2d_pes;
  cut.phi_deg = input.phi_2d_pes;
  // The Fortran multiplies the degrees by Constants_Deg2Rad (write2DPESH2SurfaceSystem).
  cut.theta = deg_to_rad(input.theta_2d_pes);
  cut.phi = deg_to_rad(input.phi_2d_pes);
  cut.r.assign(grids.r.points().begin(), grids.r.points().end());
  cut.z.assign(grids.z_scattering.points().begin(), grids.z_scattering.points().end());

  // r outer, Z inner, as the Fortran writes them; the 0-based grid indices are passed as the PES's
  // ZIndex and RIndex (no ported PES uses them).
  const int nr = grids.r.size();
  const int nz = grids.z_scattering.size();
  cut.v_ev.resize(static_cast<std::size_t>(nr) * nz);
  for (int ir = 0; ir < nr; ++ir)
    for (int iz = 0; iz < nz; ++iz)
      cut.v_ev[iz + nz * ir] = potential.value(cut.x, cut.y, cut.z[iz], cut.r[ir], cut.theta, cut.phi, iz, ir);
  return cut;
}

}  // namespace spodvr
