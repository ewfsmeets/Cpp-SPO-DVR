#include "grids/spatial_grids.hpp"

#include <cmath>

#include "core/error.hpp"
#include "core/real_compare.hpp"

namespace spodvr {

namespace {

using fortran_compare::rgt;
using fortran_compare::rlt;

// Points inside [start, end] (with the Fortran tolerant comparisons).
// Fortran: Setup.f90 nzchsp/nzchmx loops (and their nschsp, nrchsp counterparts)
IndexRange optical_range(const UniformGrid& grid, double start, double end) {
  int first = 0;
  while (first < grid.size() && rlt(grid[first], start)) ++first;
  int last = grid.size() - 1;
  while (last >= 0 && rgt(grid[last], end)) --last;
  return {first, last};
}

}  // namespace

SpatialGrids make_spatial_grids(const Config& config) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;

  // The Fortran code fills zgrid(1:nzsp) and uses zgrid(1:nz) as the scattering grid; it does
  // not check nZSpecular >= nZ (out-of-bounds access otherwise).
  if (p.nz_specular < p.nz)
    throw InputError("nZSpecular", "nZSpecular < nZ (the specular grid must contain the scattering grid)");

  SpatialGrids g{
      UniformGrid(p.z_start, p.dz, p.nz),
      UniformGrid(p.z_start, p.dz, p.nz_specular),
      UniformGrid(p.r_start, p.dr, p.nr),
      UniformGrid(0.0, d.dx, p.nx),
      UniformGrid(0.0, d.dy, p.ny),
      MomentumGrid(p.nz, p.dz, d.total_mass),
      MomentumGrid(p.nz_specular, p.dz, d.total_mass),
      MomentumGrid(p.nr, p.dr, d.reduced_mass),
      0, 0.0, 0, {0, 0}, std::nullopt, {0, 0}, std::nullopt};

  // Analysis point: first Z grid point >= ZAna; ZAna is moved onto that point.
  int iz = 0;
  while (iz < g.z_scattering.size() && !fortran_compare::rge(g.z_scattering[iz], p.z_analysis)) ++iz;
  if (iz == g.z_scattering.size())
    throw InputError("Zana", "analysis value ZAnalysis lies outside grid");
  g.iz_analysis = iz;
  g.z_analysis = g.z_scattering[iz];
  if (rgt(g.z_analysis, g.z_scattering.back()))
    throw InputError("Zana", "analysis value ZAnalysis lies outside grid");
  if (rgt(g.z_analysis, p.z_opt_start))
    throw InputError("Zana", "analysis value ZAnalysis lies inside optical potential");

  // Flux analysis surface in r.
  if (p.do_flux_analysis) {
    const int ir = static_cast<int>(std::lround((p.r_flux_analysis - p.r_start) / p.dr));
    if (ir < 0 || ir >= g.r.size() || fortran_compare::rne(g.r[ir], p.r_flux_analysis))
      throw InputError("rFluxAnalysis", "rFluxAnalysis does not coincide with a grid point on the r-grid.");
    g.ir_flux_analysis = ir;
  }

  // Number of Z points in the interaction region (Z < zasymp; zasymp itself is gas phase).
  iz = 0;
  while (iz < g.z_scattering.size() && !fortran_compare::rge(g.z_scattering[iz], d.z_asymptotic)) ++iz;
  if (iz == g.z_scattering.size())
    throw InputError("nZ", "the scattering grid in Z does not reach ZAsymptotic");
  g.nz_potential = iz;

  g.z_optical = optical_range(g.z_scattering, p.z_opt_start, p.z_opt_end);
  if (p.projection_operator == 2)
    g.z_specular_optical = optical_range(g.z_specular, p.zsp_opt_start, p.zsp_opt_end);
  g.r_optical = p.nr == 1 ? IndexRange{0, 0} : optical_range(g.r, p.r_opt_start, p.r_opt_end);
  return g;
}

}  // namespace spodvr
