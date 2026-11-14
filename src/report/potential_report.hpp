#pragma once
// Set-up report of the PES: the CRP data set (SystemID CRPlibrary, and CRPDistrtSurface through its
// CRP) and the SCM block (CRPDistrtSurface). The Fortran prints only the libraries' own messages.

#include <ostream>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential.hpp"
#include "potential/scm/scm_diagnostics.hpp"
#include "potential/scm/scm_files.hpp"
#include "report/format.hpp"

namespace spodvr {

// The report of `potential` (nothing for a PES without one). For CRPDistrtSurface it computes the
// SCM diagnostics, prints them and then throws InputError if ΔV changes V_asy (SCM D3). In an
// analysis-only run only V_asy is checked (there is no 6D potential).
void report_potential(std::ostream& out, const Style& style, const Potential& potential, const Config& config,
                      const SpatialGrids& grids, const AngularGrid& angular);

// The SCM block: the files and the snapshot, the rescaling, the displacements, the lattice check
// (SCM D9), ΔV at V_asy and on its plane (SCM D3), the cell-edge jumps (SCM D10) and ΔV below the
// CRP's Zmin (SCM D15), with a warning for each check that calls for one.
void print_scm_potential(std::ostream& out, const Style& style, const ScmInput& input, const ScmSlab& slab,
                         const ScmDiagnostics& diagnostics);

}  // namespace spodvr
