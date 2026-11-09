#pragma once

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential.hpp"
#include "symmetry/grid_symmetry.hpp"

namespace spodvr {

// The surface symmetry used for the potential of a run: the PES's point group on the
// (x, y, phi) grid with UseSymmetry, else no symmetry. Throws InputError if UseSymmetry is set but
// no operation besides the identity is usable (point group "none", or nphi not a multiple the
// group needs; the message names the multiple, see nPhisOptimal).
// Fortran: CxvSurfaceSymmetry.f90 setupCxvSurfaceSymmetry (C6v only; here every group)
[[nodiscard]] GridSymmetry make_grid_symmetry(const Config& config, const AngularGrid& angular);

// Result of comparing the PES at sample configurations q with its value at the symmetry images
// g q (not in the Fortran). The images are wrapped back into the unit cell by a lattice vector,
// as the potential grid uses them, so the check covers the PES's periodicity too.
struct SymmetryCheck {
  int n_samples = 0;             // configurations q
  int n_comparisons = 0;         // (q, g) pairs, g != identity
  double max_deviation = 0.0;    // max |V(q) - V(g q)| (eV), after VMin/VMax
  double tolerance = 0.0;        // SymmetryTolerance (eV)
  // The worst configuration (skewed x, y; Z, r, theta, phi) and operation.
  double x = 0.0, y = 0.0, z = 0.0, r = 0.0, theta = 0.0, phi = 0.0;
  SurfaceOperation operation;
  [[nodiscard]] bool passed() const { return max_deviation <= tolerance; }
};

// Samples: (x, y, phi) from a fixed low-discrepancy sequence over the cell and the circle; theta,
// Z and r at grid points (theta is invariant under every surface operation, and a PES may only be
// defined where the grid reaches, e.g. the atom range of a CRP data set; Z and r so that grid-based
// PESs get valid iz, ir).
// VMin/VMax are applied to both values (the stored potential is clipped); the switching to the
// gas phase can only reduce a difference and is ignored.
[[nodiscard]] SymmetryCheck check_potential_symmetry(const Potential& potential, const GridSymmetry& symmetry,
                                                     const Config& config, const SpatialGrids& grids,
                                                     const AngularGrid& angular, int n_samples = 1000);

// Throws InputError (parameter "SymmetryTolerance") naming SystemID, point group, the worst
// configuration and operation if the check failed.
void require_symmetric_potential(const SymmetryCheck& check, const Config& config, const GridSymmetry& symmetry);

}  // namespace spodvr
