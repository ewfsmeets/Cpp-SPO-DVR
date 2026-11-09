#pragma once
// Parts of the molecule-surface PES visible to their tests. Fortran: CRP_General.f90 CAMBCORD,
// RRINTII/SISII/FOUR6II (CRP_C3v/C6v_MoleculeSurface.f90), RRINT/SIS8x8 (CRP_C4v_MoleculeSurface.f90).

#include <array>
#include <span>

#include "crplib/atom_surface.hpp"
#include "crplib/dataset.hpp"
#include "crplib/spline.hpp"

namespace crplib::detail {

// Centre (ccdef = 0.5), bond length and orientation of the molecule from its atom positions.
// Throws EvaluationError for coinciding atoms. Fortran: CAMBCORD
struct MoleculeGeometry {
  double x, y, z, r, theta, phi;
};
MoleculeGeometry molecule_geometry(const Vec3& a, const Vec3& b);

// The interpolation function at the molecule's centre (x, y) and orientation (theta, phi) from the
// values of the cuts at (r, Z), and its derivatives. Fortran: Vinter(1..7) = value, d/dZ, d/dr,
// d/dphi, d/dtheta, d/dX, d/dY
using Interaction = std::array<double, 7>;

// C3v / C6v: the angular interpolation at the 12 Fourier sites (with their rotated images) and the
// Fourier interpolation over (x, y). Want::value: Vinter(1) only (S3). Fortran: the site part of
// RRINTII, SISII, FOUR6II
Interaction hexagonal_interaction(Symmetry symmetry, std::span<const BicubicValue> cuts, double theta, double phi,
                                  double x, double y, double lattice_constant, Want want = Want::value_and_gradient);
// C4v: 8 Fourier sites. Fortran: the site part of RRINT, SIS8x8
Interaction square_interaction(std::span<const BicubicValue> cuts, double theta, double phi, double x, double y,
                               double lattice_constant, Want want = Want::value_and_gradient);

// The 12 basis functions of the hexagonal 6D interpolation and their Cartesian x, y derivatives
// (Want::value: f only). Fortran: FOUR6II
void hexagonal_molecule_basis(double x, double y, double lattice_constant, std::span<double> f, std::span<double> dfdx,
                              std::span<double> dfdy, Want want = Want::value_and_gradient);

}  // namespace crplib::detail
