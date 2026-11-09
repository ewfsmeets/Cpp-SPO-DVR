#pragma once
// The coefficient tables of the Fourier interpolation over the surface sites, copied verbatim from
// the Fortran DATA statements and scaled as the Fortran does. Element (ksit, j) of a table is the
// weight of site ksit in the coefficient of basis function j; stored as in Fortran (column-major:
// ksit + n * j, 0-based). Tests compare them bit for bit with the tables of the Fortran library.

#include <array>

namespace crplib::detail {

// Hexagonal 6D interpolation, 12 sites (Appendix A of docs/crp_library.md), used by C3v and C6v.
// Fortran: CRP_C6v_MoleculeSurface.f90 / CRP_C3v_MoleculeSurface.f90 FOUR6D (fouri6D)
const std::array<double, 144>& fourier_table_molecule_hexagonal();

// Atom-surface interpolation of C6v, 7 sites. Fortran: CRP_C6v_AtomSurface.f90 FOUR3D (fouri3D).
// Its DATA contains the truncated literals 0.666666666666667 and 1.333333333333: the table is the
// inverse of its site basis only to ~1e-13 (kept).
const std::array<double, 49>& fourier_table_atom_c6v();

// Atom-surface interpolation of C3v, 10 sites. Fortran: CRP_C3v_AtomSurface.f90 FOUR3D (fouri3D)
const std::array<double, 100>& fourier_table_atom_c3v();

}  // namespace crplib::detail
