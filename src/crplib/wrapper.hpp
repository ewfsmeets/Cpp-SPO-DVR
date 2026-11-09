#pragma once
// The CRP library's general entry points with unit conversion, as the Fortran library offers them:
// positions in, V and gradient out, in the caller's units. Fortran: CRP_Wrapper.f90
// Calculate_CRP_MoleculeSurface_PES / Calculate_CRP_AtomSurface_PES (the *_Array variants; their
// R and F arrays become named fields here: R = (Za, Xa, Ya, Zb, Xb, Yb) but F = (dXa, dYa, dZa, dXb,
// dYb, dZb) in the Fortran).

#include "crplib/atom_surface.hpp"
#include "crplib/molecule_surface.hpp"
#include "crplib/units.hpp"

namespace crplib {

// V and gradient of the molecule at atom positions a and b, everything in `units`. Where the
// Fortran gives up silently (switch), V = 0 and the gradient is 0, as there.
MoleculeSurfaceValue calculate_molecule_surface(const MoleculeSurface& pes, const Vec3& a, const Vec3& b, Units units);

// V and gradient of one atom, everything in `units`.
AtomSurfaceValue calculate_atom_surface(const AtomSurface& pes, const Vec3& position, Units units);

}  // namespace crplib
