#pragma once
// The atom-surface PES of a CRP data set: V(x, y, z) of one atom and its gradient, the
// corrugation-reducing interpolation (a sum of 1D repulsive terms over the surface atoms plus a
// smooth interpolated remainder). Fortran: CRP_C3v_AtomSurface.f90, CRP_C6v_AtomSurface.f90 (ATII),
// CRP_C4v_AtomSurface.f90 (int_PES6CuHATOMCu100).

#include <memory>
#include <string>
#include <vector>

#include "crplib/dataset.hpp"

namespace crplib {

// A Cartesian position or gradient (bohr; Hartree/bohr): x, y along the surface (x along the first
// lattice vector), z along the surface normal (the top layer at z = 0).
struct Vec3 {
  double x, y, z;
};

// Value and gradient of the atom-surface PES (Hartree, Hartree/bohr).
struct AtomSurfaceValue {
  double v;
  Vec3 gradient;
};

class AtomSurface {
 public:
  virtual ~AtomSurface() = default;
  // V and its gradient at an atom position (internal units). Throws EvaluationError where the
  // Fortran stops (z outside the data, an atom too close to a surface atom).
  [[nodiscard]] virtual AtomSurfaceValue evaluate(const Vec3& position) const = 0;
  // Conditions the Fortran accepts but that make the PES questionable (cluster coverage).
  [[nodiscard]] virtual std::vector<std::string> warnings() const = 0;
};

// The atom-surface PES of the data set's symmetry: converts the data to internal units and builds
// the reduced cuts (data minus the repulsive sum), as Read_AtomSurface_Data. Throws DataError.
std::unique_ptr<const AtomSurface> make_atom_surface(const DataSet& data);

}  // namespace crplib
