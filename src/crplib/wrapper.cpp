// Unit-converting entry points of the CRP library. Fortran: CRP_Wrapper.f90 and the
// Calculate_CRP_*_PES_Array routines of CRP_C*v_MoleculeSurface.f90 / CRP_C*v_AtomSurface.f90.
#include "crplib/wrapper.hpp"

namespace crplib {

namespace {

Vec3 to_internal(const Vec3& p, Units units) {
  return {convert_length(p.x, units.length, LengthUnit::bohr), convert_length(p.y, units.length, LengthUnit::bohr),
          convert_length(p.z, units.length, LengthUnit::bohr)};
}

Vec3 force_from_internal(const Vec3& g, Units units) {
  return {convert_force(g.x, internal_units, units), convert_force(g.y, internal_units, units),
          convert_force(g.z, internal_units, units)};
}

}  // namespace

MoleculeSurfaceValue calculate_molecule_surface(const MoleculeSurface& pes, const Vec3& a, const Vec3& b, Units units) {
  const std::optional<MoleculeSurfaceValue> value = pes.evaluate(to_internal(a, units), to_internal(b, units));
  // Fortran: POT6D starts from V = 0 and F = 0 and returns them unchanged when switch is set.
  const MoleculeSurfaceValue v = value.value_or(MoleculeSurfaceValue{0.00, {0.00, 0.00, 0.00}, {0.00, 0.00, 0.00}});
  return {convert_energy(v.v, EnergyUnit::hartree, units.energy), force_from_internal(v.gradient_a, units),
          force_from_internal(v.gradient_b, units)};
}

AtomSurfaceValue calculate_atom_surface(const AtomSurface& pes, const Vec3& position, Units units) {
  const AtomSurfaceValue v = pes.evaluate(to_internal(position, units));
  return {convert_energy(v.v, EnergyUnit::hartree, units.energy), force_from_internal(v.gradient, units)};
}

}  // namespace crplib
