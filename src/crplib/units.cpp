// Unit conversion of the CRP library. Fortran: UnitConversion.f90
#include "crplib/units.hpp"

#include <array>
#include <string>

#include "crplib/errors.hpp"

namespace crplib {

namespace {

// Conversion factors, verbatim. Fortran: CONV_ANGSTROM (bohr -> Angstrom), CONV_EV, CONV_KCALMOL,
// CONV_KJMOL, CONV_JOULE (Hartree -> unit)
constexpr double conv_angstrom = 0.529177249;
constexpr double conv_ev = 27.2113961;
constexpr double conv_kcal_per_mol = 627.50956;
constexpr double conv_kj_per_mol = 627.50956 * 4.184;  // a Fortran PARAMETER: rounded once
constexpr double conv_joule = 4.35974434e-18;

// Factor from Hartree to each energy unit code, and its reciprocal (the Fortran's 1.0 / CONV_*,
// folded at compile time there as here). Tables indexed by the code, covering every value of the
// 3-bit code so that no index can be out of range.
// A switch with a run-time division here was miscompiled by GCC 15 with LTO at -O2/-O3: the
// division by the switch table was hoisted above the Hartree check and read the table at index -1
// (crash in the C4v atom-surface constructor; fine with -fno-tree-switch-conversion).
constexpr std::array<double, 8> hartree_to = {1.0,        1.0, 1.0, conv_ev, conv_kcal_per_mol, conv_kj_per_mol,
                                              conv_joule, 1.0};
constexpr std::array<double, 8> to_hartree = {
    1.0, 1.0, 1.0, 1.0 / conv_ev, 1.0 / conv_kcal_per_mol, 1.0 / conv_kj_per_mol, 1.0 / conv_joule, 1.0};

int code(EnergyUnit unit) { return static_cast<int>(unit) & 7; }

// The Fortran's Conversion variable: 1 or 1/factor(from), then times factor(to) unless `to` is
// the internal unit.
double length_conversion(LengthUnit from, LengthUnit to) {
  constexpr double angstrom_to_bohr = 1.0 / conv_angstrom;
  double conversion = from == LengthUnit::bohr ? 1.0 : angstrom_to_bohr;
  if (to == LengthUnit::angstrom) conversion = conversion * conv_angstrom;
  return conversion;
}

double energy_conversion(EnergyUnit from, EnergyUnit to) {
  double conversion = to_hartree[code(from)];
  if (to != EnergyUnit::hartree) conversion = conversion * hartree_to[code(to)];
  return conversion;
}

}  // namespace

LengthUnit length_unit_from_code(int code) {
  if (code == 0) return LengthUnit::bohr;
  if (code == 1) return LengthUnit::angstrom;
  throw DataError("unknown length unit code " + std::to_string(code) + " (0 = bohr, 1 = Angstrom)");
}

EnergyUnit energy_unit_from_code(int code) {
  if (code >= 2 && code <= 6) return static_cast<EnergyUnit>(code);
  throw DataError("unknown energy unit code " + std::to_string(code) +
                  " (2 = Hartree, 3 = eV, 4 = kcal/mol, 5 = kJ/mol, 6 = Joule)");
}

// Fortran: Convert_Length_Units_Int
double convert_length(double x, LengthUnit from, LengthUnit to) {
  if (from == to) return x;
  return x * length_conversion(from, to);
}

// Fortran: Convert_Energy_Units_Int
double convert_energy(double x, EnergyUnit from, EnergyUnit to) {
  if (from == to) return x;
  return x * energy_conversion(from, to);
}

// Fortran: Convert_Force_Units = Num * Convert_Energy_Units(1.0) / Convert_Length_Units(1.0)
double convert_force(double x, Units from, Units to) {
  return x * convert_energy(1.0, from.energy, to.energy) / convert_length(1.0, from.length, to.length);
}

}  // namespace crplib
