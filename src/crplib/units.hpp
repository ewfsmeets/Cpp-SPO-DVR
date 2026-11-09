#pragma once
// Units of the CRP library: the unit codes of the data-set files and the library's own conversion
// factors. Fortran: UnitConversion.f90 (M. Wijzenbroek). These factors differ from the ones the
// SPO-DVR program uses (e.g. 27.2113961 eV/Hartree); they are used only inside the library.

namespace crplib {

// Length unit codes. Fortran: UNITS_BOHR, UNITS_ANGSTROM
enum class LengthUnit { bohr = 0, angstrom = 1 };
// Energy unit codes. Fortran: UNITS_HARTREE, UNITS_EV, UNITS_KCALMOL, UNITS_KJMOL, UNITS_JOULE
enum class EnergyUnit { hartree = 2, ev = 3, kcal_per_mol = 4, kj_per_mol = 5, joule = 6 };

// A set of units. Fortran: TYPE(Units) (its mass and time units are not used by the library)
struct Units {
  LengthUnit length;
  EnergyUnit energy;
  bool operator==(const Units&) const = default;
};

// The units of all internal computations (atomic units). Fortran: InternalUnits
inline constexpr Units internal_units{LengthUnit::bohr, EnergyUnit::hartree};

// The unit for a code of the data-set file CRPData; throws DataError for an unknown code.
LengthUnit length_unit_from_code(int code);
EnergyUnit energy_unit_from_code(int code);

// x converted between units, exactly as the Fortran: x * conversion, where conversion =
// (1 / factor of `from`) * (factor of `to`), and x unchanged when the units are equal.
// Fortran: Convert_Length_Units, Convert_Energy_Units
double convert_length(double x, LengthUnit from, LengthUnit to);
double convert_energy(double x, EnergyUnit from, EnergyUnit to);

// A derivative dV/dx converted: (x * energy conversion of 1) / (length conversion of 1).
// Fortran: Convert_Force_Units
double convert_force(double x, Units from, Units to);

}  // namespace crplib
