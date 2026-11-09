// Tests: the CRP library's unit conversion and constants (src/crplib/units.*, constants.hpp).
#include <gtest/gtest.h>

#include <bit>
#include <cmath>
#include <cstdint>

#include "crplib/constants.hpp"
#include "crplib/errors.hpp"
#include "crplib/units.hpp"

using namespace crplib;

namespace {

constexpr double conv_angstrom = 0.529177249;  // Fortran: CONV_ANGSTROM
constexpr double conv_ev = 27.2113961;         // Fortran: CONV_EV
constexpr double conv_kcal = 627.50956;        // Fortran: CONV_KCALMOL

// Values whose conversion is sensitive to the order of the operations.
constexpr double samples[] = {0.1, 0.3, 1.7, 2.9, 3.3, 7.1, 13.37, 0.0123, 255.5, 1e-5};

}  // namespace

TEST(CrpUnits, CodesOfTheDataSetFiles) {
  EXPECT_EQ(length_unit_from_code(0), LengthUnit::bohr);
  EXPECT_EQ(length_unit_from_code(1), LengthUnit::angstrom);
  EXPECT_EQ(energy_unit_from_code(2), EnergyUnit::hartree);
  EXPECT_EQ(energy_unit_from_code(3), EnergyUnit::ev);
  EXPECT_EQ(energy_unit_from_code(6), EnergyUnit::joule);
  EXPECT_THROW((void)length_unit_from_code(2), DataError);
  EXPECT_THROW((void)energy_unit_from_code(1), DataError);
  EXPECT_THROW((void)energy_unit_from_code(7), DataError);
}

TEST(CrpUnits, SameUnitsAreTheIdentity) {
  for (double x : samples) {
    EXPECT_EQ(convert_length(x, LengthUnit::angstrom, LengthUnit::angstrom), x);
    EXPECT_EQ(convert_energy(x, EnergyUnit::ev, EnergyUnit::ev), x);
    EXPECT_EQ(
        convert_force(x, Units{LengthUnit::angstrom, EnergyUnit::ev}, Units{LengthUnit::angstrom, EnergyUnit::ev}), x);
  }
}

// Fortran: Val = Num * Conversion with Conversion = 1/CONV_ANGSTROM (to bohr) — a multiplication
// by the reciprocal, which differs from a division for some values.
TEST(CrpUnits, LengthMultipliesByTheReciprocal) {
  bool differs_from_division = false;
  for (double x : samples) {
    EXPECT_EQ(convert_length(x, LengthUnit::angstrom, LengthUnit::bohr), x * (1.0 / conv_angstrom));
    EXPECT_EQ(convert_length(x, LengthUnit::bohr, LengthUnit::angstrom), x * conv_angstrom);
    if (x * (1.0 / conv_angstrom) != x / conv_angstrom) differs_from_division = true;
  }
  EXPECT_TRUE(differs_from_division) << "the samples do not distinguish the Fortran form";
}

// Fortran: Conversion = (1/CONV_from) * CONV_to is formed first, then Num * Conversion.
TEST(CrpUnits, EnergyFormsTheConversionFirst) {
  bool differs_from_left_to_right = false;
  for (double x : samples) {
    EXPECT_EQ(convert_energy(x, EnergyUnit::ev, EnergyUnit::hartree), x * (1.0 / conv_ev));
    EXPECT_EQ(convert_energy(x, EnergyUnit::hartree, EnergyUnit::ev), x * conv_ev);
    EXPECT_EQ(convert_energy(x, EnergyUnit::ev, EnergyUnit::kcal_per_mol), x * ((1.0 / conv_ev) * conv_kcal));
    if (x * ((1.0 / conv_ev) * conv_kcal) != (x * (1.0 / conv_ev)) * conv_kcal) differs_from_left_to_right = true;
  }
  EXPECT_TRUE(differs_from_left_to_right) << "the samples do not distinguish the Fortran form";
  EXPECT_EQ(convert_energy(1.0, EnergyUnit::hartree, EnergyUnit::kj_per_mol), 627.50956 * 4.184);
  EXPECT_EQ(convert_energy(1.0, EnergyUnit::hartree, EnergyUnit::joule), 4.35974434e-18);
}

// Fortran: Convert_Force_Units = Num * Convert_Energy(1) / Convert_Length(1), left to right.
TEST(CrpUnits, ForceIsEnergyOverLength) {
  const Units ev_angstrom{LengthUnit::angstrom, EnergyUnit::ev};
  for (double x : samples) {
    EXPECT_EQ(convert_force(x, ev_angstrom, internal_units), x * (1.0 / conv_ev) / (1.0 / conv_angstrom));
    EXPECT_EQ(convert_force(x, internal_units, ev_angstrom), x * conv_ev / conv_angstrom);
  }
}

// The constants gfortran folded at compile time equal the run-time values of this machine.
TEST(CrpConstants, FoldedTrigonometryEqualsRunTime) {
  EXPECT_EQ(std::bit_cast<std::uint64_t>(detail::cos_pi_3), 0x3FE0000000000001u);
  EXPECT_EQ(std::bit_cast<std::uint64_t>(detail::sin_pi_3), 0x3FEBB67AE8584CAAu);
  EXPECT_EQ(std::bit_cast<std::uint64_t>(detail::cos_pi_2), 0x3C91A62633145C07u);
  EXPECT_EQ(detail::cos_pi_3, std::cos(detail::pi / 3.0));
  EXPECT_EQ(detail::sin_pi_3, std::sin(detail::pi / 3.0));
  EXPECT_EQ(detail::cos_pi_2, std::cos(detail::pi / 2.0));
  EXPECT_EQ(detail::sin_pi_2, std::sin(detail::pi / 2.0));
}

TEST(CrpConstants, LiteralsOfTheFortran) {
  EXPECT_EQ(detail::pi, 3.141592653589793);
  EXPECT_EQ(detail::sq3, 1.7320508075688772);
  EXPECT_EQ(detail::usq3, 1.0 / 1.7320508075688772);
}
