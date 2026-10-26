#pragma once

#include "core/constants.hpp"

namespace spodvr {

// Unit conversion factors plus named conversion procedures.
//
// The Fortran code mixes two sets of factors (see CLAUDE.md): `legacy` (chaev, camau, cboan,
// chacm in constants.f90, used almost everywhere) and `precise` (Constants_eV2Ha, u2aum,
// Bohr2Ang, Ha2cmInv, used in FluxAnalyzor.f90). Each factor is stored in the direction in
// which the Fortran code defines it, and every conversion multiplies or divides by it exactly
// as the Fortran code does, so results are bitwise reproducible.
class UnitConversion {
 public:
  enum class EnergyFactorKind { HartreeToEv, EvToHartree };

  constexpr UnitConversion(double energy_factor, EnergyFactorKind energy_kind, double amu_to_me,
                           double bohr_to_angstrom, double hartree_to_wavenumber)
      : energy_factor_(energy_factor),
        energy_kind_(energy_kind),
        amu_to_me_(amu_to_me),
        bohr_to_angstrom_(bohr_to_angstrom),
        hartree_to_wavenumber_(hartree_to_wavenumber) {}

  [[nodiscard]] constexpr double ev_to_hartree(double x) const {
    return energy_kind_ == EnergyFactorKind::HartreeToEv ? x / energy_factor_ : x * energy_factor_;
  }
  [[nodiscard]] constexpr double hartree_to_ev(double x) const {
    return energy_kind_ == EnergyFactorKind::HartreeToEv ? x * energy_factor_ : x / energy_factor_;
  }
  [[nodiscard]] constexpr double amu_to_me(double x) const { return x * amu_to_me_; }
  [[nodiscard]] constexpr double me_to_amu(double x) const { return x / amu_to_me_; }
  [[nodiscard]] constexpr double bohr_to_angstrom(double x) const { return x * bohr_to_angstrom_; }
  [[nodiscard]] constexpr double angstrom_to_bohr(double x) const { return x / bohr_to_angstrom_; }
  [[nodiscard]] constexpr double hartree_to_wavenumber(double x) const { return x * hartree_to_wavenumber_; }
  [[nodiscard]] constexpr double wavenumber_to_hartree(double x) const { return x / hartree_to_wavenumber_; }

  // The raw factors, for printing and for tests.
  [[nodiscard]] constexpr double energy_factor() const { return energy_factor_; }
  [[nodiscard]] constexpr EnergyFactorKind energy_factor_kind() const { return energy_kind_; }

 private:
  double energy_factor_;
  EnergyFactorKind energy_kind_;
  double amu_to_me_;
  double bohr_to_angstrom_;
  double hartree_to_wavenumber_;
};

namespace units {

// Fortran: constants.f90 chaev, camau, cboan, chacm.
inline constexpr UnitConversion legacy{27.2114, UnitConversion::EnergyFactorKind::HartreeToEv, 1822.89, 0.529177,
                                       219475.0};

// Fortran: constants.f90 Constants_eV2Ha, Constants_u2aum, Constants_Bohr2Ang, Constants_Ha2cmInv.
inline constexpr UnitConversion precise{0.036749308286, UnitConversion::EnergyFactorKind::EvToHartree, 1822.88853006,
                                        0.529177249, 219475.0};

}  // namespace units

// Angle conversions (identical in both unit sets). Fortran: Constants_Deg2Rad.
[[nodiscard]] constexpr double deg_to_rad(double degrees) { return degrees * Constants::deg_to_rad; }
[[nodiscard]] constexpr double rad_to_deg(double radians) { return radians / Constants::deg_to_rad; }

}  // namespace spodvr
