"""Unit conversions with the factors the program uses (the Fortran's "legacy" set).

Shadows: src/core/unit_conversion.hpp (UnitConversion, units::legacy, rad_to_deg, angstrom_to_bohr),
src/core/constants.hpp (pi)

The operations are those of the C++ (eV -> hartree divides by 27.2114, it does not multiply by the
inverse), so the converted values are bitwise the program's.
"""

HARTREE_IN_EV = 27.2114  # Fortran: chaev
AMU_IN_ELECTRON_MASSES = 1822.89  # Fortran: camau
BOHR_IN_ANGSTROM = 0.529177  # Fortran: cboan

PI = 3.1415926535897932384626433832795028841971  # Constants::pi (the same double as math.pi)
DEG_TO_RAD = PI / 180.0  # Constants::deg_to_rad


def rad_to_deg(radians: float) -> float:
    """Degrees from radians: divides by DEG_TO_RAD, as rad_to_deg in unit_conversion.hpp does."""
    return radians / DEG_TO_RAD


def ev_to_hartree(energy_ev: float) -> float:
    """Energy in hartree from eV (UnitConversion::ev_to_hartree)."""
    return energy_ev / HARTREE_IN_EV


def hartree_to_ev(energy_hartree: float) -> float:
    """Energy in eV from hartree (UnitConversion::hartree_to_ev)."""
    return energy_hartree * HARTREE_IN_EV


def angstrom_to_bohr(length_angstrom: float) -> float:
    """Length in bohr from angstrom (UnitConversion::angstrom_to_bohr: divides by BOHR_IN_ANGSTROM)."""
    return length_angstrom / BOHR_IN_ANGSTROM


def amu_to_electron_masses(mass_amu: float) -> float:
    """Mass in electron masses from atomic mass units (UnitConversion::amu_to_me)."""
    return mass_amu * AMU_IN_ELECTRON_MASSES
