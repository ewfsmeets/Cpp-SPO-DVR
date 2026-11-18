"""Case absorbers: L/lambda, the table ranges, the strengths and the set-up warnings (program/wave_packet_params.py,
process.py)."""

import math
from pathlib import Path

from program.config import AbsorberStrength
from program.potential_systems import ANALYTIC_SYSTEMS
from program.process import absorber_table_warnings, build_derived
from program.unit_conversion import ev_to_hartree
from program.wave_packet_params import (
    AT_FIRST_POINT,
    STRENGTH_TABLE,
    absorber_reduced_length,
    optical_potential_strength,
)
from tests.helpers import DZ_INPUT, HARDWALL_INPUT, check, check_equal, load_parameters

HARDWALL_SYSTEM = ANALYTIC_SYSTEMS["HardWallPES"]
DAI_LIGHT_SYSTEM = ANALYTIC_SYSTEMS["H2Cu111DaiLight"]
# A mass for the tests that need one: about the total mass of H2 (2.0156 amu = 3674.2 electron masses); the
# value itself does not matter, the tests build their lengths from the wavelength at this mass.
MASS = 3674.2


def test_fortran_printed_reduced_lengths(work: Path) -> None:
    """The input of TEST(Input, AbsorberOutsideTheStrengthTableIsRecorded): every energy-tuned absorber outside
    the table gets a warning, in ProcessInput's order, and L/lambda of ROpt and ZSpOpt equal the doubles Fortran
    5.24 printed (0.58847038336148938, 15.945765873846009), bitwise."""
    p = load_parameters(
        HARDWALL_INPUT,
        z_opt_energy=5.0,
        z_opt_low_start=-1.0,
        z_opt_low_end=-0.5,
        z_opt_low_energy=0.3,
        r_opt_energy=0.01,
        zsp_opt_energy=3.0,
        zsp_opt_low_start=-1.9,
        zsp_opt_low_end=-1.4,
        zsp_opt_low_energy=0.3,
    )
    absorbers = build_derived(p, HARDWALL_SYSTEM).absorbers
    check_equal([a.name for a in absorbers], ["ZOpt", "ZOptLow", "ROpt", "ZSpOpt", "ZSpOptLow"], "order")
    check_equal([a.table_range for a in absorbers], ["above", "below", "below", "above", "below"], "ranges")
    check_equal(absorbers[2].reduced_length, 0.58847038336148938, "L/lambda of ROpt")
    check_equal(absorbers[3].reduced_length, 15.945765873846009, "L/lambda of ZSpOpt")
    check_equal(absorbers[0].strength, STRENGTH_TABLE[-1] * ev_to_hartree(5.0), "ZOpt: the table's last value")
    check_equal(absorbers[1].strength, STRENGTH_TABLE[0] * ev_to_hartree(0.3), "ZOptLow: the table's first value")
    warnings = absorber_table_warnings(absorbers)
    keys = ["ZOptEnergy", "ZOptLowEnergy", "ROptEnergy", "ZSpOptEnergy", "ZSpOptLowEnergy"]
    check_equal([w.key for w in warnings], keys, "warnings")
    check("0.58847 wavelengths long" in warnings[2].message, warnings[2].message)


def test_exactly_one_wavelength_switches_the_absorber_off(work: Path) -> None:
    """L/lambda = 1 exactly matches no table interval: strength 0 and a warning (absorber D2)."""
    energy = ev_to_hartree(0.4)
    wavelength = 2.0 * math.acos(-1.0) / math.sqrt(2.0 * MASS * energy)
    reduced = absorber_reduced_length(energy, wavelength, MASS)
    check_equal((reduced.value, reduced.range), (1.0, AT_FIRST_POINT), "L/lambda")
    check_equal(optical_potential_strength(energy, reduced), 0.0, "strength")
    warnings = absorber_table_warnings((AbsorberStrength("ZOpt", 0.4, reduced.value, reduced.range, 0.0),))
    check(len(warnings) == 1 and "the absorber is off" in warnings[0].message, f"warnings {warnings}")


def test_interpolation_inside_the_table(work: Path) -> None:
    """L/lambda = 1.25 lies halfway between the first two entries: A = (2.750 + 0.25 (2.640 - 2.750) / 0.5) E
    = 2.695 E (by hand; L/lambda is 1.25 to rounding, so to 1e-14 relative)."""
    energy = ev_to_hartree(0.4)
    wavelength = 2.0 * math.acos(-1.0) / math.sqrt(2.0 * MASS * energy)
    strength = optical_potential_strength(energy, absorber_reduced_length(energy, 1.25 * wavelength, MASS))
    check(abs(strength - 2.695 * energy) <= 1e-14 * energy, f"strength {strength!r}, expected {2.695 * energy!r}")


def test_dz_strengths_match_the_program(work: Path) -> None:
    """The DZ input's three strengths, all interpolated inside the table, against "CAP strengths (Z, r, Zsp)" in
    the block INPUT of spodvr (provenance: helpers.DZ_INPUT), to half a printed unit."""
    absorbers = build_derived(load_parameters(DZ_INPUT), DAI_LIGHT_SYSTEM).absorbers
    printed = {"ZOpt": 1.27096140e-02, "ROpt": 2.45596439e-02, "ZSpOpt": 1.01936078e-02}
    for a in absorbers:
        check_equal(a.table_range, "inside", f"{a.name} range")
        half_unit = 0.5e-8 * 10 ** math.floor(math.log10(printed[a.name]))
        check(abs(a.strength - printed[a.name]) <= half_unit, f"{a.name}: {a.strength!r}, printed {printed[a.name]}")


def test_factor_and_zero_energy(work: Path) -> None:
    """A …Factor key is the strength; a target energy of 0 gives strength 0 (absorber D1)."""
    p = load_parameters(HARDWALL_INPUT, z_opt_energy=None, z_opt_factor=0.02)
    z_absorber = build_derived(p, HARDWALL_SYSTEM).absorbers[0]
    check_equal((z_absorber.name, z_absorber.strength), ("ZOpt", 0.02), "ZOptFactor")
    p = load_parameters(HARDWALL_INPUT, r_opt_energy=0.0)
    r_absorber = next(a for a in build_derived(p, HARDWALL_SYSTEM).absorbers if a.name == "ROpt")
    check_equal((r_absorber.strength, r_absorber.target_energy), (0.0, None), "ROptEnergy = 0")


TESTS = [
    test_fortran_printed_reduced_lengths,
    test_exactly_one_wavelength_switches_the_absorber_off,
    test_interpolation_inside_the_table,
    test_dz_strengths_match_the_program,
    test_factor_and_zero_energy,
]
