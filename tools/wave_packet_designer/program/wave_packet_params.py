"""Absorber strength from the Vibok-Balint-Kurti table, and the Gaussian wave-packet parameters.

Shadows: src/input/wave_packet_params.cpp (absorber_reduced_length, optical_potential_strength, wave_packet_alpha,
wave_packet_average_energy)

Atomic units: energies in hartree, lengths in bohr, masses in electron masses. The operations are in the
C++ order, so the results are bitwise the program's.
"""

import math
from dataclasses import dataclass

# Table 3 of A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. 96, 8712 (1992): the optimal strength A / E_kin of
# a quadratic absorber against its length in wavelengths L / lambda (xtab and ftab of the C++).
VIBOK_BALINT_KURTI_TABLE = (
    # (L / lambda, A / E_kin)
    (1.0, 2.750),
    (1.5, 2.640),
    (2.0, 2.510),
    (2.5, 2.380),
    (3.0, 2.230),
    (3.5, 2.090),
    (4.0, 1.940),
    (4.5, 1.800),
    (5.0, 1.670),
    (5.5, 1.550),
    (6.0, 1.440),
    (6.5, 1.340),
    (7.0, 1.250),
    (7.5, 1.170),
    (8.0, 1.100),
    (8.5, 1.040),
    (9.0, 0.986),
    (9.5, 0.949),
    (10.0, 0.903),
    (10.5, 0.870),
    (11.0, 0.840),
    (11.5, 0.813),
    (12.0, 0.788),
    (12.5, 0.764),
    (13.0, 0.741),
    (13.5, 0.718),
    (14.0, 0.697),
    (14.5, 0.678),
    (15.0, 0.663),
)
LENGTH_TABLE = tuple(length for length, _ in VIBOK_BALINT_KURTI_TABLE)  # 1.0, 1.5, ..., 15.0
STRENGTH_TABLE = tuple(strength for _, strength in VIBOK_BALINT_KURTI_TABLE)

# Where L/lambda lies relative to the table (AbsorberTableRange).
BELOW, INSIDE, ABOVE, AT_FIRST_POINT = "below", "inside", "above", "at_first_point"


@dataclass(frozen=True)
class ReducedLength:
    """L / lambda of an absorber at its target energy, and where it lies in the table."""

    value: float
    range: str  # BELOW (< 1), INSIDE, ABOVE (> 15) or AT_FIRST_POINT (exactly 1: strength 0)


def absorber_reduced_length(kinetic_energy: float, length: float, mass: float) -> ReducedLength:
    """L / lambda for a wave of the given kinetic energy (hartree) and mass (electron masses)."""
    momentum = math.sqrt(2.0 * mass * kinetic_energy)
    wavelength = 2.0 * math.acos(-1.0) / momentum  # acos(-1) is PI bitwise; written as in the C++
    reduced_length = length / wavelength
    if reduced_length < LENGTH_TABLE[0]:
        return ReducedLength(reduced_length, BELOW)
    if reduced_length > LENGTH_TABLE[-1]:
        return ReducedLength(reduced_length, ABOVE)
    if reduced_length == LENGTH_TABLE[0]:  # no table interval xtab[i] < L/lambda <= xtab[i+1] contains it
        return ReducedLength(reduced_length, AT_FIRST_POINT)
    return ReducedLength(reduced_length, INSIDE)


def optical_potential_strength(kinetic_energy: float, reduced: ReducedLength) -> float:
    """Absorber strength A (hartree), linearly interpolated in the table and clamped at its ends.

    At L/lambda = 1 exactly no interval matches and the strength is 0: the absorber is off (absorber D2).
    """
    x = reduced.value
    if reduced.range == BELOW:
        return STRENGTH_TABLE[0] * kinetic_energy
    if reduced.range == ABOVE:
        return STRENGTH_TABLE[-1] * kinetic_energy
    for i in range(len(LENGTH_TABLE) - 1):
        x_low, x_high = LENGTH_TABLE[i], LENGTH_TABLE[i + 1]
        f_low, f_high = STRENGTH_TABLE[i], STRENGTH_TABLE[i + 1]
        if x_low < x <= x_high:
            # Multiply before dividing, as the C++ does (same rounding).
            return (f_low + (x - x_low) * (f_high - f_low) / (x_high - x_low)) * kinetic_energy
    return 0.0


def wave_packet_alpha(low_energy: float, high_energy: float, mass: float) -> float:
    """alpha of psi(Z) ~ exp(-alpha (Z - Z0)^2): the momenta of EMinWP and EMaxWP lie 2 sqrt(alpha) from k0."""
    k_min = math.sqrt(2.0 * mass * low_energy)
    k_max = math.sqrt(2.0 * mass * high_energy)
    k0 = (k_min + k_max) / 2.0
    return (k0 - k_min) ** 2 / 4.0


def wave_packet_average_energy(low_energy: float, high_energy: float, mass: float) -> float:
    """Kinetic energy k0^2 / 2M at the mean momentum k0 = (k_min + k_max) / 2.

    The momenta are negative (the packet moves towards the surface), as in the C++; the sign does not change
    k0^2.
    """
    k_min = -math.sqrt(2.0 * mass * low_energy)
    k_max = -math.sqrt(2.0 * mass * high_energy)
    k0 = (k_min + k_max) / 2.0
    return k0**2 / (2.0 * mass)
