"""Reflection and transmission of an absorber, in one dimension (inputs/README.md, section 7).

A plane wave of energy E meets V = -i A (3/2) ((q - start)/L)^2 on [start, start + L]: R is what comes back,
T what passes the end (and wraps around on the periodic grid). Transfer matrix through thin slices with
constant V, integrated backwards from the transmitted wave, with the program's strength A. A low-Z absorber is
the mirror image for a wave from the other side: same R and T.
"""

import cmath
import math
from dataclasses import dataclass

from program.config import InputParameters
from program.unit_conversion import ev_to_hartree

LOWEST_CHANNEL_ENERGY = 0.05  # eV: the scattering absorber's range in section 7 ("R + T <= 4e-5 for 0.05-1.0 eV")
ENERGY_POINTS = 8  # energies per absorber, evenly spaced on a log scale
SLICE = 0.002  # bohr: halving it changes R + T of the DZ absorbers by at most 1.1e-4 relative (0.03-1.5 eV)


@dataclass(frozen=True)
class AbsorberResponse:
    """R and T of one absorber over the energies it must absorb."""

    name: str  # the key prefix, e.g. "ZSpOpt"
    energies: tuple[float, ...]  # eV
    reflection: tuple[float, ...]
    transmission: tuple[float, ...]
    total: tuple[float, ...]  # R + T: what the absorber fails to absorb


def response_window(
    name: str, p: InputParameters, available_energy: float, r_grid_energy: float
) -> tuple[float, float]:
    """The energies (eV) an absorber must absorb. The specular ones see only the incoming channel, elastically
    reflected at the packet's energies. The scattering ones see every channel, from LOWEST_CHANNEL_ENERGY up
    to the available energy (EMaxWP plus the initial state's internal and parallel energy). The r absorber
    sees dissociating fragments, up to the largest r kinetic energy the r grid holds (section 4: "pi/dR must
    hold sqrt(2 mu E_r) of the dissociating fragments")."""
    if name.startswith("ZSp"):
        return p.e_min_wave_packet, p.e_max_wave_packet
    if name == "ROpt":
        return LOWEST_CHANNEL_ENERGY, r_grid_energy
    return LOWEST_CHANNEL_ENERGY, available_energy


def absorber_mass(name: str, total_mass: float, reduced_mass: float) -> float:
    """The mass a wave in the absorber's coordinate has: the reduced mass in r, the total mass in Z."""
    return reduced_mass if name == "ROpt" else total_mass


def absorber_response(
    name: str, length: float, strength: float, mass: float, window: tuple[float, float]
) -> AbsorberResponse:
    """R and T at ENERGY_POINTS energies spread over `window` (eV)."""
    low, high = window
    energies = tuple(low * (high / low) ** (i / (ENERGY_POINTS - 1)) for i in range(ENERGY_POINTS))
    pairs = [reflection_transmission(ev_to_hartree(e), length, strength, mass) for e in energies]
    reflection, transmission = tuple(r for r, _ in pairs), tuple(t for _, t in pairs)
    total = tuple(r + t for r, t in pairs)
    return AbsorberResponse(name, energies, reflection, transmission, total)


def reflection_transmission(energy: float, length: float, strength: float, mass: float) -> tuple[float, float]:
    """R and T at kinetic energy `energy` (hartree) for an absorber of length (bohr) and strength A (hartree)."""
    k_free = math.sqrt(2.0 * mass * energy)
    n_slices = max(1, round(length / SLICE))
    h = length / n_slices
    # Beyond the end only the transmitted wave exp(i k q), amplitude 1; go back slice by slice to the start.
    psi, dpsi = complex(1.0), 1j * k_free
    for i in reversed(range(n_slices)):
        x = (i + 0.5) / n_slices  # slice centre, relative to the absorber
        potential = -1j * strength * 1.5 * x * x
        k = cmath.sqrt(2.0 * mass * (energy - potential))
        c, s = cmath.cos(k * h), cmath.sin(k * h)
        psi, dpsi = c * psi - s / k * dpsi, k * s * psi + c * dpsi
    # Split the wave at the start into incoming exp(i k q) and reflected exp(-i k q) parts.
    incoming = (psi + dpsi / (1j * k_free)) / 2.0
    reflected = (psi - dpsi / (1j * k_free)) / 2.0
    return abs(reflected / incoming) ** 2, 1.0 / abs(incoming) ** 2
