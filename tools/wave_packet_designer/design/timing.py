"""The time budget of a run: how long the packet needs to reach the surface and come back out across Zana
(inputs/README.md section 8), as numbers and as the straight-line paths Z(t) of the time panel.

A classical estimate: each part of the packet moves at its normal speed k/M from Z0 to a turning point
Z_turn, where it is reflected, and back out across Zana. Reflected elastically it leaves at the same speed:
the slowest analysed part (the momentum of EMinWP) is then back at Zana after t_elastic, a lower bound for
tEnd. Rotational excitation and diffraction take normal energy, so those channels leave slower; near every
channel threshold arbitrarily slowly, so tEnd is chosen for a lowest outgoing normal energy E_n,min, which
needs t_needed. Not included: time spent in wells or trapped, and the spreading of the packet.
"""

import math
from dataclasses import dataclass

from design.rules import DesignInput, t_end_for


@dataclass(frozen=True)
class PathZ:
    """One straight-line path Z(t) through the given points (a.u., bohr)."""

    label: str
    times: tuple[float, ...]
    z: tuple[float, ...]


@dataclass(frozen=True)
class TimeBudget:
    """The times (a.u.) the slow parts of the packet need, and the paths of the time panel."""

    z_turn: float  # bohr
    lowest_outgoing_energy: float  # eV: E_n,min
    t_elastic: float  # the slowest analysed part (EMinWP), reflected elastically, back at Zana
    t_needed: float  # an outgoing wave of normal energy E_n,min, from the slowest part's reflection, at Zana
    t_end: float  # the input's tEnd
    paths: tuple[PathZ, ...]


def time_budget(c: DesignInput) -> TimeBudget:
    """The time budget of the input with the turning point and E_n,min of `c`."""
    mass, z0, z_analysis, z_turn = c.d.total_mass, c.p.z0, c.grids.z_analysis, c.z_turn
    width = 2.0 * math.sqrt(c.d.wave_packet_alpha)  # k_min and k_max lie 2 sigma_k from k0
    paths = []
    for label, k in (("EMaxWP", c.d.kz + width), ("mean k0", c.d.kz), ("EMinWP", c.d.kz - width)):
        t_in = (z0 - z_turn) * mass / k
        t_out = t_in + (z_analysis - z_turn) * mass / k
        paths.append(PathZ(f"{label}, elastic", (0.0, t_in, t_out), (z0, z_turn, z_analysis)))
    k_slow = c.d.kz - width
    t_turn = (z0 - z_turn) * mass / k_slow
    t_elastic = t_turn + (z_analysis - z_turn) * mass / k_slow
    t_needed = t_end_for(z0, z_analysis, mass, k_slow, c.e_n_min, z_turn)
    paths.append(PathZ(f"E_n = {c.e_n_min:g} eV out", (t_turn, t_needed), (z_turn, z_analysis)))
    return TimeBudget(z_turn, c.e_n_min, t_elastic, t_needed, c.p.t_end, tuple(paths))
