"""The gas phase: the rovibrational levels on the asymptotic potential, and which of them are bound."""

from dataclasses import dataclass

import numpy as np

from readers.tables import AsymptoticPotential, GasPhaseStates


@dataclass
class Levels:
    """The levels E(v, j) of the gas phase (eV) and whether each is bound (below V_asy at the last r point)."""

    v: np.ndarray
    j: np.ndarray
    energy_ev: np.ndarray
    bound: np.ndarray  # bool


def levels(states: GasPhaseStates, potential: AsymptoticPotential, hartree_to_ev: float) -> Levels:
    """Every state's energy in eV, and bound = E < V_asy(r_end) (the potential at the grid's last point)."""
    energy = states.energy_hartree * hartree_to_ev
    return Levels(states.v, states.j, energy, energy < potential.v_ev[-1])
