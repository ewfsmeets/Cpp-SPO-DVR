"""The results against energy: scattering and reaction probabilities, P(v'), the mean final energies and
the unitarity error."""

from dataclasses import dataclass

import numpy as np

from readers.tables import Table

MEAN_ENERGY_FLOOR = 1e-6  # heuristic: below this P_scattered a mean final energy is not drawn
LEVEL_FLOOR = 1e-6  # heuristic: a level is scattered into when its P(v') exceeds this fraction of max P_scattered


@dataclass
class Results:
    """The results per analysis energy (energies in eV, probabilities dimensionless)."""

    collision_energy_ev: np.ndarray
    p_scattered: np.ndarray
    p_reaction: np.ndarray
    p_vibration: dict[int, np.ndarray]  # P(v') of the levels scattered into; empty when fewer than two are
    mean_energy_ev: dict[str, np.ndarray]  # "vibration", "rotation", "diffraction"; NaN below MEAN_ENERGY_FLOOR
    unitarity_error: float  # -min(P_reaction) when negative, else 0; NaN when P_reaction holds a NaN


def scattered_levels(summary: Table) -> dict[int, np.ndarray]:
    """P(v') of every level whose P(v') exceeds LEVEL_FLOOR of the largest P_scattered at some energy."""
    levels = sorted(int(c.removeprefix("P_v")) for c in summary.columns if c.startswith("P_v"))
    largest = np.max(summary["P_scattered"])
    return {v: summary[f"P_v{v}"] for v in levels if np.max(summary[f"P_v{v}"]) > LEVEL_FLOOR * largest}


def results(summary: Table) -> Results:
    """The results of scattering_summary.tsv, with the masks the panels need."""
    p_scattered = summary["P_scattered"]
    reliable = p_scattered >= MEAN_ENERGY_FLOOR
    mean_energy_ev = {
        "vibration": np.where(reliable, summary["E_vib_mean"], np.nan),
        "rotation": np.where(reliable, summary["E_rot_mean"], np.nan),
        "diffraction": np.where(reliable, summary["E_dif_mean"], np.nan),
    }
    lowest = float(np.min(summary["P_reaction"]))  # NaN stays NaN, so the verdict warns
    unitarity_error = 0.0 if lowest >= 0 else -lowest
    levels = scattered_levels(summary)
    p_vibration = levels if len(levels) > 1 else {}  # one level: P(v') is P_scattered
    return Results(summary["E_coll"], p_scattered, summary["P_reaction"], p_vibration, mean_energy_ev, unitarity_error)


def matched_differences(fine: Table, on_the_fly: Table, column: str) -> tuple[np.ndarray, np.ndarray]:
    """(E_coll, fine - on-the-fly) of a column at the analysis energies both tables have (the fine analysis
    leaves out energies outside its grid); energies matched to 1e-9 eV (both print the same E_coll)."""
    energies, differences = [], []
    for row, energy in enumerate(fine["E_coll"]):
        match = np.flatnonzero(np.abs(on_the_fly["E_coll"] - energy) <= 1e-9)
        if match.size:
            energies.append(energy)
            differences.append(fine[column][row] - on_the_fly[column][match[0]])
    return np.array(energies), np.array(differences)
