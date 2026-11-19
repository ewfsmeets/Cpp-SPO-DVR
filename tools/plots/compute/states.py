"""The final states of the scattered molecules (from s_matrix.tsv): P(v', j') per energy, P(mJ') per j', and
the final-state alignment A0(2) with the quantisation axis along the surface normal."""

from dataclasses import dataclass

import numpy as np

from compute.analysis import energy_of_rows
from compute.results import LEVEL_FLOOR
from readers.tables import Table


@dataclass
class StateToState:
    """Probabilities of the final rovibrational states per analysis energy (n_energies first)."""

    collision_energy_ev: np.ndarray  # (n_energies,)
    levels: list[tuple[int, int]]  # the (v', j') with P above LEVEL_FLOOR of the largest P_scattered, sorted
    p_vj: np.ndarray  # (n_energies, n_levels): summed over mJ' and the diffraction channels
    p_mj: dict[int, np.ndarray]  # j' -> (n_energies, 2 j' + 1) over mJ' = -j' .. j'; NaN at |mJ'| > MJMax
    alignment: dict[int, np.ndarray]  # j' >= 1 -> A0(2)(E); NaN where j' > MJMax or P(j') is below the floor
    mj_max: int  # the run's effective MJMax (an input MJMax < 1 means JMax): |mJ'| <= mj_max in the basis
    j_top: int  # the highest j' of the rotational basis


def alignment_a0(p_mj: np.ndarray, j: int) -> np.ndarray:
    """A0(2) = sum_mJ P(mJ) [3 mJ^2 - j(j+1)] / [j(j+1) sum_mJ P(mJ)] (dimensionless) for j >= 1; p_mj:
    (..., 2j + 1) over mJ = -j .. j. It lies between -1 (only mJ = 0: rotation in a plane containing the
    normal, "cartwheel") and (2j - 1) / (j + 1) (only |mJ| = j: rotation about the normal, "helicopter").
    NaN where the sum of P is 0."""
    m_j = np.arange(-j, j + 1)
    jj = j * (j + 1)
    total = np.sum(p_mj, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.sum(p_mj * (3 * m_j * m_j - jj), axis=-1) / (jj * total)


def helicopter_limit(j: int) -> float:
    """The largest A0(2) of a j: (2j - 1) / (j + 1)."""
    return (2 * j - 1) / (j + 1)


def summed(energy_of_row: np.ndarray, selected: np.ndarray, p: np.ndarray, n_energies: int) -> np.ndarray:
    """P of the selected rows summed per analysis energy."""
    return np.bincount(energy_of_row[selected], weights=p[selected], minlength=n_energies)


def mj_populations(table: Table, energy_of_row: np.ndarray, n_energies: int, mj_max: int) -> dict[int, np.ndarray]:
    """j' -> P(mJ') (n_energies, 2 j' + 1), summed over v' and the channels; NaN at |mJ'| > mj_max."""
    j, m_j = table["j"].astype(int), table["mJ"].astype(int)
    populations = {}
    for final_j in sorted(set(j.tolist())):
        p = np.full((n_energies, 2 * final_j + 1), np.nan)
        for k, value in enumerate(range(-final_j, final_j + 1)):
            if abs(value) <= mj_max:
                p[:, k] = summed(energy_of_row, (j == final_j) & (m_j == value), table["P"], n_energies)
        populations[final_j] = p
    return populations


def state_to_state(table: Table, collision_energy_ev: np.ndarray, mj_max: int, j_top: int) -> StateToState:
    """P(v', j'), P(mJ') and A0(2) of the final states in s_matrix.tsv (rows: energy, v, j, mJ, n, m, P);
    levels and alignments below LEVEL_FLOOR of P_scattered are rounding noise and left out."""
    energy_of_row = energy_of_rows(table, collision_energy_ev)
    n_energies = len(collision_energy_ev)
    v, j = table["v"].astype(int), table["j"].astype(int)
    p_scattered = summed(energy_of_row, np.ones(len(v), bool), table["P"], n_energies)
    floor = LEVEL_FLOOR * np.max(p_scattered)
    all_levels = sorted(set(zip(v.tolist(), j.tolist(), strict=True)))
    p_all = [summed(energy_of_row, (v == a) & (j == b), table["P"], n_energies) for a, b in all_levels]
    kept = [k for k, p in enumerate(p_all) if np.max(p) > floor]
    p_mj = mj_populations(table, energy_of_row, n_energies, mj_max)
    alignment = {}
    for final_j, p in p_mj.items():
        if final_j >= 1:
            noise = ~(np.nansum(p, axis=1) > LEVEL_FLOOR * p_scattered)  # also where P_scattered is 0
            alignment[final_j] = np.where(noise | (final_j > mj_max), np.nan, alignment_a0(p, final_j))
    levels = [all_levels[k] for k in kept]
    p_vj = np.column_stack([p_all[k] for k in kept]) if kept else np.zeros((n_energies, 0))
    return StateToState(collision_energy_ev, levels, p_vj, p_mj, alignment, mj_max, j_top)


def level_map(s: StateToState, energy_index: int) -> tuple[list[int], list[int], np.ndarray]:
    """P(v', j') at one energy as a grid: (the v' values, the j' values, grid (n_v', n_j'), NaN for a (v', j')
    without probability)."""
    v_values = sorted({v for v, _ in s.levels})
    j_values = sorted({j for _, j in s.levels})
    grid = np.full((len(v_values), len(j_values)), np.nan)
    for k, (v, j) in enumerate(s.levels):
        grid[v_values.index(v), j_values.index(j)] = s.p_vj[energy_index, k]
    return v_values, j_values, grid


def mj_fractions(s: StateToState, energy_index: int) -> dict[int, np.ndarray]:
    """j' -> P(mJ') / P(j') at one energy, for j' >= 1 whose alignment is defined there (not noise, not
    truncated); each sums to 1."""
    fractions = {}
    for j, alignment in s.alignment.items():
        if not np.isnan(alignment[energy_index]):
            p = s.p_mj[j][energy_index]
            fractions[j] = p / np.sum(p)
    return fractions
