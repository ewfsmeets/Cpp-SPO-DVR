"""The analysis of a run as the program sets it up, from run_info.json and the BK file's state list: the
analysis energies, the analysis time step, the time-independent S-matrix factors of every channel (the
formulas: copies/s_matrix.py), and for each row of s_matrix.tsv its analysis energy and BK channel."""

from dataclasses import dataclass

import numpy as np

from copies.s_matrix import channel_factors, normal_momentum
from readers.binary import BkCoefficients
from readers.errors import ReadError
from readers.run_info import RunInfo
from readers.tables import Table


@dataclass
class AnalysisEnergies:
    """The analysis energies of a run: collision (eV), total = collision + E(v_in, j_in) (eV and hartree),
    and the normal momentum k0 of the collision energy (bohr^-1)."""

    collision_ev: np.ndarray
    total_ev: np.ndarray
    total_hartree: np.ndarray
    k0_per_bohr: np.ndarray


def analysis_energies(info: RunInfo) -> AnalysisEnergies:
    """The analysis energies of run_info.json, with k0 from the packet's parallel momentum and the lattice."""
    energies, packet = info.data["analysis_energies"], info.data["wave_packet"]
    collision_ev = np.array(energies["collision_ev"])
    total_ev = np.array(energies["total_ev"])
    k0 = normal_momentum(
        collision_ev,
        total_mass_me=info.data["molecule"]["total_mass"],
        hartree_to_ev=info.hartree_to_ev,
        kx_per_bohr=packet["kx"],
        ky_per_bohr=packet["ky"],
        skew_angle_rad=info.data["lattice"]["skew_angle_rad"],
    )
    return AnalysisEnergies(collision_ev, total_ev, total_ev / info.hartree_to_ev, k0)


def analysis_time_step_au(info: RunInfo) -> float:
    """The time between two BK records, as the program's S matrix takes it: steps per analysis * dtStep
    (src/analysis/s_matrix.cpp; dtAna itself need not be a whole multiple of dtStep)."""
    return info.data["times"]["steps_per_analysis"] * info.data["times"]["dt_step"]


def s_matrix_factors(
    info: RunInfo, energies: AnalysisEnergies, bk: BkCoefficients, lines_bohr: np.ndarray | None = None
) -> np.ndarray:
    """prefactor * exp(i arg_f) * sqrt(k_f k0 / weight) of every analysis energy and channel (index
    i_state + n_states v): S = time sums * these. lines_bohr (n_channels,): the line where each channel is
    measured (compute/time_series.py analysis_lines_bohr, as s_matrix.tsv); None: every channel at Zana, the
    Fortran's phases, enough where only |S| is used. Returns: (n_energies, n_channels) complex."""
    packet = info.data["wave_packet"]
    z_ana = float(info.data["analysis"]["z_ana"])  # the grid point; the C++ uses the input ZAna (within 1e-12)
    n_channels = bk.channel_energy_hartree.size
    lines = np.full(n_channels, z_ana) if lines_bohr is None else lines_bohr
    return channel_factors(
        energies.total_ev,
        bk.channel_energy_hartree.reshape(-1),
        energies.k0_per_bohr,
        kz_per_bohr=packet["kz"],
        alpha_per_bohr2=packet["alpha"],
        total_mass_me=info.data["molecule"]["total_mass"],
        z_measured_bohr=lines,
        z0_bohr=packet["z0"],
        hartree_to_ev=info.hartree_to_ev,
    )


ENERGY_MATCH_EV = 1e-9  # a row's E_coll and run_info.json's are the same double, printed with 16 digits


def energy_of_rows(table: Table, collision_energy_ev: np.ndarray) -> np.ndarray:
    """The index of each row's analysis energy. Raises ReadError for a row whose E_coll is not an analysis
    energy of the run (a table of another run)."""
    distinct, row_to_distinct = np.unique(table["E_coll"], return_inverse=True)
    index_of_distinct = np.empty(len(distinct), dtype=int)
    for k, energy in enumerate(distinct):
        nearest = int(np.argmin(np.abs(collision_energy_ev - energy)))
        if not abs(collision_energy_ev[nearest] - energy) <= ENERGY_MATCH_EV:
            raise ReadError(f"s_matrix.tsv: E_coll {energy} is not an analysis energy of the run; another run's file?")
        index_of_distinct[k] = nearest
    return index_of_distinct[row_to_distinct]


def channel_of_rows(table: Table, bk: BkCoefficients) -> np.ndarray:
    """Per row of s_matrix.tsv (v, j, mJ, n, m): its channel index i_state + n_states v in the BK file. Raises
    ReadError for a state the BK file does not have."""
    state_of = {tuple(int(q) for q in label[:4]): i for i, label in enumerate(bk.labels)}  # (j, mJ, n, m)
    n_states = len(bk.labels)
    channels = np.empty(len(table["v"]), dtype=int)
    for row in range(len(channels)):
        key = tuple(int(table[c][row]) for c in ("j", "mJ", "n", "m"))
        if key not in state_of:
            raise ReadError(f"s_matrix.tsv: state {key} is not in bk_coefficients.bin; files of different runs?")
        channels[row] = state_of[key] + n_states * int(table["v"][row])
    return channels


def channel_energy_of_rows(table: Table, bk: BkCoefficients) -> np.ndarray:
    """Per row of s_matrix.tsv: its channel's energy (hartree; internal plus parallel kinetic)."""
    return bk.channel_energy_hartree.reshape(-1)[channel_of_rows(table, bk)]
