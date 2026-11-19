"""The analysis of a run as the program sets it up: the analysis energies, the analysis time step and the
time-independent S-matrix factors of every channel, from run_info.json and the BK file's channel energies
(the formulas: copies/s_matrix.py)."""

from dataclasses import dataclass

import numpy as np

from copies.s_matrix import channel_factors, normal_momentum
from readers.binary import BkCoefficients
from readers.run_info import RunInfo


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


def s_matrix_factors(info: RunInfo, energies: AnalysisEnergies, bk: BkCoefficients) -> np.ndarray:
    """prefactor * exp(i arg_f) * sqrt(k_f k0 / weight) of every analysis energy and channel (index
    i_state + n_states v): S = time sums * these. Returns: (n_energies, n_channels) complex."""
    packet = info.data["wave_packet"]
    return channel_factors(
        energies.total_ev,
        bk.channel_energy_hartree.reshape(-1),
        energies.k0_per_bohr,
        kz_per_bohr=packet["kz"],
        alpha_per_bohr2=packet["alpha"],
        total_mass_me=info.data["molecule"]["total_mass"],
        z_analysis_bohr=info.data["analysis"]["z_ana"],  # the grid point; the C++ uses the input ZAna (within 1e-12)
        z0_bohr=packet["z0"],
        hartree_to_ev=info.hartree_to_ev,
    )
