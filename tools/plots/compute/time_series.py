"""Quantities against time: the Balint-Kurti coefficients per channel and the convergence of the scattering
probability in time."""

from dataclasses import dataclass

import numpy as np

from compute.analysis import analysis_energies, analysis_time_step_au, s_matrix_factors
from copies.s_matrix import time_sums
from readers.binary import BkCoefficients
from readers.run_info import RunInfo


def channel_densities_per_bohr(bk: BkCoefficients) -> np.ndarray:
    """|c_f(t)|^2 of every channel (n_times, n_v * n_states), in bohr^-1 (channel index i_state + n_states v)."""
    n_times, n_v, n_states = bk.coefficients.shape
    return np.abs(bk.coefficients.reshape(n_times, n_v * n_states)) ** 2


def channel_label(bk: BkCoefficients, channel: int) -> str:
    """The quantum numbers of a channel index (i_state + n_states v): "v=0 j=2 mJ=1 (n,m)=(1,0)"."""
    n_states = bk.labels.shape[0]
    j, m_j, n, m, _order = bk.labels[channel % n_states]
    return f"v={channel // n_states} j={j} mJ={m_j} (n,m)=({n},{m})"


def strongest_channels(densities: np.ndarray, count: int) -> np.ndarray:
    """The indices of the `count` channels with the largest |c|^2 at any time, strongest first."""
    order = np.argsort(np.max(densities, axis=0))[::-1]
    return order[:count]


@dataclass
class ScatteringInTime:
    """P_scattered(E; t): the probability the S matrix gives when the time sums stop at t."""

    times_au: np.ndarray  # (n_times,)
    collision_energy_ev: np.ndarray  # (n_energies,)
    p_scattered: np.ndarray  # (n_times, n_energies)


def scattering_in_time(info: RunInfo, bk: BkCoefficients) -> ScatteringInTime:
    """sum_f |S_f(E; t)|^2 at every BK time t: the time sums stopped at t, times the program's prefactor and
    phases. The sums grow record by record (as the C++ adds them), so memory stays that of one S matrix."""
    energies = analysis_energies(info)
    factors = s_matrix_factors(info, energies, bk)  # (n_energies, n_channels), time-independent
    dt_analysis_au = analysis_time_step_au(info)
    n_times, n_v, n_states = bk.coefficients.shape
    coefficients = bk.coefficients.reshape(n_times, n_v * n_states)
    sums = np.zeros(factors.shape, dtype=complex)
    p_scattered = np.empty((n_times, len(energies.collision_ev)))
    for record, time_au in enumerate(bk.times_au):
        phases = np.exp(1j * energies.total_hartree * time_au)  # (n_energies,)
        sums += coefficients[record][None, :] * phases[:, None] * dt_analysis_au
        p_scattered[record] = np.sum(np.abs(sums * factors) ** 2, axis=1)
    return ScatteringInTime(bk.times_au.copy(), energies.collision_ev, p_scattered)


def p_scattered_stopped_at(info: RunInfo, bk: BkCoefficients, last_record: int) -> np.ndarray:
    """P_scattered(E) (n_energies,) with the time sums stopped after record `last_record` (inclusive)."""
    energies = analysis_energies(info)
    n_times, n_v, n_states = bk.coefficients.shape
    coefficients = bk.coefficients[: last_record + 1].reshape(last_record + 1, n_v * n_states)
    sums = time_sums(bk.times_au[: last_record + 1], coefficients, energies.total_hartree, analysis_time_step_au(info))
    return np.sum(np.abs(sums * s_matrix_factors(info, energies, bk)) ** 2, axis=1)


def convergence_ratio(s: ScatteringInTime) -> np.ndarray:
    """P_scattered(E; t) / P_scattered(E; tEnd), (n_energies, n_times): 1 where the time sums have converged;
    NaN where P_scattered(E; tEnd) is 0."""
    final = s.p_scattered[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(final > 0, s.p_scattered / final, np.nan).T
