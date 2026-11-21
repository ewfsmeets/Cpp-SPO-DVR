"""Quantities of the Balint-Kurti coefficients c_f(t): |c_f|^2 per channel, P_scattered when the time sums
stop at t (the convergence in time), and the time delays d arg S_f / dE with the line where each channel is
analysed."""

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


CHANNEL_FLOOR = 1e-6  # heuristic: a channel whose largest value is below this fraction of all is noise


def strongest_channels(values: np.ndarray, count: int) -> np.ndarray:
    """The indices of the `count` columns of values (n_rows, n_columns) with the largest value anywhere,
    strongest first, leaving out columns below CHANNEL_FLOOR of the largest (rounding noise)."""
    largest = np.max(values, axis=0)
    order = np.argsort(largest)[::-1]
    return np.array([k for k in order[:count] if largest[k] >= CHANNEL_FLOOR * largest[order[0]]], dtype=int)


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


def convergence_deviation(s: ScatteringInTime) -> np.ndarray:
    """|P_scattered(E; t) / P_scattered(E; tEnd) - 1|, (n_energies, n_times): how far the time sums stopped at
    t are from their final value; NaN where P_scattered(E; tEnd) is 0."""
    return np.abs(convergence_ratio(s) - 1.0)


def convergence_ratio(s: ScatteringInTime) -> np.ndarray:
    """P_scattered(E; t) / P_scattered(E; tEnd), (n_energies, n_times): 1 where the time sums have converged;
    NaN where P_scattered(E; tEnd) is 0."""
    final = s.p_scattered[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(final > 0, s.p_scattered / final, np.nan).T


def analysis_lines_bohr(info: RunInfo, bk: BkCoefficients, surface_index: int, surface_fixed: bool) -> np.ndarray:
    """Per channel (i_state + n_states v): the Z (bohr) where the program analyses it: Zana, except every v of
    the initial state (j_in, mJ_in, (0, 0)), analysed at the specular surface ZStart + index * dZ once that
    surface is fixed (src/analysis/balint_kurti.cpp; index and fixed from checkpoint.bin). Before it is fixed
    the program records those channels at Zana, where they hold no flux (~1e-32): their records are noise."""
    n_v, n_states = bk.channel_energy_hartree.shape
    lines = np.full(n_v * n_states, float(info.data["analysis"]["z_ana"]))
    if not surface_fixed:
        return lines
    state, grids = info.data["initial_state"], info.data["grids"]
    initial = [i for i, label in enumerate(bk.labels) if tuple(label[:4]) == (state["j"], state["mJ"], 0, 0)]
    for v in range(n_v):
        lines[[i + n_states * v for i in initial]] = grids["z_start"] + surface_index * grids["dz"]
    return lines


DELAY_FLOOR = 1e-3  # heuristic: a delay is shown where P_f is at least this fraction of the energy's P_scattered


@dataclass
class TimeDelay:
    """Time delays and probabilities per analysis energy and channel (n_energies first)."""

    collision_energy_ev: np.ndarray
    tau_au: np.ndarray  # (n_energies, n_channels): NaN for closed channels
    probability: np.ndarray  # (n_energies, n_channels): |S|^2


def time_delays(info: RunInfo, bk: BkCoefficients, lines_bohr: np.ndarray) -> TimeDelay:
    """tau_f(E) = Re[sum t c_f e^{iEt} / sum c_f e^{iEt}] - M Z_f / k_f - M Z0 / k0 (a.u.): the derivative
    d arg S_f / dE of the S phase referred to Z_f, the line where channel f is analysed (lines_bohr), i.e. the
    delay relative to a reflection at Z = 0. s_matrix.tsv refers each phase to that same line, so its phase
    derivative is tau + M Z_f / k_f + M Z0 / k0."""
    energies = analysis_energies(info)
    n_times, n_v, n_states = bk.coefficients.shape
    coefficients = bk.coefficients.reshape(n_times, n_v * n_states)
    dt_au = analysis_time_step_au(info)
    sums = time_sums(bk.times_au, coefficients, energies.total_hartree, dt_au)
    time_weighted_sums = time_sums(bk.times_au, coefficients * bk.times_au[:, None], energies.total_hartree, dt_au)
    mass_me, z0_bohr = info.data["molecule"]["total_mass"], info.data["wave_packet"]["z0"]
    kf_sq = 2.0 * mass_me * (energies.total_hartree[:, None] - bk.channel_energy_hartree.reshape(-1)[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        kf_per_bohr = np.sqrt(np.where(kf_sq > 0, kf_sq, np.nan))  # NaN: closed channels
        arrival_time_au = np.real(time_weighted_sums / sums)  # the mean time the flux crosses Z_f
        outgoing_flight_au = mass_me * lines_bohr[None, :] / kf_per_bohr  # from Z = 0 out to Z_f
        incoming_flight_au = mass_me * z0_bohr / energies.k0_per_bohr[:, None]  # from Z0 in to Z = 0
    probability = np.abs(sums * s_matrix_factors(info, energies, bk)) ** 2
    return TimeDelay(energies.collision_ev, arrival_time_au - outgoing_flight_au - incoming_flight_au, probability)


def masked_delays(delays: TimeDelay) -> np.ndarray:
    """tau_f(E) with NaN where P_f is below DELAY_FLOOR of the energy's P_scattered (sum over the channels):
    the phase of a small S is rounding noise, or jumps near a zero of S."""
    p_scattered = np.sum(delays.probability, axis=1, keepdims=True)
    has_weight = (delays.probability > 0) & (delays.probability >= DELAY_FLOOR * p_scattered)
    return np.where(has_weight, delays.tau_au, np.nan)
