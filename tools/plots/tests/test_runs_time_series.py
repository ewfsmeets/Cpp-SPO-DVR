"""Time delays on the oblique fixture: the exact tau_f(E) from the BK sums equals the energy derivative of the
S-matrix phase (computed with the S-matrix copy, which equals s_matrix.tsv: tests/test_runs_s_matrix.py),
and the initial state is referenced to the specular line."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from compute.analysis import AnalysisEnergies, analysis_energies, analysis_time_step_au, s_matrix_factors
from compute.time_series import analysis_lines_bohr, time_delays
from copies.s_matrix import normal_momentum, time_sums
from readers.binary import BkCoefficients, read_bk_coefficients, read_checkpoint
from readers.run_info import RunInfo, read_run_info

RUN = "dailight_oblique"
WEIGHT_FLOOR = 1e-3  # energies where P_f is above this fraction of the channel's largest P_f
STEP_HARTREE = 1e-7  # the central-difference step: arg S changes by tau * 2 step ~ 3e-4 rad, far above rounding


def shifted(info: RunInfo, energies: AnalysisEnergies, shift_hartree: float) -> AnalysisEnergies:
    """The analysis energies moved by shift_hartree (total and collision energy; k0 recomputed)."""
    collision_ev = energies.collision_ev + shift_hartree * info.hartree_to_ev
    packet, mass = info.data["wave_packet"], info.data["molecule"]["total_mass"]
    k0 = normal_momentum(
        collision_ev,
        total_mass_me=mass,
        hartree_to_ev=info.hartree_to_ev,
        kx_per_bohr=packet["kx"],
        ky_per_bohr=packet["ky"],
        skew_angle_rad=info.data["lattice"]["skew_angle_rad"],
    )
    total = energies.total_hartree + shift_hartree
    return replace(
        energies, collision_ev=collision_ev, total_ev=total * info.hartree_to_ev, total_hartree=total, k0_per_bohr=k0
    )


def s_phase_derivative(info: RunInfo, bk: BkCoefficients, step_hartree: float) -> np.ndarray:
    """d arg S / dE by a central difference of the S-matrix copy (all channels referenced to Zana, as in the
    Fortran's fine analysis and its on-the-fly analysis up to r1719), (n_energies, n_channels) in a.u."""
    energies = analysis_energies(info)
    n_times = len(bk.times_au)
    coefficients = bk.coefficients.reshape(n_times, -1)
    phases = []
    for sign in (+1, -1):
        moved = shifted(info, energies, sign * step_hartree)
        sums = time_sums(bk.times_au, coefficients, moved.total_hartree, analysis_time_step_au(info))
        phases.append(np.angle(sums * s_matrix_factors(info, moved, bk)))
    return np.angle(np.exp(1j * (phases[0] - phases[1]))) / (2 * step_hartree)  # the difference wrapped to (-pi, pi]


def rounding_floor_au(n_records: int) -> float:
    """The rounding error of the extrapolated difference: arg S is good to about n_records * eps relative,
    times 1 / sqrt(WEIGHT_FLOOR) where |S| is small; divided by the 2 * STEP_HARTREE of the difference and
    times (4 + 1) / 3 for the extrapolation."""
    return n_records * np.finfo(float).eps / np.sqrt(WEIGHT_FLOOR) / (2 * STEP_HARTREE) * 5 / 3


def test_exact_time_delay_is_the_derivative_of_the_s_phase(fixture_runs: Path) -> None:
    """Where P_f carries weight, the exact tau (all channels at Zana) equals the Richardson extrapolation
    (4 D(h) - D(2h)) / 3 of the central differences with h = 1e-7 hartree, within the estimate |D(h) - D(2h)| / 3
    of D(h)'s error plus the rounding floor, at more than 1000 points; negative: the Z0 term left out of tau
    (a realistic slip, M Z0 / k0 ~ 5000 a.u.) fails at every one of them."""
    info, bk = read_run_info(fixture_runs / RUN), read_bk_coefficients(fixture_runs / RUN / "bk_coefficients.bin")
    delays = time_delays(info, bk, np.full(bk.channel_energy_hartree.size, float(info.data["analysis"]["z_ana"])))
    fine, coarse = s_phase_derivative(info, bk, STEP_HARTREE), s_phase_derivative(info, bk, 2 * STEP_HARTREE)
    extrapolated = (4 * fine - coarse) / 3
    tolerance = np.abs(fine - coarse) / 3 + rounding_floor_au(len(bk.times_au))
    weight = (delays.probability > 0) & (delays.probability >= WEIGHT_FLOOR * np.max(delays.probability, axis=0))
    assert np.sum(weight) > 1000
    assert np.max(np.abs(delays.tau_au - extrapolated)[weight] - tolerance[weight]) <= 0.0
    mass, z0 = info.data["molecule"]["total_mass"], info.data["wave_packet"]["z0"]
    without_z0 = delays.tau_au + mass * z0 / analysis_energies(info).k0_per_bohr[:, None]
    assert np.all((np.abs(without_z0 - extrapolated) > tolerance)[weight])


def test_initial_state_is_referenced_to_the_specular_line(fixture_runs: Path) -> None:
    """The plotted delays of every v of the initial state move by -M (Z_sp - Zana) / k_f against the
    Zana-referenced ones (and no other channel moves); Z_sp is where the program fixed the specular line in
    this run (10.88 bohr, its stdout); a surface never fixed leaves every channel at Zana."""
    info, bk = read_run_info(fixture_runs / RUN), read_bk_coefficients(fixture_runs / RUN / "bk_coefficients.bin")
    z_ana = float(info.data["analysis"]["z_ana"])
    checkpoint = read_checkpoint(fixture_runs / RUN / "checkpoint.bin")
    assert checkpoint.surface_fixed
    lines = analysis_lines_bohr(info, bk, checkpoint.surface_index, checkpoint.surface_fixed)
    assert np.all(analysis_lines_bohr(info, bk, checkpoint.surface_index, False) == z_ana)  # never fixed: Zana
    moved = lines != z_ana
    state = info.data["initial_state"]
    initial = [
        tuple(label[:4]) == (state["j"], state["mJ"], 0, 0) for label in bk.labels
    ] * bk.channel_energy_hartree.shape[0]
    assert list(moved) == initial and np.all(np.abs(lines[moved] - 10.88) <= 1e-12)
    assert "fixed at Z = 10.880000" in (fixture_runs / RUN / "stdout.txt").read_text(encoding="utf-8")
    at_zana = time_delays(info, bk, np.full(lines.size, z_ana)).tau_au
    program = time_delays(info, bk, lines).tau_au
    mass = info.data["molecule"]["total_mass"]
    kf_sq = 2 * mass * (analysis_energies(info).total_hartree[:, None] - bk.channel_energy_hartree.reshape(-1))
    k_f = np.sqrt(np.where(kf_sq > 0, kf_sq, np.nan))  # NaN: closed channels
    expected = np.where(moved, -mass * (lines - z_ana) / k_f, 0.0)
    open_channels = np.isfinite(program)
    assert np.max(np.abs((program - at_zana) - expected)[open_channels]) <= 1e-9
