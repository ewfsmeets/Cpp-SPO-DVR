"""The run verdict (layer 3): reads what a run's files allow and returns the findings of
compute/diagnostics.py, worst first. A file the run has not written is a finding, never an error: the program
writes averages.dat, checkpoint.bin and the result tables only after the propagation (src/app/run.cpp), and
appends to bk_coefficients.bin every analysis step, so a run that stopped early has only the BK file."""

from pathlib import Path

import numpy as np

from compute.analysis import analysis_energies, analysis_time_step_au, channel_energy_of_rows, energy_of_rows
from compute.densities import density_grids, region_norms
from compute.diagnostics import (
    Finding,
    diverged,
    lowest_converged_normal_energy_ev,
    moving_window_start,
    norm_by_region,
    norm_inside_analysis_lines,
    norm_left,
    not_checked,
    outer_diffraction_channels,
    r_grid_ends,
    scattering_density_per_zr,
    slow_probability,
    slowest_converged,
    sorted_findings,
    still_moving,
    stopped_early,
    top_rotational_level,
    unanalysed_flux,
    unitarity,
)
from compute.results import Results, results
from compute.time_series import p_scattered_stopped_at
from readers.arrays import read_series
from readers.binary import BkCoefficients, read_bk_coefficients, read_checkpoint
from readers.run_info import RunInfo
from readers.segments import read_joined_averages
from readers.tables import Averages, Table, read_gas_phase_states, read_s_matrix, read_scattering_summary

PROPAGATING_RUN_MODES = (0, 1)  # new, restart (run_info.json run_mode.code)
END_OF_RUN_FILES = ("averages.dat", "scattering_summary.tsv", "s_matrix.tsv")  # written after the propagation


def run_verdict(runs: list[Path], info: RunInfo) -> list[Finding]:
    """The findings of a run (its segments in time order; `info` the joined run_info), worst first."""
    mode = info.data["run_mode"]
    if mode["code"] not in PROPAGATING_RUN_MODES:
        return [not_checked("the propagation", f"a {mode['name']} run does not propagate")]
    last = runs[-1]
    bk_path = last / "bk_coefficients.bin"
    bk = read_bk_coefficients(bk_path) if bk_path.is_file() else None
    missing = [name for name in END_OF_RUN_FILES if not (last / name).is_file()]
    if missing:
        reached = float(bk.times_au[-1]) if bk is not None and len(bk.times_au) else None
        return sorted_findings([stopped_early(missing, reached, info.t_end_au), diverged(first_bad_time(None, bk))])
    averages, _ = read_joined_averages(runs)
    run_results = results(read_scattering_summary(last / "scattering_summary.tsv"))
    findings = [  # within a severity in this order: completeness, convergence in time, results, basis edges
        diverged(first_bad_time(averages, bk)),
        flux_inside(last, info),
        slowest(last, info, bk),
        moving(info, bk, run_results),
        unitarity(run_results.unitarity_error, worst_unitarity_energy(run_results)),
        *basis_edges(last, info),
        norm_left(averages.norm[-1]),
        regions_at_end(runs, info),
    ]
    return sorted_findings(findings)


def regions_at_end(runs: list[Path], info: RunInfo) -> Finding | None:
    """The norm by region at the last densities/ record (writeDensities), or nothing without one."""
    if not (runs[-1] / "densities" / "manifest.json").is_file():
        return None
    series = read_series(runs[-1] / "densities")
    if len(series.times_au) == 0:
        return None
    regions = region_norms(series, info, density_grids(series))
    at_end = {label: float(values[-1]) for label, values in regions.regions.items()}
    return norm_by_region(float(series.times_au[-1]), at_end)


def worst_unitarity_energy(run_results: Results) -> float:
    """The collision energy (eV) of the lowest P_reaction (a NaN counts as the lowest)."""
    ranking = np.where(np.isnan(run_results.p_reaction), -np.inf, run_results.p_reaction)
    return float(run_results.collision_energy_ev[int(np.argmin(ranking))])


def basis_edges(run: Path, info: RunInfo) -> list[Finding | None]:
    """The flux at the top j' of the rotational basis and in the outermost diffraction channels
    (s_matrix.tsv), and the gas-phase states at the ends of the r grid (gas_phase_states.dat)."""
    table = read_s_matrix(run / "s_matrix.tsv")
    energies_ev = analysis_energies(info).collision_ev
    energy_of_row = energy_of_rows(table, energies_ev)
    j_top = max(j for j, _ in info.data["rotational_basis"]["states"])
    top = table["j"].astype(int) == j_top
    return [
        top_rotational_level(fraction_of_p_scattered(table, energy_of_row, top, len(energies_ev)), energies_ev, j_top),
        outer_diffraction_channels(
            fraction_of_p_scattered(table, energy_of_row, outermost(table), len(energies_ev)), energies_ev
        ),
        gas_phase_edges(run, info),
    ]


def fraction_of_p_scattered(
    table: Table, energy_of_row: np.ndarray, selected: np.ndarray, n_energies: int
) -> np.ndarray:
    """Per analysis energy: P of the selected rows / P of all rows (NaN where nothing is scattered)."""
    part = np.bincount(energy_of_row[selected], weights=table["P"][selected], minlength=n_energies)
    total = np.bincount(energy_of_row, weights=table["P"], minlength=n_energies)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0, part / total, np.nan)


def outermost(table: Table) -> np.ndarray:
    """The rows whose diffraction channel lies at an end of the analysed n or m range, (0, 0) excepted (the
    range of a dimension without freedom is one value, which is no edge)."""
    n, m = table["n"].astype(int), table["m"].astype(int)
    on_n_edge = ((n == n.min()) | (n == n.max())) & (n.min() < n.max())
    on_m_edge = ((m == m.min()) | (m == m.max())) & (m.min() < m.max())
    return (on_n_edge | on_m_edge) & ~((n == 0) & (m == 0))


def gas_phase_edges(run: Path, info: RunInfo) -> Finding | None:
    """The largest |psi_v,j| at either end of the r grid, relative to its maximum, of the levels of the
    rotational basis's j with v up to the highest analysed v' (s_matrix.tsv)."""
    states = read_gas_phase_states(run / "gas_phase_states.dat")
    j_values = {j for j, _ in info.data["rotational_basis"]["states"]}
    v_top = int(np.max(read_s_matrix(run / "s_matrix.tsv")["v"]))
    worst, worst_level = 0.0, (0, 0)
    for psi, v, j in zip(states.wave_functions, states.v, states.j, strict=True):
        if v <= v_top and j in j_values:
            ratio = max(abs(psi[0]), abs(psi[-1])) / np.max(np.abs(psi))
            if not ratio <= worst:  # a NaN is the worst
                worst, worst_level = ratio, (int(v), int(j))
    return r_grid_ends(worst, worst_level)


def first_bad_time(averages: Averages | None, bk: BkCoefficients | None) -> float | None:
    """The earliest time at which the averages or the BK coefficients hold a NaN or inf (a.u.), or None."""
    times = []
    if averages is not None and averages.first_non_finite is not None:
        times.append(float(averages.times_au[averages.first_non_finite]))
    if bk is not None and bk.first_non_finite is not None:
        times.append(float(bk.times_au[bk.first_non_finite]))
    return min(times) if times else None


def slowest_converged_ev(info: RunInfo, bk: BkCoefficients | None) -> float:
    """E_n,min of the whole run (eV), for the slowest incident normal energy EMinWP. The packet started one
    analysis step before the first BK record (the BK file of a restart holds the whole run's records)."""
    packet = info.data["wave_packet"]
    has_records = bk is not None and len(bk.times_au) > 0
    run_start_au = float(bk.times_au[0]) - analysis_time_step_au(info) if has_records else info.t_start_au
    return lowest_converged_normal_energy_ev(
        info.t_end_au - run_start_au,
        z0_bohr=packet["z0"],
        z_analysis_bohr=info.data["analysis"]["z_ana"],
        slowest_normal_energy_ev=packet["e_min_ev"],
        total_mass_me=info.data["molecule"]["total_mass"],
        hartree_to_ev=info.hartree_to_ev,
    )


def slowest(run: Path, info: RunInfo, bk: BkCoefficients | None) -> Finding:
    """E_n,min and the probability already in the channels that leave slower than it (s_matrix.tsv with the
    BK file's channel energies)."""
    if bk is None or not (run / "s_matrix.tsv").is_file():
        return not_checked("slowest converged energy", "no bk_coefficients.bin or s_matrix.tsv")
    e_n_min_ev = slowest_converged_ev(info, bk)
    table = read_s_matrix(run / "s_matrix.tsv")
    collision_ev = analysis_energies(info).collision_ev
    energy_index, outgoing_ev = outgoing_normal_energies(table, bk, collision_ev, info.hartree_to_ev)
    slow_p = slow_probability(energy_index, outgoing_ev, table["P"], len(collision_ev), e_n_min_ev)
    return slowest_converged(e_n_min_ev, slow_p, collision_ev)


def outgoing_normal_energies(
    table: Table, bk: BkCoefficients, collision_ev: np.ndarray, hartree_to_ev: float
) -> tuple[np.ndarray, np.ndarray]:
    """Per row of s_matrix.tsv: the index of its analysis energy, and its channel's outgoing normal energy
    E_total - E_channel (eV)."""
    return energy_of_rows(table, collision_ev), table["E_total"] - channel_energy_of_rows(table, bk) * hartree_to_ev


def flux_inside(run: Path, info: RunInfo) -> Finding | None:
    """The unanalysed-flux finding from checkpoint.bin, or "not checked" without it (CheckpointFormat
    "fortran" writes fort.9 instead)."""
    path = run / "checkpoint.bin"
    if not path.is_file():
        return not_checked("unanalysed flux", "no checkpoint.bin (CheckpointFormat fortran writes fort.9 instead)")
    checkpoint = read_checkpoint(path)
    grids, iz_analysis = info.data["grids"], info.data["analysis"]["iz_ana"]
    inside = norm_inside_analysis_lines(
        scattering_density_per_zr(checkpoint.scattering[: iz_analysis * checkpoint.n_r]),  # only Z < Zana
        np.abs(checkpoint.specular) ** 2,
        iz_analysis,
        checkpoint.surface_index,
        checkpoint.n_r,
        dz_bohr=grids["dz"],
        dr_bohr=grids["dr"],
        xy_volume_element_bohr2=grids["xy_volume_element"],
    )
    return unanalysed_flux(inside)


def moving(info: RunInfo, bk: BkCoefficients | None, run_results: Results) -> Finding | None:
    """Whether P_scattered still changed over the last part of the run (BK time sums stopped earlier)."""
    if bk is None or len(bk.times_au) < 2:
        return not_checked("P_scattered still moving", "fewer than two BK records")
    earlier = p_scattered_stopped_at(info, bk, moving_window_start(bk.times_au))
    end = p_scattered_stopped_at(info, bk, len(bk.times_au) - 1)
    return still_moving(end, earlier, run_results.collision_energy_ev)
