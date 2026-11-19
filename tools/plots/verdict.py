"""The run verdict (layer 3): reads what a run's files allow and returns the findings of
compute/diagnostics.py, worst first. A file the run has not written is a finding, never an error: the program
writes averages.dat, checkpoint.bin and the result tables only after the propagation (src/app/run.cpp), and
appends to bk_coefficients.bin every analysis step, so a run that stopped early has only the BK file."""

from pathlib import Path

import numpy as np

from compute.analysis import analysis_energies, analysis_time_step_au
from compute.diagnostics import (
    Finding,
    diverged,
    lowest_converged_normal_energy_ev,
    moving_window_start,
    norm_inside_analysis_lines,
    norm_left,
    not_checked,
    scattering_density_per_zr,
    slow_probability,
    slowest_converged,
    sorted_findings,
    still_moving,
    stopped_early,
    unanalysed_flux,
    unitarity,
)
from compute.results import Results, results
from compute.time_series import p_scattered_stopped_at
from readers.binary import BkCoefficients, read_bk_coefficients, read_checkpoint
from readers.errors import ReadError
from readers.run_info import RunInfo
from readers.segments import read_joined_averages
from readers.tables import Averages, Table, read_s_matrix, read_scattering_summary

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
    findings = [  # within a severity in this order: completeness, convergence in time, then the results
        diverged(first_bad_time(averages, bk)),
        flux_inside(last, info),
        slowest(last, info, bk),
        moving(info, bk, run_results),
        unitarity(run_results.unitarity_error),
        norm_left(averages.norm[-1]),
    ]
    return sorted_findings(findings)


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
    E_total - E_channel (eV), with the channel energy of the BK file's state list."""
    state_of = {tuple(int(q) for q in label[:4]): i for i, label in enumerate(bk.labels)}  # (j, mJ, n, m)
    energy_index = np.array([int(np.argmin(np.abs(collision_ev - e))) for e in table["E_coll"]], dtype=int)
    channel_hartree = np.empty(len(energy_index))
    for row in range(len(energy_index)):
        key = tuple(int(table[c][row]) for c in ("j", "mJ", "n", "m"))
        if key not in state_of:
            raise ReadError(f"s_matrix.tsv: state {key} is not in bk_coefficients.bin; files of different runs?")
        channel_hartree[row] = bk.channel_energy_hartree[int(table["v"][row]), state_of[key]]
    return energy_index, table["E_total"] - channel_hartree * hartree_to_ev


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
