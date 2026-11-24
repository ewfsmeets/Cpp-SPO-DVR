"""The run verdict (layer 3): reads what a run's files allow and returns the findings of
compute/diagnostics.py, worst first. A file the run has not written is a finding, never an error: the program
writes averages.dat, checkpoint.bin and the result tables only after the propagation (src/app/run.cpp), and
appends to bk_coefficients.bin every analysis step, so a run that stopped early has only the BK file."""

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from compute.absorber_reflection import returned_amplitude
from compute.analysis import analysis_energies, analysis_time_step_au, channel_energy_of_rows, energy_of_rows
from compute.densities import combined_probability, region_norms, specular_inside_line
from compute.diagnostics import (
    Finding,
    absorber_reflection,
    budget_open,
    diverged,
    grid_edge,
    lowest_converged_normal_energy_ev,
    momentum_edge,
    moving_window_start,
    norm_by_region,
    norm_inside_analysis_lines,
    norm_left,
    not_checked,
    outer_diffraction_channels,
    r_grid_ends,
    rotational_edge,
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
from compute.frame_extras import momentum_edges
from compute.frames import edge_shares, grid_edge_rises, rotational_edges, usable_frames
from compute.norm_budget import budget_check, summed_budget, usable_budget
from compute.results import Results, results
from compute.series import usable_records
from compute.time_series import p_scattered_stopped_at
from readers.binary import BkCoefficients, read_bk_coefficients, read_checkpoint
from readers.densities import DIRECTORY, read_densities
from readers.errors import ReadError
from readers.frames import DIRECTORY as FRAMES
from readers.norm_budget import DIRECTORY as BUDGET
from readers.run_info import RunInfo
from readers.segments import read_joined_averages
from readers.tables import Averages, Table, read_gas_phase_states, read_s_matrix, read_scattering_summary

PROPAGATING_RUN_MODES = (0, 1)  # new, restart (run_info.json run_mode.code)
END_OF_RUN_FILES = ("averages.dat", "scattering_summary.tsv", "s_matrix.tsv")  # written after the propagation
INITIAL_CHANNEL_ABSORBER = "ZSpOpt"  # every v of the initial channel is measured on the specular grid
OTHER_CHANNELS_ABSORBER = "ZOpt"
LOW_Z_ABSORBERS = ("ZOptLow", "ZSpOptLow")
CHANNEL_COLUMNS = ("v", "j", "mJ", "n", "m")


def run_verdict(runs: list[Path], info: RunInfo) -> list[Finding]:
    """The findings of a run (its segments in time order; `info` the joined run_info), worst first."""
    mode = info.data["run_mode"]
    if mode["code"] not in PROPAGATING_RUN_MODES:
        return [not_checked("the propagation", f"a {mode['name']} run does not propagate")]
    if "pes_cut" in info.data:  # run_info.json of a write2DPES run, which stops after the cut of the PES
        return [not_checked("the propagation", "a write2DPES run stops after the cut of the PES")]
    last = runs[-1]
    bk_path = last / "bk_coefficients.bin"
    bk = read_bk_coefficients(bk_path) if bk_path.is_file() else None
    missing = [name for name in END_OF_RUN_FILES if not (last / name).is_file()]
    if missing:
        reached = float(bk.times_au[-1]) if bk is not None and len(bk.times_au) else None
        return sorted_findings([stopped_early(missing, reached, info.t_end_au), diverged(first_bad_time(None, bk))])
    averages, _ = read_joined_averages(runs)
    run_results = results(read_scattering_summary(last / "scattering_summary.tsv"))
    rows = s_matrix_rows(last, info, bk)
    findings = [  # within a severity in this order: completeness, convergence in time, results, basis edges
        diverged(first_bad_time(averages, bk)),
        flux_inside(last, info),
        slowest(rows, info, bk),
        moving(info, bk, run_results),
        unitarity(run_results.unitarity_error, worst_unitarity_energy(run_results)),
        absorber_reflection_finding(rows, info),
        *basis_edges(rows, last, info),
        *frames_at_edges(runs),
        norm_budget_finding(runs),
        norm_left(averages.norm[-1]),
        regions_at_end(runs, info),
    ]
    return sorted_findings(findings)


def norm_budget_finding(runs: list[Path]) -> Finding | None:
    """The norm budget (norm_budget/, NormBudget) before any NaN or inf: a warning where its residual exceeds its
    error estimate by more than BUDGET_LIMIT of the initial norm; nothing without norm_budget/, a "not checked" when
    it cannot be read or holds no finite record."""
    if not (runs[-1] / BUDGET / "manifest.json").is_file():
        return None
    try:
        records, _ = usable_budget(runs)
    except ReadError as error:
        return not_checked("the norm budget", str(error))
    budget = summed_budget(records)
    check = budget_check(budget)
    k = check.window
    return budget_open(
        excess=check.excess,
        residual=float(budget.residual[k]),
        error_estimate=float(budget.error_estimate[k]),
        t_au=float(budget.times_au[k]),
    )


def regions_at_end(runs: list[Path], info: RunInfo) -> Finding | None:
    """The norm by region at the last densities/ record before any NaN or inf (writeDensities), or nothing
    without densities/; a densities/ that cannot be read or holds no finite record is a finding, not an error."""
    if not (runs[-1] / DIRECTORY / "manifest.json").is_file():
        return None
    try:
        series, grids, _ = read_densities([runs[-1]])
        series, _ = usable_records(series)
    except ReadError as error:
        return not_checked("norm by region", str(error))
    last = replace(series, times_au=series.times_au[-1:], arrays={k: v[-1:] for k, v in series.arrays.items()})
    regions = region_norms(last, info, grids, combined_probability(last))
    at_end = {label: float(values[0]) for label, values in regions.regions.items()}
    return norm_by_region(float(last.times_au[0]), at_end, float(specular_inside_line(last, info)[0]))


def frames_at_edges(runs: list[Path]) -> list[Finding | None]:
    """The basis edges in the frames (frames/, dtFrames) before any NaN or inf: information for the rotational
    edges, a warning for a grid edge that rises above its first frame's share and for a momentum grid's edge (with
    the momentum or diffraction extra) (compute.diagnostics rotational_edge, grid_edge, momentum_edge); nothing
    without frames/, a "not checked" when it cannot be read or holds no finite record."""
    if not (runs[-1] / FRAMES / "manifest.json").is_file():
        return []
    try:
        frames, _ = usable_frames(runs)
    except ReadError as error:
        return [not_checked("the basis edges of frames/", str(error))]
    edges = edge_shares(frames)
    rotational = [rotational_edge(s.label, s.values, edges.times_au) for s in rotational_edges(edges)]
    grid = [grid_edge(s.label, s.values, edges.times_au) for s in grid_edge_rises(edges)]
    momenta = [momentum_edge(s.label, s.values, edges.times_au) for s in momentum_edges(frames, edges.reference)]
    return rotational + grid + momenta


def worst_unitarity_energy(run_results: Results) -> float:
    """The collision energy (eV) of the lowest P_reaction (a NaN counts as the lowest)."""
    ranking = np.where(np.isnan(run_results.p_reaction), -np.inf, run_results.p_reaction)
    return float(run_results.collision_energy_ev[int(np.argmin(ranking))])


@dataclass(frozen=True, kw_only=True)
class SMatrixRows:
    """s_matrix.tsv, read once for the findings that use it, with per row the index of its analysis energy into
    collision_ev and its channel's outgoing normal energy (None without a BK file, whose channel energies it
    needs)."""

    table: Table
    collision_ev: np.ndarray  # the analysis energies (eV)
    energy_index: np.ndarray
    outgoing_ev: np.ndarray | None  # E_total - E_channel (eV)


def s_matrix_rows(run: Path, info: RunInfo, bk: BkCoefficients | None) -> SMatrixRows:
    """s_matrix.tsv of the run with each row's analysis energy and, with the BK file, outgoing normal energy."""
    table = read_s_matrix(run / "s_matrix.tsv")
    collision_ev = analysis_energies(info).collision_ev
    outgoing_ev = None if bk is None else outgoing_normal_energies(table, bk, info.hartree_to_ev)
    return SMatrixRows(
        table=table,
        collision_ev=collision_ev,
        energy_index=energy_of_rows(table, collision_ev),
        outgoing_ev=outgoing_ev,
    )


def basis_edges(rows: SMatrixRows, run: Path, info: RunInfo) -> list[Finding | None]:
    """The flux at the top j' of the rotational basis and in the outermost diffraction channels
    (s_matrix.tsv), and the gas-phase states at the ends of the r grid (gas_phase_states.dat)."""
    table, energy_of_row, energies_ev = rows.table, rows.energy_index, rows.collision_ev
    j_top = max(j for j, _ in info.data["rotational_basis"]["states"])
    top = table["j"].astype(int) == j_top
    return [
        top_rotational_level(fraction_of_p_scattered(table, energy_of_row, top, len(energies_ev)), energies_ev, j_top),
        outer_diffraction_channels(
            fraction_of_p_scattered(table, energy_of_row, outermost(table), len(energies_ev)), energies_ev
        ),
        gas_phase_edges(run, info, int(np.max(table["v"]))),
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


def gas_phase_edges(run: Path, info: RunInfo, v_top: int) -> Finding | None:
    """The largest |psi_v,j| at either end of the r grid, relative to its maximum, of the levels of the
    rotational basis's j with v up to v_top, the highest analysed v' (s_matrix.tsv)."""
    states = read_gas_phase_states(run / "gas_phase_states.dat")
    j_values = {j for j, _ in info.data["rotational_basis"]["states"]}
    worst, worst_level, inner = 0.0, (0, 0), False
    for psi, v, j in zip(states.wave_functions, states.v, states.j, strict=True):
        if v <= v_top and j in j_values:
            ratio = max(abs(psi[0]), abs(psi[-1])) / np.max(np.abs(psi))
            if not ratio <= worst:  # a NaN is the worst
                inner = None if np.isnan(ratio) else bool(abs(psi[0]) >= abs(psi[-1]))  # NaN: neither end
                worst, worst_level = ratio, (int(v), int(j))
    return r_grid_ends(worst, worst_level, inner)


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


def slowest(rows: SMatrixRows, info: RunInfo, bk: BkCoefficients | None) -> Finding:
    """E_n,min and the probability already in the channels that leave slower than it (s_matrix.tsv with the
    BK file's channel energies)."""
    if bk is None or rows.outgoing_ev is None:
        return not_checked("slowest converged energy", "no bk_coefficients.bin or s_matrix.tsv")
    e_n_min_ev = slowest_converged_ev(info, bk)
    n_energies = len(rows.collision_ev)
    slow_p = slow_probability(rows.energy_index, rows.outgoing_ev, rows.table["P"], n_energies, e_n_min_ev)
    return slowest_converged(e_n_min_ev, slow_p, rows.collision_ev)


def absorber_reflection_finding(rows: SMatrixRows, info: RunInfo) -> Finding:
    """The error of each channel's P from the waves the Z absorbers return (s_matrix.tsv, the BK file's channel
    energies, run_info.json's absorbers): ZSpOpt for every v of the initial channel (subtracted from the
    scattering grid and measured on the specular grid), ZOpt for the others, at the channel's outgoing normal
    energy with the total mass. "not checked" without the BK file or for an absorber power other than 2."""
    if rows.outgoing_ev is None:
        return not_checked("absorber reflection", "no bk_coefficients.bin or s_matrix.tsv")
    absorbers = info.data["absorbers"]
    if absorbers["power"] != 2:
        return not_checked("absorber reflection", f"power {absorbers['power']}, not 2")
    table, outgoing_ev = rows.table, rows.outgoing_ev
    by_name = {absorber["name"]: absorber for absorber in absorbers["list"]}
    initial = initial_channel_rows(table, info)
    reflection, amplitude = np.full(len(outgoing_ev), np.inf), np.full(len(outgoing_ev), np.inf)  # closed: inf
    for name, selected in ((INITIAL_CHANNEL_ABSORBER, initial), (OTHER_CHANNELS_ABSORBER, ~initial)):
        open_rows = selected & (outgoing_ev > 0.0)
        reflection[open_rows], amplitude[open_rows] = returned_by_absorber(
            by_name[name], outgoing_ev[open_rows] / info.hartree_to_ev, info.data["molecule"]["total_mass"]
        )
    off = tuple(
        (name, by_name[name]["strength"])
        for name in (INITIAL_CHANNEL_ABSORBER, OTHER_CHANNELS_ABSORBER)
        if not by_name[name]["strength"] > 0.0
    )
    return absorber_reflection(
        reflection=reflection,
        amplitude=amplitude,
        probability=table["P"],
        collision_energy_ev=rows.collision_ev,
        channels=np.column_stack([table[c].astype(int) for c in CHANNEL_COLUMNS]),
        energy_index=rows.energy_index,
        low_z_absorbers=any(name in by_name for name in LOW_Z_ABSORBERS),
        absorbers_off=off,
    )


def initial_channel_rows(table: Table, info: RunInfo) -> np.ndarray:
    """The rows of the initial channel (j_in, mJ_in, specular diffraction (0, 0)), every v."""
    state = info.data["initial_state"]
    return (
        (table["j"].astype(int) == state["j"])
        & (table["mJ"].astype(int) == state["mJ"])
        & (table["n"].astype(int) == 0)
        & (table["m"].astype(int) == 0)
    )


def returned_by_absorber(absorber: dict, energies_hartree: np.ndarray, mass_me: float) -> tuple[np.ndarray, np.ndarray]:
    """R and the returned amplitude a = sqrt(R) + T of an absorber of run_info.json at the given kinetic energies
    (hartree). An absorber that is off (strength <= 0) reflects nothing and returns everything: R = 0, a = 1."""
    if not absorber["strength"] > 0.0:
        return np.zeros(len(energies_hartree)), np.ones(len(energies_hartree))
    return returned_amplitude(
        energies_hartree,
        length_bohr=absorber["end"] - absorber["start"],
        strength_hartree=absorber["strength"],
        mass_me=mass_me,
    )


def outgoing_normal_energies(table: Table, bk: BkCoefficients, hartree_to_ev: float) -> np.ndarray:
    """Per row of s_matrix.tsv: its channel's outgoing normal energy E_total - E_channel (eV)."""
    return table["E_total"] - channel_energy_of_rows(table, bk) * hartree_to_ev


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
