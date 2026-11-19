"""The figures (layer 5): one figure_<name>(runs, info, options) per entry of PLOTS, composing panels. runs: the
run's segment directories in time order; info: the run's RunInfo (readers.segments.joined_info); options:
default_options(), where None means "the figure chooses". Each returns a list of (file stem, Figure) pairs.
PLOTS is the only list of entries: `list`, `plot`, `report` and the tests all go through it."""

import argparse
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from compute.analysis import analysis_energies
from compute.diffraction import diffraction, out_of_plane_fraction, p_per_order, scattering_angles
from compute.gas_phase import levels
from compute.grids import incident_parallel_momentum
from compute.results import matched_differences, results
from compute.states import level_map, mj_fractions, state_to_state
from compute.time_series import (
    analysis_lines_bohr,
    channel_densities_per_bohr,
    channel_label,
    masked_delays,
    scattering_in_time,
    strongest_channels,
    time_delays,
)
from panels.diagnostics import VERDICT_LINE_IN, draw_verdict, wrapped_verdict
from panels.diffraction import draw_orders, draw_out_of_plane, draw_scattering_angles, draw_spot_map
from panels.gas_phase import draw_potential_and_levels, draw_wave_functions
from panels.results import (
    draw_fine_comparison,
    draw_fine_difference,
    draw_mean_energies,
    draw_reaction,
    draw_scattering,
    draw_vibrational_probabilities,
)
from panels.states import draw_alignment, draw_level_map, draw_level_probabilities, draw_mj_distribution
from panels.time_series import (
    draw_channel_densities,
    draw_channel_probabilities,
    draw_convergence,
    draw_mean_positions,
    draw_norm,
    draw_time_delays,
)
from readers.binary import read_bk_coefficients, read_checkpoint
from readers.errors import ReadError
from readers.fine_analysis import read_fine_grid
from readers.run_info import RunInfo
from readers.segments import read_joined_averages
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_s_matrix, read_scattering_summary
from style import add_info_box, new_figure
from verdict import run_verdict

Figures = list[tuple[str, Figure]]
STRONGEST = 6  # channels or levels drawn against energy or time unless options.channels says otherwise
WAVE_FUNCTIONS = 3  # levels drawn by the gas_phase figure


def default_options() -> argparse.Namespace:
    """The options the figures read; None means "the figure chooses": channels (how many channels or levels),
    energy_ev (the analysis energy of the one-energy panels), azimuth_deg (the in-plane direction at normal
    incidence)."""
    return argparse.Namespace(channels=None, energy_ev=None, azimuth_deg=None)


def chosen_energy(info: RunInfo, options: argparse.Namespace) -> int:
    """The index of the analysis energy the one-energy panels show: the one nearest options.energy_ev, else the
    middle one."""
    energies_ev = analysis_energies(info).collision_ev
    if options.energy_ev is None:
        return len(energies_ev) // 2
    return int(np.argmin(np.abs(energies_ev - options.energy_ev)))


def figure_verdict(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The verdict as a figure (the report prints it as a table)."""
    lines = wrapped_verdict(run_verdict(runs, info))
    fig = new_figure(height_in=0.5 + VERDICT_LINE_IN * len(lines))
    draw_verdict(fig.add_subplot(), lines)
    add_info_box(fig, info, ["the run's files (inputs/README.md sections)", "run_info.json"])
    return [("verdict", fig)]


def figure_results(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The scattering and reaction probabilities, P(v') (two or more levels) and the mean final energies
    against the collision energy."""
    run_results = results(read_scattering_summary(runs[-1] / "scattering_summary.tsv"))
    bottom_right = "vibration" if run_results.p_vibration else "energies"
    notes = () if run_results.p_vibration else ("P(v') not drawn: fewer than two levels are scattered into",)
    fig = new_figure(n_notes=len(notes), height_in=7.0)
    axes = fig.subplot_mosaic([["scattered", "reaction"], ["energies", bottom_right]])
    draw_scattering(axes["scattered"], run_results)
    draw_reaction(axes["reaction"], run_results)
    draw_mean_energies(axes["energies"], run_results, info.data["wave_packet"]["e_parallel_ev"])
    if run_results.p_vibration:
        draw_vibrational_probabilities(axes["vibration"], run_results)
    add_info_box(fig, info, ["scattering_summary.tsv", "run_info.json"], notes=notes)
    return [("results", fig)]


def figure_averages(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The norms and mean positions against time, the segments joined."""
    averages, starts_au = read_joined_averages(runs)
    fig = new_figure()
    left, right = fig.subplots(1, 2)
    draw_norm(left, averages.times_au, averages.norm, averages.rest_norm, starts_au)
    draw_mean_positions(right, averages.times_au, averages.z_mean_bohr, averages.r_mean_bohr, starts_au)
    add_info_box(fig, info, [f"averages.dat ({len(runs)} segment{'s' if len(runs) > 1 else ''})", "run_info.json"])
    return [("averages", fig)]


def figure_gas_phase(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The asymptotic potential with the levels of the initial j, and the wave functions."""
    states = read_gas_phase_states(runs[-1] / "gas_phase_states.dat")
    potential = read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat")
    gas_levels = levels(states, potential, info.hartree_to_ev)
    j = info.data["initial_state"]["j"]
    r_opt_start_bohr = next(a["start"] for a in info.data["absorbers"]["list"] if a["name"] == "ROpt")
    fig = new_figure()
    left, right = fig.subplots(1, 2)
    draw_potential_and_levels(left, potential, gas_levels, r_opt_start_bohr, j)
    draw_wave_functions(right, states, j, count=WAVE_FUNCTIONS, r_opt_start_bohr=r_opt_start_bohr)
    add_info_box(fig, info, ["gas_phase_states.dat", "asymptotic_potential.dat", "run_info.json"])
    return [("gas_phase", fig)]


def figure_bk(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """|c(t)|^2 of the strongest channels and the convergence of P_scattered in time."""
    path = runs[-1] / "bk_coefficients.bin"  # a restart's file holds the whole run
    bk = read_bk_coefficients(path)
    if len(bk.times_au) == 0:
        raise ReadError(f"{path}: no BK record yet; the run has not finished its first analysis step")
    densities = channel_densities_per_bohr(bk)
    chosen = strongest_channels(densities, options.channels or STRONGEST)
    fig = new_figure()
    left, right = fig.subplots(1, 2)
    draw_channel_densities(left, bk.times_au, densities[:, chosen], [channel_label(bk, k) for k in chosen])
    draw_convergence(right, scattering_in_time(info, bk))
    add_info_box(fig, info, ["bk_coefficients.bin", "run_info.json"])
    return [("bk", fig)]


def figure_state_to_state(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The final rovibrational states: P(v', j') against energy and at one energy, P(mJ'), the alignment."""
    j_top = max(j for j, _ in info.data["rotational_basis"]["states"])
    energies_ev = analysis_energies(info).collision_ev
    table = read_s_matrix(runs[-1] / "s_matrix.tsv")
    states = state_to_state(table, energies_ev, info.data["grids"]["mj_max"], j_top)
    index = chosen_energy(info, options)
    fig = new_figure(height_in=7.0)
    axes = fig.subplot_mosaic([["levels", "map"], ["mj", "alignment"]])
    draw_level_probabilities(axes["levels"], states, strongest_channels(states.p_vj, options.channels or STRONGEST))
    draw_level_map(axes["map"], *level_map(states, index), energies_ev[index], j_top)
    draw_mj_distribution(axes["mj"], mj_fractions(states, index), energies_ev[index])
    draw_alignment(axes["alignment"], states)
    add_info_box(fig, info, ["s_matrix.tsv", "run_info.json"])
    return [("state_to_state", fig)]


def figure_diffraction(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """Diffraction: the spots at one energy, the orders and the out-of-plane fraction against energy."""
    energies_ev = analysis_energies(info).collision_ev
    azimuth_rad = None if options.azimuth_deg is None else np.radians(options.azimuth_deg)
    channels = diffraction(read_s_matrix(runs[-1] / "s_matrix.tsv"), info, energies_ev, azimuth_rad)
    in_plane = [c for c, inside in zip(channels.channels, channels.in_plane, strict=True) if inside]
    only_specular = in_plane == [(0, 0)]
    fig = new_figure(height_in=7.0)
    axes = fig.subplot_mosaic([["spots", "orders"], ["spots", "plane"]])
    draw_spot_map(axes["spots"], channels, chosen_energy(info, options), incident_parallel_momentum(info))
    draw_orders(axes["orders"], energies_ev, p_per_order(channels))
    draw_out_of_plane(axes["plane"], energies_ev, out_of_plane_fraction(channels), channels.plane, only_specular)
    add_info_box(fig, info, ["s_matrix.tsv", "run_info.json"])
    return [("diffraction", fig)]


def figure_scattering_angles(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The polar and azimuthal angles of the final states at one energy, and the specular direction."""
    bk = read_bk_coefficients(runs[-1] / "bk_coefficients.bin")
    index = chosen_energy(info, options)
    energies = analysis_energies(info)
    angles = scattering_angles(read_s_matrix(runs[-1] / "s_matrix.tsv"), bk, info, index)
    incident = incident_parallel_momentum(info)
    specular = (np.arctan2(incident[1], incident[0]), np.arctan2(np.hypot(*incident), energies.k0_per_bohr[index]))
    fig = new_figure(height_in=5.5)
    draw_scattering_angles(fig.add_subplot(projection="polar"), angles, energies.collision_ev[index], specular)
    add_info_box(fig, info, ["s_matrix.tsv", "bk_coefficients.bin", "run_info.json"])
    return [("scattering_angles", fig)]


def figure_time_delay(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The time delays of the strongest channels, and their probabilities."""
    bk = read_bk_coefficients(runs[-1] / "bk_coefficients.bin")
    checkpoint = read_checkpoint(runs[-1] / "checkpoint.bin")
    delays = time_delays(info, bk, analysis_lines_bohr(info, bk, checkpoint.surface_index, checkpoint.surface_fixed))
    chosen = strongest_channels(delays.probability, options.channels or STRONGEST)
    labels = [channel_label(bk, k) for k in chosen]
    notes = (
        () if checkpoint.surface_fixed else ("the specular line was never fixed: the initial state is not analysed",)
    )
    fig = new_figure(n_notes=len(notes))
    left, right = fig.subplots(1, 2)
    draw_time_delays(left, delays.collision_energy_ev, masked_delays(delays)[:, chosen], labels)
    draw_channel_probabilities(right, delays.collision_energy_ev, delays.probability[:, chosen], labels)
    add_info_box(fig, info, ["bk_coefficients.bin", "checkpoint.bin (the specular line)", "run_info.json"], notes=notes)
    return [("time_delay", fig)]


def figure_fine_analysis(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """P_scattered of the on-the-fly and the fine analysis, and their difference."""
    fine_dir = runs[-1] / "fine_analysis"
    on_the_fly = read_scattering_summary(runs[-1] / "scattering_summary.tsv")
    fine = read_scattering_summary(fine_dir / "scattering_summary.tsv")
    spectrum = read_scattering_summary(fine_dir / "spectrum" / "scattering_summary.tsv")
    grid = read_fine_grid(fine_dir / "grid.json")
    energies_ev, difference = matched_differences(fine, on_the_fly, "P_scattered")
    outside = ", ".join(f"{e_z:.3g}" for e_z, _ in grid.skipped) or "none"
    notes = (
        f"fine grid: spacing {1000 * grid.spectral_resolution_ev:.3g} meV, effective time {grid.effective_time_au:g} "
        f"a.u.; analysis energies (E_n, eV) outside it: {outside}",
    )
    fig = new_figure(n_notes=len(notes))
    left, right = fig.subplots(1, 2)
    draw_fine_comparison(left, results(on_the_fly), results(fine), results(spectrum))
    draw_fine_difference(right, energies_ev, difference)
    add_info_box(fig, info, ["scattering_summary.tsv", "fine_analysis/", "run_info.json"], notes=notes)
    return [("fine_analysis", fig)]


# name, figure function, the files it needs in the last segment's directory
FINE_ANALYSIS_FILES = ["fine_analysis/scattering_summary.tsv", "fine_analysis/spectrum/scattering_summary.tsv"]
PLOTS = [
    ("verdict", figure_verdict, []),
    ("results", figure_results, ["scattering_summary.tsv"]),
    ("state_to_state", figure_state_to_state, ["s_matrix.tsv"]),
    ("diffraction", figure_diffraction, ["s_matrix.tsv"]),
    ("scattering_angles", figure_scattering_angles, ["s_matrix.tsv", "bk_coefficients.bin"]),
    ("time_delay", figure_time_delay, ["bk_coefficients.bin", "checkpoint.bin"]),
    (
        "fine_analysis",
        figure_fine_analysis,
        ["scattering_summary.tsv", *FINE_ANALYSIS_FILES, "fine_analysis/grid.json"],
    ),
    ("bk", figure_bk, ["bk_coefficients.bin"]),
    ("averages", figure_averages, ["averages.dat"]),
    ("gas_phase", figure_gas_phase, ["gas_phase_states.dat", "asymptotic_potential.dat"]),
]


def has_files(run: Path, needed: list[str]) -> bool:
    """Whether the run directory holds every file a figure needs."""
    return all((run / name).is_file() for name in needed)
