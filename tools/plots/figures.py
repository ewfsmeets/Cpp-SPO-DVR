"""The figures (layer 5): one figure_<name>(runs, info, options) per entry of PLOTS, composing panels. runs: the
run's segment directories in time order; info: the run's RunInfo (readers.segments.joined_info); options:
default_options(), where None means "the figure chooses". Each returns a list of (file stem, Figure) pairs.
PLOTS is the only list of entries: `list`, `plot`, `report` and the tests all go through it."""

import argparse
from pathlib import Path

from matplotlib.figure import Figure

from compute.gas_phase import levels
from compute.results import results
from compute.time_series import channel_densities_per_bohr, channel_label, scattering_in_time, strongest_channels
from panels.diagnostics import VERDICT_LINE_IN, draw_verdict, wrapped_verdict
from panels.gas_phase import draw_potential_and_levels, draw_wave_functions
from panels.results import draw_mean_energies, draw_probabilities, draw_vibrational_probabilities
from panels.time_series import draw_channel_densities, draw_convergence, draw_mean_positions, draw_norm
from readers.binary import read_bk_coefficients
from readers.errors import ReadError
from readers.run_info import RunInfo
from readers.segments import read_joined_averages
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_scattering_summary
from style import add_info_box, new_figure
from verdict import run_verdict

Figures = list[tuple[str, Figure]]
STRONGEST_CHANNELS = 6  # channels drawn by the bk figure unless options.channels says otherwise
WAVE_FUNCTIONS = 3  # levels drawn by the gas_phase figure


def default_options() -> argparse.Namespace:
    """The options the figures read; None means "the figure chooses"."""
    return argparse.Namespace(channels=None)


def figure_verdict(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The verdict as a figure (the report prints it as a table)."""
    lines = wrapped_verdict(run_verdict(runs, info))
    fig = new_figure(height_in=0.5 + VERDICT_LINE_IN * len(lines))
    draw_verdict(fig.add_subplot(), lines)
    add_info_box(fig, info, ["the run's files (inputs/README.md sections)", "run_info.json"])
    return [("verdict", fig)]


def figure_results(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The probabilities and mean final energies against the collision energy; P(v') in a panel of its own
    when two or more levels are scattered into."""
    run_results = results(read_scattering_summary(runs[-1] / "scattering_summary.tsv"))
    if run_results.p_vibration:  # its values span decades: a panel (and log axis) of its own
        notes = ()
        fig = new_figure(height_in=7.0)
        axes = fig.subplot_mosaic([["probabilities", "vibration"], ["energies", "energies"]])
        draw_vibrational_probabilities(axes["vibration"], run_results)
    else:
        notes = ("P(v') not drawn: fewer than two levels are scattered into",)
        fig = new_figure(n_notes=len(notes))
        axes = dict(zip(["probabilities", "energies"], fig.subplots(1, 2), strict=True))
    draw_probabilities(axes["probabilities"], run_results)
    draw_mean_energies(axes["energies"], run_results)
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
    chosen = strongest_channels(densities, options.channels or STRONGEST_CHANNELS)
    fig = new_figure()
    left, right = fig.subplots(1, 2)
    draw_channel_densities(left, bk.times_au, densities[:, chosen], [channel_label(bk, k) for k in chosen])
    draw_convergence(right, scattering_in_time(info, bk))
    add_info_box(fig, info, ["bk_coefficients.bin", "run_info.json"])
    return [("bk", fig)]


# name, figure function, the files it needs in the last segment's directory
PLOTS = [
    ("verdict", figure_verdict, []),
    ("results", figure_results, ["scattering_summary.tsv"]),
    ("averages", figure_averages, ["averages.dat"]),
    ("gas_phase", figure_gas_phase, ["gas_phase_states.dat", "asymptotic_potential.dat"]),
    ("bk", figure_bk, ["bk_coefficients.bin"]),
]


def has_files(run: Path, needed: list[str]) -> bool:
    """Whether the run directory holds every file a figure needs."""
    return all((run / name).is_file() for name in needed)
