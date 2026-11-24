"""The figures (layer 5): one figure_<name>(runs, info, options) per entry of PLOTS, composing panels. runs: the
run's segment directories in time order; info: the run's RunInfo (readers.segments.joined_info); options:
default_options(), where None means "the figure chooses". Each returns a list of (file stem, Figure or Movie)
pairs.
PLOTS is the only list of entries: `list`, `plot`, `report` and the tests all go through it."""

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.analysis import analysis_energies
from compute.densities import (
    Mark,
    Marks,
    combined_probability,
    density_marks,
    r_profile_per_bohr,
    region_norms,
    specular_line_bohr,
    z_profile_per_bohr,
)
from compute.diagnostics import FRAME_EDGE_LIMIT
from compute.diffraction import diffraction, out_of_plane_fraction, p_per_order, scattering_angles
from compute.energy_resolved import (
    area_density,
    in_valid_z,
    profile_change,
    stationary_state,
    valid_region,
    z_profile,
)
from compute.frame_extras import (
    NYQUIST_FRACTION,
    bound_levels,
    densest_z,
    diffraction_profiles,
    momentum_distribution,
    momentum_edges,
    rotation_moments,
    steric_distribution_at,
    vibration_profiles,
)
from compute.frames import (
    EDGE_POINTS,
    MAX_FRAME_SITES,
    FrameContours,
    contours_at,
    edge_shares,
    grid_edge_rises,
    j_populations,
    largest_edge,
    most_at_surface,
    rotational_edges,
    state_populations_z,
    usable_frames,
    window_probability,
    worst_record,
    xy_probability,
    zr_probability_at,
)
from compute.gas_phase import levels
from compute.grids import (
    cell_geometry,
    cell_volume_bohr4,
    density_per_volume,
    incident_parallel_momentum,
    other_images,
    periodic,
    r_absorber,
    r_absorber_first_index,
)
from compute.norm_budget import Budget, budget_check, summed_budget, usable_budget
from compute.potential import (
    CONTOUR_LINE_EVERY,
    CONTOUR_STEP_EV,
    ENERGY_CAP_EV,
    FRAME_CONTOUR_STEP_EV,
    MAX_ELBOW_SITES,
    RAW_PES,
    PointOfInterest,
    barrier_on,
    barriers_on_grid,
    bond_length_index,
    count_not_finite,
    gas_phase_minimum_hartree,
    lowest_barrier,
    lowest_site,
    map_point,
    node_edges_deg,
    orientation_at,
    positions_bohr,
    relative_ev,
    saddle_heights_hartree,
    site_points,
)
from compute.results import matched_differences, packet_spectrum_per_ev, results
from compute.series import evenly_spaced_records, usable_records
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
from panels.densities import (
    add_density_colorbar,
    draw_carpet,
    draw_density_map,
    draw_region_norms,
    draw_specular_line,
)
from panels.diagnostics import VERDICT_LINE_IN, draw_verdict, wrapped_verdict
from panels.diffraction import draw_orders, draw_out_of_plane, draw_scattering_angles, draw_spot_map
from panels.energy_resolved import add_rho_e_colorbar, draw_profiles, draw_rho_e_map, mark_valid_region
from panels.frame_extras import (
    add_steric_colorbar,
    draw_moment,
    draw_momentum_distribution,
    draw_steric_distribution,
    draw_z_profiles,
)
from panels.frames import (
    add_cell_density_colorbar,
    draw_cell_density,
    draw_edge_series,
    draw_j_populations,
    draw_point_density,
    draw_state_populations,
)
from panels.gas_phase import draw_potential_and_levels, draw_wave_functions
from panels.norm_budget import draw_budget_residual, draw_where_the_norm_went
from panels.potential import (
    add_potential_colorbar,
    barrier_title,
    draw_cell_map,
    draw_elbow,
    draw_orientation,
    point_label,
    potential_color_scale,
)
from panels.results import (
    draw_fine_comparison,
    draw_fine_difference,
    draw_mean_energies,
    draw_packet_spectrum,
    draw_reaction,
    draw_reaction_log,
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
from readers.densities import read_densities
from readers.energy_resolved import EnergyResolved, read_energy_resolved
from readers.errors import ReadError
from readers.fine_analysis import read_fine_grid
from readers.frames import Frames
from readers.pes_cut import read_pes_cut
from readers.potential import DIRECTORY as POTENTIAL
from readers.potential import PotentialData, read_potential
from readers.run_info import RunInfo, lattice_of
from readers.segments import read_joined_averages
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_s_matrix, read_scattering_summary
from style import (
    DENSITY_FLOOR,
    Movie,
    add_info_box,
    axis_label,
    density_color_scale,
    info_lines,
    linear_color_scale,
    new_figure,
)
from verdict import run_verdict

Figures = list[tuple[str, Figure | Movie]]
STRONGEST = 6  # channels or levels drawn against energy or time unless options.channels says otherwise
MOVIE_FRAMES = 30  # frames of a movie at most, evenly spaced over the records
WAVE_FUNCTIONS = 3  # levels drawn by the gas_phase figure


def default_options() -> argparse.Namespace:
    """The options the figures read; None means "the figure chooses": channels (how many channels or levels),
    energy_ev (the analysis energy of the one-energy panels), azimuth_deg (the in-plane direction at normal
    incidence), z_bohr and r_bohr (the point of the potential maps)."""
    return argparse.Namespace(channels=None, energy_ev=None, azimuth_deg=None, z_bohr=None, r_bohr=None)


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
    """The reaction probability (log and linear), the packet's energy spectrum, the mean final energies and P(v')
    (two or more levels) against the collision energy."""
    run_results = results(read_scattering_summary(runs[-1] / "scattering_summary.tsv"))
    energies = run_results.collision_energy_ev
    spectrum = packet_spectrum_per_ev(info, energies)
    layout = [["reaction_log", "reaction"], ["spectrum", "energies"]]
    if run_results.p_vibration:
        layout.append(["vibration", "vibration"])
    notes = () if run_results.p_vibration else ("P(v') not drawn: fewer than two levels are scattered into",)
    fig = new_figure(n_notes=len(notes), height_in=3.4 * len(layout))
    axes = fig.subplot_mosaic(layout)
    draw_reaction_log(axes["reaction_log"], run_results)
    draw_reaction(axes["reaction"], run_results)
    draw_packet_spectrum(axes["spectrum"], energies, spectrum)
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
    r_opt_start_bohr = r_absorber(info)["start"]
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
    resolution_mev = 1000 * 2 * np.pi / grid.effective_time_au * info.hartree_to_ev  # 2 pi hbar / T_eff
    notes = (
        f"fine grid: spacing {1000 * grid.spectral_resolution_ev:.3g} meV (pi hbar / T_eff), effective time "
        f"T_eff {grid.effective_time_au:g} a.u. (energy resolution 2 pi hbar / T_eff = {resolution_mev:.3g} meV); "
        f"analysis energies (E_n, eV) outside it: {outside}",
    )
    fig = new_figure(n_notes=len(notes))
    left, right = fig.subplots(1, 2)
    draw_fine_comparison(left, results(on_the_fly), results(fine), results(spectrum))
    draw_fine_difference(right, energies_ev, difference)
    add_info_box(fig, info, ["scattering_summary.tsv", "fine_analysis/", "run_info.json"], notes=notes)
    return [("fine_analysis", fig)]


def figure_densities(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The packet's Z and r profiles against time (the specular analysis line on the Z one), the norm per region,
    1 - norm and the probability inside each absorber."""
    joined, grids, starts_au = read_densities(runs)
    series, notes = usable_records(joined)
    combined = combined_probability(series)
    z_profile, r_profile = z_profile_per_bohr(combined, grids), r_profile_per_bohr(combined, grids)
    z_scale, r_scale = density_color_scale(z_profile), density_color_scale(r_profile)
    marks = density_marks(info, grids)
    times_au = series.times_au
    fig = new_figure(n_notes=len(notes), height_in=7.0)
    axes = fig.subplot_mosaic([["z", "r"], ["regions", "legend"]])
    draw_carpet(axes["z"], times_au, grids.z_bohr, z_profile, z_scale, "Z", marks.z_lines, starts_au)
    add_density_colorbar(fig, axes["z"], z_scale, "probability per bohr of Z")
    draw_specular_line(axes["z"], times_au, specular_line_bohr(series, grids))
    draw_carpet(axes["r"], times_au, grids.r_bohr, r_profile, r_scale, "r", marks.r_lines, starts_au)
    add_density_colorbar(fig, axes["r"], r_scale, "probability per bohr of r")
    regions = region_norms(series, info, grids, combined)
    draw_region_norms(axes["regions"], axes["legend"], times_au, regions, starts_au)
    segments = f"{len(runs)} segment{'s' if len(runs) > 1 else ''}"
    add_info_box(fig, info, [f"densities/ ({segments})", "run_info.json"], notes=notes)
    return [("densities", fig)]


def figure_density_movie(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The density over (Z, r) against time, as a movie (at most MOVIE_FRAMES frames)."""
    joined, grids, _ = read_densities(runs)
    series, notes = usable_records(joined)
    records = evenly_spaced_records(len(series.times_au), MOVIE_FRAMES)
    density = density_per_volume(combined_probability(series, records), grids.dz_bohr * grids.dr_bohr)
    color_scale = density_color_scale(density)  # one scale for all frames
    marks = density_marks(info, grids)
    times_au, lines_bohr = series.times_au[records], specular_line_bohr(series, grids)[records]
    sources = ["densities/", "run_info.json"]
    fig = new_figure(n_notes=len(notes))
    ax = fig.add_subplot()
    add_density_colorbar(fig, ax, color_scale, "rho(Z, r) (bohr^-2)")
    box = add_info_box(fig, info, sources, time_au=times_au[-1], notes=notes)

    def update(k: int) -> None:
        """Redraws the map, the specular line and the info box's time for frame k."""
        ax.clear()
        draw_density_map(ax, grids, density[k], color_scale, marks, times_au[k], lines_bohr[k])
        box.set_text("\n".join(info_lines(info, sources, time_au=times_au[k], notes=notes)))

    update(len(records) - 1)  # complete as a figure (the last frame); save_movie redraws every frame
    return [("density_movie", Movie(fig=fig, update=update, n_frames=len(records)))]


MIN_OVER_ORIENTATIONS = "V_min = the lower envelope of V over the DVR nodes at fixed (x, y); V relative to min V_asy"
CONTOURS = (CONTOUR_STEP_EV, CONTOUR_LINE_EVERY)


def figure_potential_maps(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The barrier of each (x, y) point, and V_min over the unit cell at the lowest barrier's (Z, r) and the bond length."""
    potential = read_potential(runs[-1])
    asymptotic = read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat")
    reference, to_ev = gas_phase_minimum_hartree(asymptotic), info.hartree_to_ev
    barriers = barriers_on_grid(potential, r_absorber_first_index(info, len(potential.r_bohr)))
    lowest = lowest_barrier(potential, barriers)
    iz, ir = map_point(potential, lowest, options.z_bohr, options.r_bohr)
    ir_bond = bond_length_index(asymptotic, potential.r_bohr)
    v = relative_ev(potential.v_min_hartree, reference, to_ev)
    heights = relative_ev(saddle_heights_hartree(barriers), reference, to_ev)
    lattice = lattice_of(info)
    cell = cell_geometry(potential.nx, potential.ny, lattice)
    marks = [(site.name, x, y) for site, x, y in positions_bohr(potential.sites, lattice)]
    marks += [(point_label(p), x, y) for p, x, y in positions_bohr([lowest], lattice)]
    images = other_images([(s.fraction_x, s.fraction_y) for s in potential.sites], info)  # the sites' only
    shown = (heights, v[iz, ir], v[iz, ir_bond])
    scale = potential_color_scale(ENERGY_CAP_EV, list(shown))
    notes = (MIN_OVER_ORIENTATIONS, *not_finite_notes(potential.v_min_hartree, "V_min"))
    fig = new_figure(n_notes=len(notes), width_in=10.0)
    axes = fig.subplots(1, 3)
    z, r, r_bond = potential.z_bohr[iz], potential.r_bohr[ir], potential.r_bohr[ir_bond]
    default = " (lowest barrier)" if options.z_bohr is None and options.r_bohr is None else ""
    titles = (
        "Barrier on V_min(Z, r) of each\n(x, y) point (hatched: none)",
        f"V_min at Z = {z:.2f},\nr = {r:.2f} bohr{default}",
        f"V_min at Z = {z:.2f},\nr = {r_bond:.2f} bohr (bond)",
    )
    for ax, values, title in zip(axes, shown, titles, strict=True):
        draw_cell_map(ax, cell, periodic(values), scale, marks, title, images)
    add_potential_colorbar(fig, list(axes), scale)
    add_info_box(fig, info, ["potential/", "asymptotic_potential.dat", "run_info.json"], notes=notes)
    return [("potential_maps", fig)]


PES_ITSELF = "V of the PES itself (no switching to the gas phase, no VMin/VMax), relative to min V_asy"


def figure_pes_cut(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """V(Z, r) of the PES at the fixed x, y, theta, phi of a write2DPES run, with its minimax barrier."""
    cut = read_pes_cut(runs[-1], info)
    reference = gas_phase_minimum_hartree(read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat"))
    barrier = barrier_on(cut.v_hartree, r_absorber_first_index(info, len(cut.r_bohr)), RAW_PES)
    height_ev = float(relative_ev(barrier.height_hartree, reference, info.hartree_to_ev))
    x, y, theta, phi = (value + 0.0 for value in (cut.x_bohr, cut.y_bohr, cut.theta_deg, cut.phi_deg))  # no -0
    where = f"x {x:.4g}, y {y:.4g} bohr, theta {theta:.4g}, phi {phi:.4g} deg"
    notes = (PES_ITSELF, *not_finite_notes(cut.v_hartree, "V"))
    fig = new_figure(n_notes=len(notes), height_in=4.5)
    ax = fig.subplots()
    v_ev = relative_ev(cut.v_hartree, reference, info.hartree_to_ev)
    scale = potential_color_scale(ENERGY_CAP_EV, [v_ev])
    title = barrier_title(where, barrier, cut.z_bohr, cut.r_bohr, height_ev, "at this orientation")
    draw_elbow(ax, cut.z_bohr, cut.r_bohr, v_ev, scale, CONTOURS, barrier, title)
    add_potential_colorbar(fig, [ax], scale)
    add_info_box(fig, info, ["pes_cut_2d.dat", "asymptotic_potential.dat", "run_info.json"], notes=notes)
    return [("pes_cut", fig)]


def not_finite_notes(values: np.ndarray, quantity: str) -> tuple[str, ...]:
    """A note when a potential (named `quantity`) holds NaN or inf (drawn hatched, like the walls above the cap)."""
    n = count_not_finite(values)
    return (f"{n} values of {quantity} are NaN or inf: hatched like the walls",) if n else ()


def figure_elbows(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """V_min(Z, r) at the sites (or the lowest barrier's point) with the minimax path, and V at the barrier's node."""
    potential = read_potential(runs[-1])
    reference = gas_phase_minimum_hartree(read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat"))
    ir_c = r_absorber_first_index(info, len(potential.r_bohr))
    sites = site_points(potential, ir_c)
    points = sites[:MAX_ELBOW_SITES] if sites else [lowest_barrier(potential, barriers_on_grid(potential, ir_c))]
    more = (f"elbows of the first {len(points)} of {len(sites)} sites",) if len(sites) > len(points) else ()
    notes = (MIN_OVER_ORIENTATIONS, *more, *not_finite_notes(potential.v_min_hartree, "V_min"))
    fig = new_figure(n_notes=len(notes), height_in=6.5 if sites else 4.0)
    axes = fig.subplots(2 if sites else 1, len(points), squeeze=False)
    shown = [potential.v_min_hartree[:, :, p.ix, p.iy] for p in points]  # the top row; the bottom row's nodes:
    shown += [potential.v_sites_hartree[p.isite] for p in points if p.isite is not None]  # all of the site's nodes
    scale = potential_color_scale(ENERGY_CAP_EV, [relative_ev(v, reference, info.hartree_to_ev) for v in shown])
    for k, point in enumerate(points):
        bottom = axes[1, k] if sites else None
        elbow_column(axes[0, k], bottom, potential, point, (reference, info.hartree_to_ev), ir_c, scale)
    add_potential_colorbar(fig, list(axes.ravel()), scale)
    add_info_box(fig, info, ["potential/", "asymptotic_potential.dat", "run_info.json"], notes=notes)
    return [("elbows", fig)]


def elbow_column(
    top: Axes,
    bottom: Axes | None,
    potential: PotentialData,
    point: PointOfInterest,
    energy: tuple[float, float],
    ir_c: int,
    scale: Normalize,
) -> None:
    """One point's elbows: V_min(Z, r) with its minimax pass (top) and, for a site, V at that pass's highest point's
    node with its own pass (bottom). energy: the reference (hartree) and the hartree-to-eV factor."""
    reference_hartree, to_ev = energy
    z, r, label = potential.z_bohr, potential.r_bohr, point_label(point)
    barrier = point.barrier
    v_min = relative_ev(potential.v_min_hartree[:, :, point.ix, point.iy], reference_hartree, to_ev)
    title = barrier_title(label, barrier, z, r, float(relative_ev(barrier.height_hartree, reference_hartree, to_ev)))
    draw_elbow(top, z, r, v_min, scale, CONTOURS, barrier, title)
    if bottom is None:
        return
    itheta, iphi = orientation_at(potential, point, barrier.iz, barrier.ir)
    v_node_hartree = potential.v_sites_hartree[point.isite, :, :, itheta, iphi]
    node_barrier = barrier_on(v_node_hartree, ir_c, potential.limits)
    node = f"{label}, theta {np.degrees(potential.theta_rad[itheta]):.0f} phi {np.degrees(potential.phi_rad[iphi]):.0f}"
    height = float(relative_ev(node_barrier.height_hartree, reference_hartree, to_ev))
    v_node = relative_ev(v_node_hartree, reference_hartree, to_ev)
    draw_elbow(bottom, z, r, v_node, scale, CONTOURS, node_barrier, barrier_title(node, node_barrier, z, r, height))


def figure_orientation(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """V over the DVR nodes (theta, phi) at the barrier point of the site with the lowest barrier."""
    potential = read_potential(runs[-1])
    reference = gas_phase_minimum_hartree(read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat"))
    site = lowest_site(site_points(potential, r_absorber_first_index(info, len(potential.r_bohr))))
    barrier = site.barrier
    values = relative_ev(potential.v_sites_hartree[site.isite, barrier.iz, barrier.ir], reference, info.hartree_to_ev)
    edges = node_edges_deg(potential.theta_rad, potential.phi_rad, potential.homonuclear)
    nodes = (np.degrees(potential.theta_rad), np.degrees(potential.phi_rad))
    scale = potential_color_scale(ENERGY_CAP_EV, [values])
    fig = new_figure()
    ax = fig.add_subplot()
    z, r = potential.z_bohr[barrier.iz], potential.r_bohr[barrier.ir]
    draw_orientation(
        ax,
        edges,
        nodes,
        values,
        scale,
        f"{site.name} (the site with the lowest barrier), at Z = {z:.2f}, r = {r:.2f} bohr",
    )
    add_potential_colorbar(fig, [ax], scale)
    add_info_box(fig, info, ["potential/", "asymptotic_potential.dat", "run_info.json"])
    return [("orientation", fig)]


FRAME_DENSITY_NOTE = "psi_full = psi_scat + s2s chi at one (x, y) grid point, summed over the rotational states"
EDGE_SHARE_LABEL = "Probability / initial norm"
EDGE_RISE_LABEL = "Rise above the first frame / initial norm"


@dataclass(frozen=True, kw_only=True)
class FramePoint:
    """An (x, y) grid point a frame movie shows, and the movie's file stem (numbered: a site's name may hold any
    character)."""

    stem: str
    name: str
    ix: int
    iy: int


def frame_points(frames: Frames) -> list[FramePoint]:
    """The points the frame movies show: the first MAX_FRAME_SITES sites, else the cell's origin."""
    if not frames.sites:
        return [FramePoint(stem="frame_movie_origin", name="the cell's origin", ix=0, iy=0)]
    sites = frames.sites[:MAX_FRAME_SITES]
    return [
        FramePoint(stem=f"frame_movie_site{k + 1}", name=site.name, ix=site.ix, iy=site.iy)
        for k, site in enumerate(sites)
    ]


def potential_with_reference(runs: list[Path]) -> tuple[PotentialData, float] | None:
    """potential/ of the last segment and min V_asy (hartree, the energy reference of its contours), or None
    without potential/."""
    if not (runs[-1] / POTENTIAL / "manifest.json").is_file():
        return None
    asymptotic = read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat")
    return read_potential(runs[-1]), gas_phase_minimum_hartree(asymptotic)


def figure_frame_movie(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The density over (Z, r) above each of the first MAX_FRAME_SITES sites (or the cell's origin) against time, as
    movies, with V_min's contours."""
    frames, notes = usable_frames(runs)
    records = evenly_spaced_records(len(frames.series.times_au), MOVIE_FRAMES)
    k_surface = most_at_surface(frames, int(info.data["grids"]["nz_potential"]))
    still_frame = int(np.argmin(np.abs(records - k_surface)))  # the shown frame nearest k_surface
    volume_bohr4 = cell_volume_bohr4(info.data["grids"])
    source = potential_with_reference(runs)
    movies = []
    for point in frame_points(frames):
        contours = None
        if source is not None:
            potential, reference_hartree = source
            contours = contours_at(
                potential, reference_hartree, info.hartree_to_ev, point.ix, point.iy, info.data["z_asymptotic"]
            )
        movie = point_movie(
            info,
            times_au=frames.series.times_au[records],
            z_bohr=frames.z_bohr[: frames.nz_scattering],
            r_bohr=frames.r_bohr,
            density=density_per_volume(zr_probability_at(frames, point.ix, point.iy)[records], volume_bohr4),
            name=point.name,
            contours=contours,
            notes=notes,
            still_frame=still_frame,
        )
        movies.append((point.stem, movie))
    return movies


def point_movie(
    info: RunInfo,
    times_au: np.ndarray,
    z_bohr: np.ndarray,
    r_bohr: np.ndarray,
    density: np.ndarray,
    name: str,
    contours: FrameContours | None,
    notes: tuple[str, ...],
    still_frame: int,
) -> Movie:
    """A movie of rho(Z, r) (density (frames, nz, nr), bohr^-4) at the (x, y) point `name`, with V_min's contours
    when given; one colour scale for all frames; its still (saved as the PNG) is frame still_frame."""
    color_scale = density_color_scale(density)
    sources = ["frames/", "run_info.json"] + (["potential/"] if contours is not None else [])
    lines = (*notes, FRAME_DENSITY_NOTE)
    if contours is not None:
        lines += (f"V_min's contours every {FRAME_CONTOUR_STEP_EV:g} eV relative to min V_asy, to zasymp",)
    fig = new_figure(n_notes=len(lines))
    ax = fig.add_subplot()
    add_density_colorbar(fig, ax, color_scale, "rho(Z, r) at (x, y) (bohr^-4)")
    box = add_info_box(fig, info, sources, time_au=times_au[still_frame], notes=lines)

    def update(k: int) -> None:
        """Redraws the map and the info box's time for frame k."""
        ax.clear()
        title = f"rho(Z, r) above {name} at t = {times_au[k]:g} a.u."
        draw_point_density(ax, z_bohr, r_bohr, density[k], color_scale, contours, title)
        box.set_text("\n".join(info_lines(info, sources, time_au=times_au[k], notes=lines)))

    update(still_frame)
    return Movie(fig=fig, update=update, n_frames=len(times_au), still=still_frame)


def figure_xy_maps(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """Where over the unit cell the packet is: below zasymp when most of it is there, and inside the r absorber."""
    frames, notes = usable_frames(runs)
    nz_surface = int(info.data["grids"]["nz_potential"])  # the points below zasymp
    ir_c = r_absorber_first_index(info, len(frames.r_bohr))
    in_absorber = window_probability(frames, slice(0, nz_surface), slice(ir_c, None))
    k_near, k_absorber = most_at_surface(frames, nz_surface), int(np.argmax(in_absorber))
    area_bohr2 = float(info.data["grids"]["xy_volume_element"])
    surface = density_per_volume(xy_probability(frames, k_near, slice(0, nz_surface), slice(None)), area_bohr2)
    absorber = density_per_volume(
        xy_probability(frames, k_absorber, slice(0, nz_surface), slice(ir_c, None)), area_bohr2
    )
    lattice = lattice_of(info)
    cell = cell_geometry(frames.nx, frames.ny, lattice)
    marks = [(p.name, x, y) for p, x, y in positions_bohr(frames.sites, lattice)]
    images = other_images([(s.fraction_x, s.fraction_y) for s in frames.sites], info)
    times_au = frames.series.times_au
    titles = [
        f"Over the cell below zasymp at t = {times_au[k_near]:g} a.u.\n(most of the packet there)",
        absorber_title(ir_c < len(frames.r_bohr), float(in_absorber[k_absorber]), times_au[k_absorber]),
    ]
    fig = new_figure(n_notes=len(notes) + 1, width_in=9.0)
    axes = fig.subplots(1, 2)
    floor_per_bohr2 = DENSITY_FLOOR * float(np.max(surface))  # both maps: noise below this is not drawn
    for ax, values, title in zip(axes, (surface, absorber), titles, strict=True):
        scale = linear_color_scale(values)
        draw_cell_density(ax, cell, periodic(values), scale, marks, floor_per_bohr2, title, images)
        if np.max(values) > floor_per_bohr2:  # a map of noise gets no colour bar
            add_cell_density_colorbar(fig, ax, scale)
    note = "summed over Z, r and the rotational states in the window, divided by dx dy sin(gamma)"
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, note))
    return [("xy_maps", fig)]


def absorber_title(has_absorber: bool, probability: float, t_au: float) -> str:
    """The title of the map inside the r absorber (Z < zasymp, r >= r_c = ROptStart): the frame with the most
    probability there and how much that is, or why there is none."""
    if not has_absorber:
        return f"No r absorber: no window r >= r_c\n(t = {t_au:g} a.u.)"
    return (
        f"Inside the r absorber (Z < zasymp, r >= r_c)\nat t = {t_au:g} a.u.: P in the window {probability:.1e}\n"
        "(a density at this time, not P_reaction)"
    )


def figure_rotational_states(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The rotational states over Z when most of the packet is below zasymp, and the J populations against time."""
    frames, notes = usable_frames(runs)
    k = most_at_surface(frames, int(info.data["grids"]["nz_potential"]))
    j_values, per_j = j_populations(frames)
    populations = state_populations_z(frames, k)
    color_scale = density_color_scale(populations)
    fig = new_figure(n_notes=len(notes), width_in=9.0)
    left, right = fig.subplots(1, 2)
    z = frames.z_bohr[: frames.nz_scattering]
    title = f"P(J, mJ; Z) at t = {frames.series.times_au[k]:g} a.u.,\nwhen most of the packet is below zasymp"
    draw_state_populations(left, z, populations, frames.states, color_scale, title)
    add_density_colorbar(fig, left, color_scale, "probability per Z point and state")
    draw_j_populations(right, frames.series.times_au, j_values, per_j, frames.starts_au)
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=notes)
    return [("rotational_states", fig)]


def figure_basis_edges(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """Does the basis hold the wave? The rotational edges (top J and |mJ|) and how much the grid edges (first and last
    r, first Z points) gain during the run, against the verdict's limit."""
    frames, notes = usable_frames(runs)
    edges = edge_shares(frames)
    rotational, rises = rotational_edges(edges), grid_edge_rises(edges)
    largest = largest_edge(rises)
    rise = float(largest.values[worst_record(largest)])
    answer = "Yes" if rise <= FRAME_EDGE_LIMIT else "No"  # the verdict's grid_edge test (NaN: No)
    fig = new_figure(n_notes=len(notes) + 2, width_in=9.0)
    left, right = fig.subplots(1, 2)
    title = f"Does the grid hold the wave? {answer}:\nlargest rise {rise:.1e} ({largest.label})"
    draw_edge_series(right, edges.times_au, rises, FRAME_EDGE_LIMIT, frames.starts_au, EDGE_RISE_LABEL, title)
    if rotational:
        top = largest_edge(rotational)
        title = f"Rotational edges: {top.label} up to {top.values[worst_record(top)]:.1e}\n(normal while it passes)"
        draw_edge_series(left, edges.times_au, rotational, None, frames.starts_au, EDGE_SHARE_LABEL, title)
    else:
        left.axis("off")  # one J and every |mJ| of it: no rotational edge
    lines = (f"relative to the initial norm ({edges.reference:.3g}); the first and last {EDGE_POINTS} grid points",)
    lines += ("grid edges: their share minus the first frame's; a rise <= 0 is not drawn on a log axis",)
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, *lines))
    return [("basis_edges", fig)]


def figure_norm_budget(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """Where the norm went (the absorbers' losses summed over the analysis windows) and whether the budget closes."""
    records, notes = usable_budget(runs)
    budget = summed_budget(records)
    fig = new_figure(n_notes=len(notes) + 1, width_in=10.0)
    axes = fig.subplot_mosaic([["went", "legend", "closes"]], width_ratios=[3, 1.2, 3])
    draw_where_the_norm_went(axes["went"], axes["legend"], budget)
    draw_budget_residual(axes["closes"], budget, budget_title(budget))
    note = f"each window {records.steps_per_window} time steps (dtAna); losses: the program's estimate per cell"
    add_info_box(fig, info, ["norm_budget/", "run_info.json"], notes=(*notes, note))
    return [("norm_budget", fig)]


def budget_title(budget: Budget) -> str:
    """The residual panel's title, answering as the verdict does (compute.norm_budget.budget_check), with the worst
    window's residual and error estimate."""
    check = budget_check(budget)
    k = check.window
    answer = "Yes" if check.closes else "No"
    return (
        f"Does the budget close? {answer}:\n"
        f"|residual| {abs(budget.residual[k]):.1e} at t = {budget.times_au[k]:g} a.u.\n"
        f"(error estimate {budget.error_estimate[k]:.1e})"
    )


def surface_and_last(frames: Frames, info: RunInfo) -> list[tuple[str, int]]:
    """The two records the extras' Z panels show: the one with the most probability below zasymp (named "at the
    surface") and the last one."""
    k_surface = most_at_surface(frames, int(info.data["grids"]["nz_potential"]))
    return [("at the surface", k_surface), ("in the last frame", len(frames.series.times_au) - 1)]


def z_marks(info: RunInfo) -> Marks:
    """The Z lines the extras' Z panels mark: zasymp and the analysis line Zana (labelled together where equal)."""
    lines = [
        Mark(label="zasymp", position_bohr=info.data["z_asymptotic"], grid=None),
        Mark(label="Zana", position_bohr=info.data["analysis"]["z_ana"], grid=None),
    ]
    return Marks(z_lines=lines, r_lines=[])


def figure_steric(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The distribution of the molecular axis per steradian at the densest Z below zasymp when most of the packet is
    there, and at the densest Z of the last frame (angular extra)."""
    frames, notes = usable_frames(runs)
    (_, k_surface), (_, k_last) = surface_and_last(frames, info)
    iz_surface = densest_z(frames, k_surface, int(info.data["grids"]["nz_potential"]))
    iz_last = densest_z(frames, k_last, frames.nz_scattering)
    at_surface = steric_distribution_at(frames, k_surface, iz_surface)
    last = steric_distribution_at(frames, k_last, iz_last)
    color_scale = linear_color_scale(np.array([at_surface.per_sr, last.per_sr]))
    fig = new_figure(n_notes=len(notes) + 1, width_in=9.0)
    left, right = fig.subplots(1, 2)
    surface_title = (
        f"At the surface: t = {at_surface.t_au:g} a.u.,\nZ = {at_surface.z_bohr:.2f} bohr (densest below zasymp)"
    )
    draw_steric_distribution(left, at_surface, color_scale, surface_title)
    last_title = f"Last frame: t = {last.t_au:g} a.u.,\nZ = {last.z_bohr:.2f} bohr (densest)"
    draw_steric_distribution(right, last, color_scale, last_title)
    add_steric_colorbar(fig, [left, right], color_scale)
    note = "per steradian at the DVR nodes, normalised to 1 over the sphere at that Z (prob_angle_z / w_t dphi)"
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, note))
    return [("steric", fig)]


def figure_diffraction_z(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The probability per parallel-momentum order against Z at the surface and in the last frame, the same orders in
    both (diffraction extra)."""
    frames, notes = usable_frames(runs)
    shown = surface_and_last(frames, info)
    profiles = diffraction_profiles(frames, [k for _, k in shown])
    z = frames.z_bohr[: frames.nz_scattering]
    fig = new_figure(n_notes=len(notes) + 1, width_in=9.0)
    for ax, (name, _), profile in zip(fig.subplots(1, 2), shown, profiles, strict=True):
        draw_z_profiles(ax, z, profile, z_marks(info), f"P(n, m; Z) {name},\nt = {profile.t_au:g} a.u.")
    note = "parallel momentum K = k_par + n b1 + m b2 of psi_full; only open orders at Z >= zasymp are diffraction"
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, note))
    return [("diffraction_z", fig)]


def analysed_levels(run: Path) -> int:
    """How many vibrational levels the analysis resolves: one more than the largest v of s_matrix.tsv (1 without it)."""
    s_matrix = run / "s_matrix.tsv"
    if not s_matrix.is_file():
        return 1
    return int(np.max(read_s_matrix(s_matrix)["v"])) + 1


def figure_vibration_z(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The probability per vibrational level (the analysed v, the other levels below V_asy at the last r point, those
    above it) against Z at the surface and in the last frame (vibration extra)."""
    frames, notes = usable_frames(runs)
    states = read_gas_phase_states(runs[-1] / "gas_phase_states.dat")
    potential = read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat")
    gas_levels = levels(states, potential, info.hartree_to_ev)
    bound = bound_levels(gas_levels, frames.extras.vibration_j_values, len(frames.r_bohr))
    n_analysed = analysed_levels(runs[-1])
    z = frames.z_bohr[: frames.nz_scattering]
    fig = new_figure(n_notes=len(notes) + 1, width_in=9.0)
    for ax, (name, k) in zip(fig.subplots(1, 2), surface_and_last(frames, info), strict=True):
        profiles = vibration_profiles(frames, k, bound, n_analysed)
        draw_z_profiles(ax, z, profiles, z_marks(info), f"P(v; Z) {name},\nt = {profiles.t_au:g} a.u.")
    note = "projection on the gas-phase levels phi_vJ(r), summed over J: also where the molecule is not free"
    sources = ["frames/", "gas_phase_states.dat", "asymptotic_potential.dat", "s_matrix.tsv", "run_info.json"]
    add_info_box(fig, info, sources, notes=(*notes, note))
    return [("vibration_z", fig)]


def figure_rotation_moments(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """<P2(cos theta)> and the rotational kinetic energy <J(J+1) / (2 mu r^2)> against Z in the first frame, at the
    surface and in the last frame (rotation_moments extra)."""
    frames, notes = usable_frames(runs)
    records = [("first frame", 0), *surface_and_last(frames, info)]
    moments = [(name, rotation_moments(frames, k, info.hartree_to_ev)) for name, k in records]
    z = frames.z_bohr[: frames.nz_scattering]
    fig = new_figure(n_notes=len(notes) + 1, width_in=9.0)
    left, right = fig.subplots(1, 2)
    p2 = [(f"{name}, t = {m.t_au:g} a.u.", m.p2, m.probability) for name, m in moments]
    references = [("isotropic 0", 0.0), ("upright +1", 1.0), ("flat -1/2", -0.5)]
    draw_moment(left, z, p2, references, "<P2(cos theta)>", "Alignment of the axis against Z")
    energy = [(f"{name}, t = {m.t_au:g} a.u.", m.rotational_energy_ev, m.probability) for name, m in moments]
    energy_label = axis_label("<J(J+1) / (2 mu r^2)>", "eV")
    draw_moment(right, z, energy, [], energy_label, "Rotational kinetic energy against Z")
    note = (
        f"moments of psi_full at each Z; not drawn where the probability there is below {DENSITY_FLOOR:g} of the most"
    )
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, note))
    return [("rotation_moments", fig)]


def figure_momentum(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """P(k_r) of psi_full and P(k_Z) of the specular wave in the first frame, at the surface and in the last frame,
    and whether the momentum grids hold the wave: what the run brings to their edges against time (momentum extra;
    the outer diffraction orders with the diffraction extra)."""
    frames, notes = usable_frames(runs)
    records = [("first frame", 0), *surface_and_last(frames, info)]
    fig = new_figure(n_notes=len(notes) + 1, width_in=11.0)
    left, middle, right = fig.subplots(1, 3)
    for ax, dataset, quantity, title in (
        (left, "prob_kr", "k_r", "P(k_r) of psi_full"),
        (middle, "prob_kz_spec", "k_Z", "P(k_Z) of the specular wave"),
    ):
        distributions = [momentum_distribution(frames, k, dataset) for _, k in records]
        series = [(f"{name}, t = {d.t_au:g} a.u.", d) for (name, _), d in zip(records, distributions, strict=True)]
        draw_momentum_distribution(ax, series, quantity, title)
    edges = momentum_edges(frames, edge_shares(frames).reference)
    largest = largest_edge(edges)
    rise = float(largest.values[worst_record(largest)])
    answer = "Yes" if rise <= FRAME_EDGE_LIMIT else "No"  # the verdict's momentum_edge test (NaN: No)
    title = f"Does the momentum grid\nhold the wave? {answer}: up to\n{rise:.1e} ({largest.label})"
    draw_edge_series(right, frames.series.times_au, edges, FRAME_EDGE_LIMIT, frames.starts_au, EDGE_RISE_LABEL, title)
    note = (
        f"edges: |k| above {NYQUIST_FRACTION:g} of the grid's largest (dashed); with the diffraction extra the orders "
        f"with |n| >= {NYQUIST_FRACTION:g} nX/2 or |m| >= {NYQUIST_FRACTION:g} nY/2"
    )
    add_info_box(fig, info, ["frames/", "run_info.json"], notes=(*notes, note))
    return [("momentum", fig)]


FINE_ANALYSIS_FILES = ["fine_analysis/scattering_summary.tsv", "fine_analysis/spectrum/scattering_summary.tsv"]

MAX_SUM_ENERGIES = 3  # the energy_resolved figure shows the first this many energies (one row each)
MAX_SUM_SITES = 6  # ... and the profiles of the first this many sites (six lines stay legible: own colour and marker)
DENSITY_NOT_PROBABILITY = "A rho_E: a stationary density per unit-cell area, not a reaction probability"


def sum_figure_notes(sums: EnergyResolved, rows: int, sites: range) -> tuple[str, ...]:
    """The energy_resolved figure's notes: what A rho_E is, and which energies and sites it leaves out."""
    notes = [DENSITY_NOT_PROBABILITY]
    if len(sums.collision_ev) > rows:
        notes.append(f"the first {rows} of {len(sums.collision_ev)} energies")
    if len(sums.sites) > len(sites):
        notes.append(f"profiles of the first {len(sites)} of {len(sums.sites)} sites")
    return tuple(notes)


def figure_energy_resolved(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The stationary scattering state at the first MAX_SUM_ENERGIES energies of the sums: A rho_E(Z, r) at the
    first site and its r integral against Z at the first MAX_SUM_SITES sites, in the valid region."""
    sums = read_energy_resolved(runs[-1])
    region = valid_region(info)
    rows = min(len(sums.collision_ev), MAX_SUM_ENERGIES)
    sites = range(min(len(sums.sites), MAX_SUM_SITES))
    notes = sum_figure_notes(sums, rows, sites)
    fig = new_figure(n_notes=len(notes), width_in=10.0, height_in=3.4 * rows)
    axes = fig.subplots(rows, 2, squeeze=False)
    valid_z = in_valid_z(sums.z_bohr, region)
    inside = valid_z[:, None] & (sums.r_bohr < region.r_high_bohr)[None, :]
    for e in range(rows):
        densities = [area_density(stationary_state(sums, info, k, e), info) for k in sites]
        profiles = [
            (sums.sites[k].name, np.where(valid_z, z_profile(d, sums.r_bohr, region), np.nan))
            for k, d in zip(sites, densities, strict=True)
        ]
        shown = np.where(inside, densities[0], np.nan)
        scale = linear_color_scale(shown[inside])
        changed = max(profile_change(sums, info, k, e, region) for k in sites)
        energy = f"E_coll {sums.collision_ev[e]:g} eV"
        draw_rho_e_map(
            axes[e, 0], sums.z_bohr, sums.r_bohr, shown, scale, f"{energy}: A rho_E above {sums.sites[0].name}"
        )
        mark_valid_region(axes[e, 0], region, horizontal_r=True)
        add_rho_e_colorbar(fig, axes[e, 0], scale)
        title = f"{energy}: integral over r\n(changed after t = {sums.partial_time_au:g} a.u. by up to {changed:.1e})"
        draw_profiles(axes[e, 1], sums.z_bohr, profiles, region, title)
    add_info_box(fig, info, ["energy_resolved/", "run_info.json"], notes=notes)
    return [("energy_resolved", fig)]


# name, figure function, the files it needs in the last segment's directory
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
    ("densities", figure_densities, ["densities/manifest.json"]),
    ("density_movie", figure_density_movie, ["densities/manifest.json"]),
    ("gas_phase", figure_gas_phase, ["gas_phase_states.dat", "asymptotic_potential.dat"]),
    ("potential_maps", figure_potential_maps, ["potential/manifest.json", "asymptotic_potential.dat"]),
    ("elbows", figure_elbows, ["potential/manifest.json", "asymptotic_potential.dat"]),
    ("orientation", figure_orientation, ["potential/v_sites.bin", "asymptotic_potential.dat"]),
    ("pes_cut", figure_pes_cut, ["pes_cut_2d.dat", "asymptotic_potential.dat"]),
    ("frame_movie", figure_frame_movie, ["frames/manifest.json"]),
    ("xy_maps", figure_xy_maps, ["frames/manifest.json"]),
    ("rotational_states", figure_rotational_states, ["frames/manifest.json"]),
    ("basis_edges", figure_basis_edges, ["frames/manifest.json"]),
    ("steric", figure_steric, ["frames/prob_angle_z.bin"]),
    ("diffraction_z", figure_diffraction_z, ["frames/pop_nm_z.bin"]),
    ("vibration_z", figure_vibration_z, ["frames/pop_vj_z.bin", "gas_phase_states.dat", "asymptotic_potential.dat"]),
    ("rotation_moments", figure_rotation_moments, ["frames/p2_numerator_z.bin"]),
    ("momentum", figure_momentum, ["frames/prob_kr.bin"]),
    ("norm_budget", figure_norm_budget, ["norm_budget/manifest.json"]),
    ("energy_resolved", figure_energy_resolved, ["energy_resolved/manifest.json"]),
]


def has_files(run: Path, needed: list[str]) -> bool:
    """Whether the run directory holds every file a figure needs."""
    return all((run / name).is_file() for name in needed)
