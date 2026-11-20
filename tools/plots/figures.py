"""The figures (layer 5): one figure_<name>(runs, info, options) per entry of PLOTS, composing panels. runs: the
run's segment directories in time order; info: the run's RunInfo (readers.segments.joined_info); options:
default_options(), where None means "the figure chooses". Each returns a list of (file stem, Figure or Movie)
pairs.
PLOTS is the only list of entries: `list`, `plot`, `report` and the tests all go through it."""

import argparse
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.analysis import analysis_energies
from compute.densities import (
    combined_probability,
    density_marks,
    r_profile_per_bohr,
    region_norms,
    specular_line_bohr,
    z_profile_per_bohr,
)
from compute.diffraction import diffraction, out_of_plane_fraction, p_per_order, scattering_angles
from compute.gas_phase import levels
from compute.grids import (
    cell_geometry,
    density_per_volume,
    incident_parallel_momentum,
    periodic,
    r_absorber,
    r_absorber_first_index,
)
from compute.potential import (
    CONTOUR_LINE_EVERY,
    CONTOUR_STEP_EV,
    ENERGY_CAP_EV,
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
from compute.results import matched_differences, results
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
from panels.gas_phase import draw_potential_and_levels, draw_wave_functions
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
from readers.densities import read_densities
from readers.errors import ReadError
from readers.fine_analysis import read_fine_grid
from readers.pes_cut import read_pes_cut
from readers.potential import PotentialData, read_potential
from readers.run_info import RunInfo, lattice_of
from readers.segments import read_joined_averages
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_s_matrix, read_scattering_summary
from style import Movie, add_info_box, density_color_scale, info_lines, new_figure
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
    return [("density_movie", Movie(fig, update, len(records)))]


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
    scale = potential_color_scale(ENERGY_CAP_EV)
    notes = (MIN_OVER_ORIENTATIONS, *not_finite_notes(potential.v_min_hartree, "V_min"))
    fig = new_figure(n_notes=len(notes), width_in=10.0)
    axes = fig.subplots(1, 3)
    z, r, r_bond = potential.z_bohr[iz], potential.r_bohr[ir], potential.r_bohr[ir_bond]
    default = " (lowest barrier)" if options.z_bohr is None and options.r_bohr is None else ""
    draw_cell_map(
        axes[0], cell, periodic(heights), scale, marks, "Barrier on V_min(Z, r) of each\n(x, y) point (hatched: none)"
    )
    draw_cell_map(axes[1], cell, periodic(v[iz, ir]), scale, marks, f"V_min at Z = {z:.2f},\nr = {r:.2f} bohr{default}")
    draw_cell_map(
        axes[2], cell, periodic(v[iz, ir_bond]), scale, marks, f"V_min at Z = {z:.2f},\nr = {r_bond:.2f} bohr (bond)"
    )
    add_potential_colorbar(fig, list(axes), scale)
    add_info_box(fig, info, ["potential/", "asymptotic_potential.dat", "run_info.json"], notes=notes)
    return [("potential_maps", fig)]


PES_ITSELF = "V of the PES itself (no switching to the gas phase, no VMin/VMax), relative to min V_asy"


def figure_pes_cut(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """V(Z, r) of the PES at the fixed x, y, theta, phi of a write2DPES run, with its minimax barrier."""
    cut = read_pes_cut(runs[-1], info.data["grids"]["nz"], info.data["grids"]["nr"])
    reference = gas_phase_minimum_hartree(read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat"))
    barrier = barrier_on(cut.v_hartree, r_absorber_first_index(info, len(cut.r_bohr)), RAW_PES)
    height_ev = float(relative_ev(barrier.height_hartree, reference, info.hartree_to_ev))
    x, y, theta, phi = (value + 0.0 for value in (cut.x_bohr, cut.y_bohr, cut.theta_deg, cut.phi_deg))  # no -0
    where = f"x {x:.4g}, y {y:.4g} bohr, theta {theta:.4g}, phi {phi:.4g} deg"
    notes = (PES_ITSELF, *not_finite_notes(cut.v_hartree, "V"))
    fig = new_figure(n_notes=len(notes), height_in=4.5)
    ax = fig.subplots()
    scale = potential_color_scale(ENERGY_CAP_EV)
    v_ev = relative_ev(cut.v_hartree, reference, info.hartree_to_ev)
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
    scale = potential_color_scale(ENERGY_CAP_EV)
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
    scale = potential_color_scale(ENERGY_CAP_EV)
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
    ("densities", figure_densities, ["densities/manifest.json"]),
    ("density_movie", figure_density_movie, ["densities/manifest.json"]),
    ("gas_phase", figure_gas_phase, ["gas_phase_states.dat", "asymptotic_potential.dat"]),
    ("potential_maps", figure_potential_maps, ["potential/manifest.json", "asymptotic_potential.dat"]),
    ("elbows", figure_elbows, ["potential/manifest.json", "asymptotic_potential.dat"]),
    ("orientation", figure_orientation, ["potential/v_sites.bin", "asymptotic_potential.dat"]),
    ("pes_cut", figure_pes_cut, ["pes_cut_2d.dat", "asymptotic_potential.dat"]),
]


def has_files(run: Path, needed: list[str]) -> bool:
    """Whether the run directory holds every file a figure needs."""
    return all((run / name).is_file() for name in needed)
