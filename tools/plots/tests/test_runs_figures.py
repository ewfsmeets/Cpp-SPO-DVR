"""Every figure on every kind of run: the figures each kind gives, drawn, labelled, the info box clear of the
axes, saved (a figure as PNG and PDF, a movie as GIF and PNG); the content each panel must show; every PLOTS
entry has a title and an explanation."""

import shutil
from pathlib import Path

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.contour import ContourSet
from matplotlib.figure import Figure

from compute.analysis import analysis_energies
from compute.diffraction import diffraction
from compute.frame_extras import MAX_ORDERS
from compute.frames import most_at_surface, usable_frames, window_probability
from compute.grids import cell_geometry, r_absorber_first_index
from compute.potential import RAW_PES, barrier_on, gas_phase_minimum_hartree, relative_ev
from compute.series import evenly_spaced_records
from explanations import FIGURES, TITLES
from figures import (
    MOVIE_FRAMES,
    PLOTS,
    default_options,
    figure_averages,
    figure_densities,
    figure_density_movie,
    figure_elbows,
    figure_frame_movie,
    figure_potential_maps,
    has_files,
)
from panels.potential import barrier_title
from readers.arrays import read_series
from readers.densities import read_densities
from readers.errors import ReadError
from readers.pes_cut import read_pes_cut
from readers.potential import read_potential
from readers.run_info import RunInfo, lattice_of, read_run_info
from readers.segments import joined_info, ordered_segments
from readers.tables import read_asymptotic_potential, read_s_matrix
from report import save
from style import SURFACE, Movie, add_info_box, save_movie

KINDS = {  # fixture directories of each kind of run ("stopped": the stopped_run fixture)
    "oblique": ["dailight_oblique"],
    "hard wall": ["hardwall_tiny"],
    "split": ["hardwall_tiny_split/first", "hardwall_tiny_split/second"],
    "stopped": [],
    "potential only": ["potential_only"],
    "pes cut": ["pes_cut"],
    "energy resolved": ["dailight_rho_e"],
}
ALL = {name for name, _, _ in PLOTS}
POTENTIAL = {"potential_maps", "elbows", "orientation"}  # potential/ is written at set-up
SUMS = {"energy_resolved"}  # energy_resolved/: only the runs with EnergyResolvedSites
PROPAGATING = ALL - {"pes_cut"} - SUMS  # pes_cut_2d.dat: only a write2DPES run, which stops after it
EXPECTED = {  # the PLOTS entries each kind of run allows
    "oblique": PROPAGATING,
    "hard wall": PROPAGATING,
    "split": PROPAGATING,
    "stopped": {"verdict", "gas_phase", "bk", "densities", "density_movie", *POTENTIAL},
    "potential only": {"verdict", "gas_phase", *POTENTIAL},
    "pes cut": {"verdict", "pes_cut"},
    # the energy-resolved sums only: no fine analysis, densities, potential, frames or norm budget
    "energy resolved": {
        "verdict",
        "results",
        "state_to_state",
        "diffraction",
        "scattering_angles",
        "time_delay",
        "bk",
        "averages",
        "gas_phase",
        "energy_resolved",
    },
}
PIXEL_SLACK = 1  # bounding boxes are rounded to whole pixels
XY_MAPS_COLORBARS = 2  # the oblique run: both maps hold more than the floor


def run_of(kind: str, fixture_runs: Path, stopped_run: Path) -> tuple[list[Path], RunInfo]:
    """The segment directories (in time order) and the joined run_info of a kind of run."""
    directories = [fixture_runs / r for r in KINDS[kind]] or [stopped_run]
    infos, _ = ordered_segments(directories)
    return [i.directory for i in infos], joined_info(infos)


def check_layout(fig: Figure) -> None:
    """Every visible axes has a title and both axis labels; the info box (the figure's only figure-level text)
    lies inside the figure and below every axes and figure legend (tick and axis labels included)."""
    renderer = FigureCanvasAgg(fig).get_renderer()
    fig.draw_without_rendering()
    for ax in fig.axes:
        if ax.axison and ax.get_label() != "<colorbar>":
            assert ax.get_title() and ax.get_xlabel() and ax.get_ylabel(), f"unlabelled axes {ax.get_title()!r}"
    (info_box,) = fig.texts
    info = info_box.get_window_extent(renderer)
    boxes = [ax.get_tightbbox(renderer) for ax in fig.axes if ax.get_visible()]
    boxes += [legend.get_window_extent(renderer) for legend in fig.legends]
    assert info.y1 <= min(box.y0 for box in boxes) + PIXEL_SLACK, "the info box overlaps an axes or a legend"
    assert info.x0 >= 0 and info.y0 >= 0 and info.x1 <= fig.bbox.width + PIXEL_SLACK, "the info box leaves the figure"


@pytest.mark.parametrize("kind", sorted(KINDS))
def test_every_figure_each_kind_of_run_allows(fixture_runs: Path, stopped_run: Path, tmp_path: Path, kind: str) -> None:
    """The entries whose files the run has are those EXPECTED; each draws at least one figure, laid out and
    saved (report.save: PNG and PDF, a movie GIF and PNG)."""
    runs, info = run_of(kind, fixture_runs, stopped_run)
    allowed = {name for name, _, needed in PLOTS if has_files(runs[-1], needed)}
    assert allowed == EXPECTED[kind]
    for name, figure, _ in PLOTS:
        if name in allowed:
            figures = figure(runs, info, default_options())
            assert figures, name
            for stem, fig in figures:
                check_layout(fig.fig if isinstance(fig, Movie) else fig)
                if isinstance(fig, Movie):  # two frames: the frame loop is test_a_movie_draws_its_frames_into_a_gif's
                    fig = Movie(fig=fig.fig, update=fig.update, n_frames=min(fig.n_frames, 2))
                assert all(p.stat().st_size > 0 for p in save(fig, tmp_path / stem))


def test_an_info_box_with_more_lines_than_its_strip_is_flagged(fixture_runs: Path, stopped_run: Path) -> None:
    """Negative test of check_layout: notes the figure did not make room for (new_figure(n_notes=...))
    push the info box into the axes."""
    runs, info = run_of("hard wall", fixture_runs, stopped_run)
    ((_, fig),) = figure_averages(runs, info, default_options())
    fig.texts[0].remove()
    add_info_box(fig, info, ["averages.dat"], notes=("one", "two", "three"))
    with pytest.raises(AssertionError, match="info box overlaps"):
        check_layout(fig)


def drawn(name: str, kind: str, fixture_runs: Path, stopped_run: Path) -> Figure:
    """The one figure of a PLOTS entry on a kind of run."""
    runs, info = run_of(kind, fixture_runs, stopped_run)
    figure = next(f for p, f, _ in PLOTS if p == name)
    ((_, fig),) = figure(runs, info, default_options())
    return fig


def test_split_averages_mark_the_restart(fixture_runs: Path, stopped_run: Path) -> None:
    """The averages of a two-segment run mark where the second segment starts, on both panels."""
    fig = drawn("averages", "split", fixture_runs, stopped_run)
    assert sum(t.get_text().strip() == "restart" for ax in fig.axes for t in ax.texts) == 2


def test_bk_channels_differ_by_style_and_marker(fixture_runs: Path, stopped_run: Path) -> None:
    """The six strongest channels, each named in the legend with its own line style and marker pair; the
    convergence map's colour bar is open at the top (ratios above its end)."""
    fig = drawn("bk", "oblique", fixture_runs, stopped_run)
    lines = fig.axes[0].get_lines()
    assert len(fig.axes[0].get_legend().get_texts()) == len(lines) == 6
    assert len({(line.get_linestyle(), line.get_marker()) for line in lines}) == 6
    assert fig.axes[2].get_label() == "<colorbar>" and fig.axes[2].collections


def test_results_panels(fixture_runs: Path, stopped_run: Path) -> None:
    """Two levels scattered into (oblique): P_reaction on a log axis (filled points above the unitarity error and
    open points inside it, both in this run) and on a linear one, the packet's spectrum on a linear axis from 0, a
    log P(v') panel; one level (hard wall): no P(v') panel and a note."""
    oblique = drawn("results", "oblique", fixture_runs, stopped_run)
    titles = {ax.get_title(): ax for ax in oblique.axes}
    reaction_log = titles["Reaction probability (log)\n(includes flux not yet analysed and basis loss)"]
    assert reaction_log.get_yscale() == "log"
    points = {line.get_label(): line for line in reaction_log.get_lines() if line.get_linestyle() == "None"}
    assert set(points) == {"above the unitarity error", "inside the unitarity error"}
    assert legend_texts(reaction_log) == ["above the unitarity error", "inside the unitarity error"]
    assert points["inside the unitarity error"].get_markerfacecolor() == "none"
    assert points["above the unitarity error"].get_markerfacecolor() != "none"
    assert titles["Reaction probability (linear)\n(includes flux not yet analysed and basis loss)"]
    spectrum = titles["The packet's energy spectrum\n(the S matrix divides by it)"]
    assert spectrum.get_yscale() == "linear" and spectrum.get_ylim()[0] == 0.0
    assert titles["Final vibrational state"].get_yscale() == "log"
    hard_wall = drawn("results", "hard wall", fixture_runs, stopped_run)
    assert "Final vibrational state" not in {ax.get_title() for ax in hard_wall.axes}
    assert "P(v') not drawn" in hard_wall.texts[0].get_text()


def test_every_entry_has_a_title_and_an_explanation() -> None:
    """Every PLOTS entry has a title and a text in explanations, and every text an entry."""
    assert sorted(FIGURES) == sorted(TITLES) == sorted(ALL)


def axes_titled(fig: Figure, start: str) -> list:
    """The axes of a figure whose title starts with `start`."""
    return [ax for ax in fig.axes if ax.get_title().startswith(start)]


def test_state_to_state_panels(fixture_runs: Path, stopped_run: Path) -> None:
    """Four panels and the map's colour bar; the alignment of j' = 2 is drawn (MJMax 2 holds all its mJ');
    each j''s mJ' bars sum to 1; on the hard wall (no coupling) no alignment is drawn: P(j' = 2) is noise."""
    fig = drawn("state_to_state", "oblique", fixture_runs, stopped_run)
    (alignment,) = axes_titled(fig, "Final-state alignment")
    assert len(fig.axes) == 5 and [t.get_text() for t in alignment.get_legend().get_texts()] == ["j' = 2"]
    (bars,) = axes_titled(fig, "mJ' distribution")
    assert abs(sum(patch.get_height() for patch in bars.patches) - 1.0) < 1e-12
    hard_wall = drawn("state_to_state", "hard wall", fixture_runs, stopped_run)
    (alignment,) = axes_titled(hard_wall, "Final-state alignment")
    assert alignment.get_lines() == [] or all(np.isnan(line.get_ydata()).all() for line in alignment.get_lines())


def test_diffraction_panels(fixture_runs: Path, stopped_run: Path) -> None:
    """The spot map holds the analysed channels at their K and the probabilities (two scatter layers) and
    names the plane of incidence; at 17 degrees only (0,0) is in it."""
    fig = drawn("diffraction", "oblique", fixture_runs, stopped_run)
    (spots,) = axes_titled(fig, "Diffraction at E_coll")
    assert len(spots.collections) == 2
    assert "the plane of incidence" in [t.get_text() for t in spots.get_legend().get_texts()]
    assert axes_titled(fig, "Out of the plane of incidence? (only (0,0) is in it)")
    run = fixture_runs / "dailight_oblique"
    info = read_run_info(run)
    channels = diffraction(read_s_matrix(run / "s_matrix.tsv"), info, analysis_energies(info).collision_ev, None)
    assert np.array_equal(spots.collections[0].get_offsets(), channels.k_cart_per_bohr)


def test_scattering_angles_and_time_delays(fixture_runs: Path, stopped_run: Path) -> None:
    """Polar axes with the final states and the specular marker; six channels in both time-delay panels, and
    on the hard wall only the channels that carry flux."""
    angles = drawn("scattering_angles", "oblique", fixture_runs, stopped_run)
    assert angles.axes[0].name == "polar" and len(angles.axes[0].collections[0].get_offsets()) > 1
    assert [t.get_text() for t in angles.axes[0].get_legend().get_texts()] == ["final states (area: P)", "specular"]
    delays = drawn("time_delay", "oblique", fixture_runs, stopped_run)
    assert [len(ax.get_lines()) for ax in delays.axes] == [6, 6]
    hard_wall = drawn("time_delay", "hard wall", fixture_runs, stopped_run)
    assert len(hard_wall.axes[0].get_lines()) < 6


def test_fine_analysis_panels(fixture_runs: Path, stopped_run: Path) -> None:
    """On the fly, the fine grid and the interpolated points; the box notes the fine grid."""
    fig = drawn("fine_analysis", "oblique", fixture_runs, stopped_run)
    assert [line.get_label() for line in fig.axes[0].get_lines()] == [
        "on the fly",
        "fine, its grid",
        "fine, interpolated",
    ]
    assert "fine grid: spacing" in fig.texts[0].get_text()


def annotations(ax: object) -> list[str]:
    """The stripped texts of an axes' annotations (the labels of its marked lines)."""
    return [t.get_text().strip() for t in ax.texts]


def test_densities_panels(fixture_runs: Path, stopped_run: Path) -> None:
    """The split hard-wall run: the Z carpet marks zasymp, Zana, the scattering grid's end, the Z absorbers and
    the specular analysis line, the r carpet the r absorber, each with its own colour bar; the regions' legend
    lists r >= r_c, the four Z bands, 1 - norm and what is inside each absorber; every panel marks the restart."""
    fig = drawn("densities", "split", fixture_runs, stopped_run)
    by_title = {ax.get_title(): ax for ax in fig.axes if ax.get_title()}
    z_carpet, r_carpet = by_title["The packet in Z against time"], by_title["The packet in r against time"]
    marked = set(annotations(z_carpet))
    assert {"zasymp", "Zana", "end of the scattering grid", "specular analysis line"} <= marked
    assert {"ZOpt starts (scattering wave)", "ZSpOpt starts (specular wave)"} <= marked
    assert "ROpt starts (both waves)" in annotations(r_carpet)
    assert sum(ax.get_label() == "<colorbar>" for ax in fig.axes) == 2
    legend = next(ax.get_legend() for ax in fig.axes if not ax.axison and ax.get_legend())
    labels = [t.get_text() for t in legend.get_texts()]
    assert labels[0].startswith("r >= r_c") and labels[5] == "1 - norm (removed)" and len(labels) == 9
    assert labels[6:] == ["inside ZOpt", "inside ROpt", "inside ZSpOpt"]
    assert sum(text == "restart" for ax in by_title.values() for text in annotations(ax)) == 3


def test_density_movie_marks(fixture_runs: Path, stopped_run: Path) -> None:
    """The movie's still (its last frame) marks the lines of the map, the specular analysis line and the r
    absorber, on one colour bar, and names its time."""
    _, info = run_of("hard wall", fixture_runs, stopped_run)
    fig = drawn("density_movie", "hard wall", fixture_runs, stopped_run).fig
    ax = fig.axes[0]
    marked = set(annotations(ax))
    assert {"zasymp", "Zana", "end of the scattering grid", "specular analysis line"} <= marked
    assert {"ZOpt starts (scattering wave)", "ZSpOpt starts (specular wave)", "ROpt starts (both waves)"} <= marked
    assert sum(a.get_label() == "<colorbar>" for a in fig.axes) == 1
    assert ax.get_title().endswith(f"at t = {info.t_end_au:g} a.u.")


def test_density_figures_need_a_finite_record(no_density_record_run: Path) -> None:
    """A run stopped before its first densities record: both density figures refuse with a ReadError (the report
    lists them as not drawn), not a traceback."""
    info = read_run_info(no_density_record_run)
    for figure in (figure_densities, figure_density_movie):
        with pytest.raises(ReadError, match="no finite record"):
            figure([no_density_record_run], info, default_options())


def test_density_figures_stop_at_a_non_finite_record(fixture_runs: Path, tmp_path: Path) -> None:
    """A NaN in record 3 of the specular cells (a diverged run): the figures draw records 0..2 and say so."""
    run = tmp_path / "hardwall_tiny_nan"
    shutil.copytree(fixture_runs / "hardwall_tiny", run, ignore=shutil.ignore_patterns("checkpoint.bin"))
    cells = run / "densities" / "prob_spec_zr.bin"
    values = np.fromfile(cells, "<f8")
    values[3 * read_series(run / "densities").arrays["prob_spec_zr"][0].size] = np.nan
    values.tofile(cells)
    info = read_run_info(run)
    ((_, fig),) = figure_densities([run], info, default_options())
    times_au = read_series(run / "densities").times_au
    assert f"the records from t = {times_au[3]:g} a.u. on hold a NaN or inf" in fig.texts[0].get_text()
    assert max(fig.axes[0].collections[0].get_coordinates()[:, :, 0].ravel()) < times_au[3]  # records 0..2 only
    ((_, movie),) = figure_density_movie([run], info, default_options())
    assert movie.n_frames == 3 and "NaN or inf" in movie.fig.texts[0].get_text()


def test_a_movie_draws_its_frames_into_a_gif(fixture_runs: Path, stopped_run: Path, tmp_path: Path) -> None:
    """The density movie cut to 3 frames: update is called once per frame, the GIF and the still PNG are
    written, each frame's title names its record's time, and the colour scale and the map's place stay the
    same in every frame (fixed scale, layout fixed after the first frame)."""
    runs, info = run_of("oblique", fixture_runs, stopped_run)
    ((_, movie),) = figure_density_movie(runs, info, default_options())
    series, _, _ = read_densities(runs)
    times_au = series.times_au[evenly_spaced_records(len(series.times_au), MOVIE_FRAMES)]
    frames, positions = [], []

    def recorded(k: int) -> None:
        """The movie's update, recording what frame k draws (title, colour scale) and where the map sat after the
        previous frame was drawn."""
        if k > 0:
            positions.append(movie.fig.axes[0].get_position().bounds)
        movie.update(k)
        ax = movie.fig.axes[0]
        frames.append((ax.get_title(), ax.collections[0].norm.vmin, ax.collections[0].norm.vmax))

    gif, png = save_movie(Movie(fig=movie.fig, update=recorded, n_frames=3), tmp_path / "movie")
    assert gif.stat().st_size > 0 and png.stat().st_size > 0
    assert len(frames) == 3
    for k, (title, vmin, vmax) in enumerate(frames):
        assert f"t = {times_au[k]:g} a.u." in title
        assert (vmin, vmax) == frames[0][1:]
    positions.append(movie.fig.axes[0].get_position().bounds)  # after the last frame and the PNG
    assert len(set(positions)) == 1, "the map moved between frames"


def test_potential_figures_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's sites: the maps (the barrier of each point, V_min at the lowest barrier's (Z, r), V_min at
    the bond length) draw the unit cell's outline and mark every site and the lowest barrier's point; the elbows
    show each site's V_min with its barrier at the run's nodes above V at that barrier's node with its own; the
    orientation map marks every DVR node."""
    runs, info = run_of("oblique", fixture_runs, stopped_run)
    potential = read_potential(runs[-1])
    names = [site.name for site in potential.sites]
    maps = drawn("potential_maps", "oblique", fixture_runs, stopped_run)
    panels = [ax for ax in maps.axes if ax.get_label() != "<colorbar>"]
    assert len(panels) == 3 and panels[0].get_title().startswith("Barrier")
    assert panels[1].get_title().endswith("(lowest barrier)")
    outline = cell_geometry(potential.nx, potential.ny, lattice_of(info))
    for ax in panels:
        marked = set(annotations(ax))
        assert set(names) <= marked and any(m.startswith("lowest barrier (x") for m in marked)
        assert np.allclose(ax.get_lines()[0].get_xdata(), outline.outline_x_bohr)
    elbows = drawn("elbows", "oblique", fixture_runs, stopped_run)
    titles = [ax.get_title() for ax in elbows.axes if ax.get_title()]
    top, bottom = titles[: len(names)], titles[len(names) :]
    assert [t.split(":")[0] for t in top] == names and all("at the run's nodes" in t for t in top)
    assert [t.split(", theta")[0] for t in bottom] == names
    assert all(len(ax.get_lines()) == 2 for ax in elbows.axes[: len(names)])  # the path and the barrier's mark
    orientation = drawn("orientation", "oblique", fixture_runs, stopped_run)
    (nodes,) = orientation.axes[0].get_lines()
    assert len(nodes.get_xdata()) == len(potential.theta_rad) * len(potential.phi_rad)


def test_potential_map_follows_z_and_r(fixture_runs: Path, stopped_run: Path) -> None:
    """--z and --r move the middle map to the nearest grid point (and drop "(lowest barrier)")."""
    runs, info = run_of("oblique", fixture_runs, stopped_run)
    options = default_options()
    options.z_bohr, options.r_bohr = 3.01, 1.4  # the grids: Z = -2 + 0.08 k, so 3.04; r = 0.8 + 0.15 k, so 1.40
    ((_, fig),) = figure_potential_maps(runs, info, options)
    assert fig.axes[1].get_title() == "V_min at Z = 3.04,\nr = 1.40 bohr"


def test_potential_figures_without_sites(no_sites_run: Path) -> None:
    """A run without PlotSites (the hard wall: V only rises up to r_c): one row of elbows at the lowest pass's (x, y)
    point, titled as having no saddle, that point marked on the maps, and no orientation figure (it needs a
    site)."""
    info = read_run_info(no_sites_run)
    ((_, elbows),) = figure_elbows([no_sites_run], info, default_options())
    panels = [ax for ax in elbows.axes if ax.get_label() != "<colorbar>"]
    title = panels[0].get_title()
    assert len(panels) == 1 and title.startswith("lowest barrier (x") and "no saddle before r_c" in title
    ((_, maps),) = figure_potential_maps([no_sites_run], info, default_options())
    assert any(m.startswith("lowest barrier (x") for m in annotations(maps.axes[0]))
    assert not has_files(no_sites_run, next(needed for name, _, needed in PLOTS if name == "orientation"))


def test_pes_cut_figure_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The cut of the write2DPES fixture: one panel titled with the cut's position and orientation and what the
    minimax pass found on the PES itself over the whole Z grid (the title of that barrier, recomputed here), the
    pass drawn (its path, and the mark of a saddle) in white on the sequential scale (V stays above the gas-phase
    minimum here), one colour bar, the note that V is the PES itself, and no analysis energies (it printed nan once)."""
    runs, info = run_of("pes cut", fixture_runs, stopped_run)
    cut = read_pes_cut(runs[-1], info)
    barrier = barrier_on(cut.v_hartree, r_absorber_first_index(info, len(cut.r_bohr)), RAW_PES)
    reference = gas_phase_minimum_hartree(read_asymptotic_potential(runs[-1] / "asymptotic_potential.dat"))
    height_ev = float(relative_ev(barrier.height_hartree, reference, info.hartree_to_ev))
    where = f"x 2.42, y 0 bohr, theta {cut.theta_deg:.4g}, phi 0 deg"
    expected = barrier_title(where, barrier, cut.z_bohr, cut.r_bohr, height_ev, "at this orientation")
    fig = drawn("pes_cut", "pes cut", fixture_runs, stopped_run)
    (panel,) = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
    assert panel.get_title() == expected
    assert len(panel.get_lines()) == (2 if barrier.saddle else 1)
    assert {line.get_color() for line in panel.get_lines()} == {SURFACE}
    assert len([ax for ax in fig.axes if ax.get_label() == "<colorbar>"]) == 1
    assert "the PES itself" in fig.texts[0].get_text()
    assert "no analysis energies" in fig.texts[0].get_text() and "nan" not in fig.texts[0].get_text()


def test_frame_movies_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's frame movies: one per site, numbered in PlotSites order (site1 top, site2 bridge: a name may
    hold any character), each titled with its site; V_min's contours (a ContourSet) and the zasymp mark over the
    map; the still (the figure as saved) is the shown frame nearest the one with the most probability below zasymp,
    not the last."""
    runs, info = run_of("oblique", fixture_runs, stopped_run)
    movies = figure_frame_movie(runs, info, default_options())
    assert [stem for stem, _ in movies] == ["frame_movie_site1", "frame_movie_site2"]
    ax = movies[0][1].fig.axes[0]
    assert any(isinstance(c, ContourSet) for c in ax.collections) and "zasymp" in annotations(ax)
    frames, _ = usable_frames(runs)
    records = evenly_spaced_records(len(frames.series.times_au), MOVIE_FRAMES)
    still_frame = int(np.argmin(np.abs(records - most_at_surface(frames, info.data["grids"]["nz_potential"]))))
    assert movies[0][1].still == still_frame != len(records) - 1  # save_movie saves this frame as the PNG
    assert ax.get_title() == f"rho(Z, r) above top at t = {frames.series.times_au[records[still_frame]]:g} a.u."
    assert movies[1][1].fig.axes[0].get_title().startswith("rho(Z, r) above bridge at t = ")


def test_xy_maps_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's maps over the cell, both sites marked: the left at the frame with the most probability below
    zasymp, the right inside the r absorber (Z < zasymp, r >= r_c) at the frame with the most there, its title giving
    that time and probability (the map a density, not P_reaction)."""
    runs, info = run_of("oblique", fixture_runs, stopped_run)
    maps = drawn("xy_maps", "oblique", fixture_runs, stopped_run)
    panels = [ax for ax in maps.axes if ax.get_label() != "<colorbar>"]
    frames, _ = usable_frames(runs)
    nz_surface = info.data["grids"]["nz_potential"]
    k = most_at_surface(frames, nz_surface)
    window = window_probability(
        frames, slice(0, nz_surface), slice(r_absorber_first_index(info, len(frames.r_bohr)), None)
    )
    k_absorber = int(np.argmax(window))
    assert len(panels) == 2 and len(maps.axes) == XY_MAPS_COLORBARS + 2
    assert panels[0].get_title().startswith(f"Over the cell below zasymp at t = {frames.series.times_au[k]:g} a.u.")
    assert panels[1].get_title() == (
        f"Inside the r absorber (Z < zasymp, r >= r_c)\nat t = {frames.series.times_au[k_absorber]:g} a.u.: "
        f"P in the window {window[k_absorber]:.1e}\n(a density at this time, not P_reaction)"
    )
    assert all({"top", "bridge"} <= set(annotations(ax)) for ax in panels)


def test_rotational_states_and_edges_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's state map with every state labelled (sorted by J) and one line per J; the edges: left the
    top J alone (no |mJ| edge, since MJMax = JMax), its title giving how high it gets; right the rises of the first
    and last r and the first Z points with the limit line, the title answering the question as the verdict does (yes:
    no grid edge gains more than the limit)."""
    states = drawn("rotational_states", "oblique", fixture_runs, stopped_run)
    labels = [t.get_text() for t in states.axes[0].get_yticklabels()]
    assert labels[:2] == ["J 0, mJ 0", "J 2, mJ -2"] and len(labels) == 6  # sorted by J, then mJ; every state
    assert len(states.axes[1].get_lines()) == len({label.split(",")[0] for label in labels})  # one line per J
    edges = drawn("basis_edges", "oblique", fixture_runs, stopped_run)
    rotational, grid = edges.axes[0], edges.axes[1]
    assert rotational.get_title().startswith("Rotational edges: top J up to ")
    assert [t.get_text() for t in rotational.get_legend().get_texts()] == ["top J"]
    assert grid.get_title().startswith("Does the grid hold the wave? Yes:")
    legend = [t.get_text() for t in grid.get_legend().get_texts()]
    assert legend == ["first r points", "last r points", "first Z points", "limit 1e-06 (heuristic)"]


def test_norm_budget_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's norm budget: left, the norm on the grids, the norm that left them and each absorber with its
    grid, in norm_budget/'s order (its legend beside the panel); right, the residual against its error band, the
    title answering as the verdict does (no: the fixture loses about 4e-3 of the norm in the angular transform at the
    surface)."""
    fig = drawn("norm_budget", "oblique", fixture_runs, stopped_run)
    went, legend_ax, closes = fig.axes
    assert went.get_title() == "Where did the norm go?"
    labels = [t.get_text() for t in legend_ax.get_legend().get_texts()]
    assert labels == [
        "on the grids",
        "left the grids",
        "taken by ZOpt (scattering grid)",
        "taken by ROpt (scattering grid)",
        "taken by ZSpOpt (specular grid)",
        "taken by ROpt (specular grid)",
    ]
    assert closes.get_title().startswith("Does the budget close? No:")
    assert [t.get_text() for t in closes.get_legend().get_texts()] == ["error estimate", "residual"]


def legend_texts(ax: object) -> list[str]:
    """The texts of an axes' legend."""
    return [t.get_text() for t in ax.get_legend().get_texts()]


def test_frame_extras_content(fixture_runs: Path, stopped_run: Path) -> None:
    """The oblique run's FrameExtras figures: the axis distribution at the surface and in the last frame with one
    colour bar naming the isotropic value, the homonuclear mirror noted, one square per node and theta from 0
    to past the last node; P(n, m; Z) with MAX_ORDERS orders and the rest (4 x 4 orders), the same orders in both panels,
    zasymp and Zana (equal here) marked as one line; P(v; Z) with the fixture's two analysed levels, the other levels
    below V_asy(r_end) and those above; the moments in three frames with the reference levels of <P2>; the momentum
    panels, and the edges with the limit, the title answering as the verdict does (no: the outer diffraction orders of
    the 4 x 4 grid)."""
    steric = drawn("steric", "oblique", fixture_runs, stopped_run)
    maps = [ax for ax in steric.axes if ax.get_label() != "<colorbar>"]
    (bar,) = [ax for ax in steric.axes if ax.get_label() == "<colorbar>"]
    assert [ax.get_title().split(":")[0] for ax in maps] == ["At the surface", "Last frame"]
    assert bar.get_ylabel() == f"probability per steradian (isotropic {1 / (4 * np.pi):.3f})"
    assert all("theta > 90: the mirror (180 - theta, phi + 180)" in annotations(ax) for ax in maps)
    for ax in maps:  # theta from 0 down past the last node, which lies at 90 degrees here; one square per node
        (squares,) = ax.collections
        bottom, top = ax.get_ylim()
        assert top == 0.0 and bottom > 90.0 and len(squares.get_paths()) == len(squares.get_array())
        assert np.max([path.vertices[:, 1].max() for path in squares.get_paths()]) < bottom
    diffraction = drawn("diffraction_z", "oblique", fixture_runs, stopped_run).axes
    labels = [legend_texts(ax) for ax in diffraction]
    assert labels[0] == labels[1] and len(labels[0]) == MAX_ORDERS + 1 and labels[0][-1] == "other orders"
    assert all("zasymp = Zana" in annotations(ax) for ax in diffraction)
    levels = ["v = 0", "v = 1", "other levels below V_asy(r_end)", "levels above V_asy(r_end)"]
    assert all(legend_texts(ax) == levels for ax in drawn("vibration_z", "oblique", fixture_runs, stopped_run).axes)
    alignment, energy = drawn("rotation_moments", "oblique", fixture_runs, stopped_run).axes
    frames = ["first frame", "at the surface", "in the last frame"]
    assert all([label.split(",")[0] for label in legend_texts(ax)] == frames for ax in (alignment, energy))
    assert {"isotropic 0", "upright +1", "flat -1/2"} <= set(annotations(alignment))
    k_r, k_z, edges = drawn("momentum", "oblique", fixture_runs, stopped_run).axes
    assert len(k_r.get_lines()) == len(k_z.get_lines()) == 3 + 2  # three frames, two dashed edge lines
    title = edges.get_title()
    assert title.startswith("Does the momentum grid\nhold the wave? No:") and "(outer diffraction orders)" in title
    edge_labels = ["k_r", "k_Z of the specular wave", "outer diffraction orders", "limit 1e-06 (heuristic)"]
    assert legend_texts(edges) == edge_labels
