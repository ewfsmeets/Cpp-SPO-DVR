"""Every figure on every kind of run: the figures each kind gives, drawn, labelled, the info box clear of the
axes, saved as PNG and PDF; the content each panel must show; every PLOTS entry has a title and an
explanation."""

from pathlib import Path

import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from compute.analysis import analysis_energies
from compute.diffraction import diffraction
from explanations import FIGURES, TITLES
from figures import PLOTS, default_options, figure_averages, has_files
from readers.run_info import RunInfo, read_run_info
from readers.segments import joined_info, ordered_segments
from readers.tables import read_s_matrix
from style import add_info_box, save_figure

KINDS = {  # fixture directories of each kind of run ("stopped": the stopped_run fixture)
    "oblique": ["dailight_oblique"],
    "hard wall": ["hardwall_tiny"],
    "split": ["hardwall_tiny_split/first", "hardwall_tiny_split/second"],
    "stopped": [],
    "potential only": ["potential_only"],
}
ALL = {name for name, _, _ in PLOTS}
EXPECTED = {  # the PLOTS entries each kind of run allows
    "oblique": ALL,
    "hard wall": ALL,
    "split": ALL,
    "stopped": {"verdict", "gas_phase", "bk"},
    "potential only": {"verdict", "gas_phase"},
}
PIXEL_SLACK = 1  # bounding boxes are rounded to whole pixels


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
    saved (PNG and PDF)."""
    runs, info = run_of(kind, fixture_runs, stopped_run)
    allowed = {name for name, _, needed in PLOTS if has_files(runs[-1], needed)}
    assert allowed == EXPECTED[kind]
    for name, figure, _ in PLOTS:
        if name in allowed:
            figures = figure(runs, info, default_options())
            assert figures, name
            for stem, fig in figures:
                check_layout(fig)
                assert all(p.stat().st_size > 0 for p in save_figure(fig, tmp_path / stem))


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
    """Two levels scattered into (oblique): a log P(v') panel; one level (hard wall): no such panel and a note."""
    oblique = drawn("results", "oblique", fixture_runs, stopped_run)
    titles = {ax.get_title(): ax for ax in oblique.axes}
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
