"""The viewer (viewer.py, view.py) on the fixture runs, without a window: built on an Agg figure and driven through
the bound methods its widgets call; each state the test checks is drawn for real. Without a GUI backend view.py
writes the frame movies instead."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import LogNorm
from matplotlib.contour import ContourSet

import view
from compute.frames import densest_zr, xy_probability
from compute.grids import grid_point_bohr
from explanations import CONTROLS, PANELS
from readers.run_info import AU_TIME_FS
from readers.segments import joined_info, ordered_segments
from style import scale_and_limits
from viewer import HELP_PROMPT, Click, FrameViewer, click_on, open_window, self_test_steps


def new_viewer(fixture_runs: Path, run: str = "dailight_oblique") -> FrameViewer:
    """The viewer on a fixture run. Agg has no event loop: its draw_idle redraws the whole figure at once, on every
    slider set_val (about 0.4 s); skipped here for speed (as the wave-packet designer's tests do), and draw_viewer
    draws the states a test checks."""
    infos, _ = ordered_segments([fixture_runs / run])
    figure, window = open_window([i.directory for i in infos], joined_info(infos))
    figure.canvas.draw_idle = lambda *args, **kwargs: None
    return window


def draw_viewer(window: FrameViewer) -> None:
    """Draws the window for real: an error that shows only when drawn fails the test."""
    window.figure.canvas.draw()


def map_texts(window: FrameViewer) -> set[str]:
    """The labels on the cell map."""
    return {t.get_text().strip() for t in window.map_ax.texts}


@pytest.fixture(autouse=True)
def close_figures() -> None:
    """Closes the viewers' pyplot figures after each test."""
    yield
    plt.close("all")


def test_viewer_opens_on_the_first_site(fixture_runs: Path) -> None:
    """The oblique run opens at t = 0 above its first site (top, at the origin): the right title names it, its
    position and the contours; the info box the time. At t = 0 the window below zasymp holds nothing (the packet
    starts at Z0 = 11.44): the map says so and draws no cell. Its lines: Zana and zasymp as one line (both at 8
    bohr), not the specular absorber (beyond the scattering grid)."""
    window = new_viewer(fixture_runs)
    draw_viewer(window)
    assert window.point.name == "top" and window.record() == 0
    assert window.point_ax.get_title().startswith(
        "rho(Z, r) above top, (X, Y) = (0.00, 0.00) bohr, t = 0 a.u. (0.0 fs)\n"
    )
    assert "contours: V_min - min V_asy every" in window.point_ax.get_title()
    assert "| frame t = 0 a.u. (0.0 fs)" in window.info_box.get_text()
    assert window.map_ax.get_title().startswith("t = 0 a.u. (0.0 fs): nothing in the window\nZ in [-2.00, 7.92]")
    assert np.ma.count(window.map_ax.collections[0].get_array()) == 0  # every cell masked
    assert {"top (picked)", "bridge"} <= map_texts(window)
    open_markers = [line for line in window.map_ax.get_lines() if line.get_markerfacecolor() == "none"]
    assert len(open_markers) == 3 + 4  # C6v: top's three other corners, bridge's other four positions on the cell
    lines = {t.get_text().strip() for t in window.point_ax.texts}
    assert "zasymp = Zana" in lines and not any("ZSpOpt" in line for line in lines)


def test_time_slider_play_and_scale(fixture_runs: Path) -> None:
    """The time slider shows its time; a step moves one frame on, from the last back to the first; Play and Pause
    swap; the scale button switches the map to one scale over all frames (a higher top than a late frame's own)."""
    window = new_viewer(fixture_runs)
    times = window.frames.series.times_au
    window.time_slider.set_val(times[5])
    fs = times[5] * AU_TIME_FS
    assert window.point_ax.get_title().split("\n")[0].endswith(f"t = {times[5]:g} a.u. ({fs:.1f} fs)")
    window.step()
    assert window.record() == 6
    window.time_slider.set_val(times[-1])
    window.step()
    assert window.record() == 0
    window.toggle_play()
    assert window.playing and window.buttons["Play"].label.get_text() == "Pause"
    window.toggle_play()
    assert not window.playing and window.buttons["Play"].label.get_text() == "Play"
    window.time_slider.set_val(times[-1])
    own_top = window.map_ax.collections[0].norm.vmax
    window.toggle_scale()
    assert window.buttons["Scale"].label.get_text() == "Scale: all frames"
    assert window.map_ax.collections[0].norm.vmax > own_top
    draw_viewer(window)


def test_click_picks_the_cell_drawn_under_it(fixture_runs: Path) -> None:
    """A click at the centre of the drawn cell of grid point (1, 2) (iy > 0: the 60-degree skew matters) picks that
    point, named by its indices and marked; the same Cartesian point read with orthogonal axes would be another
    point; a click on a site's point picks the site (bridge); a click outside the map changes nothing."""
    window = new_viewer(fixture_runs)
    corners_x, corners_y = window.cell.corners_x_bohr, window.cell.corners_y_bohr
    # cell (ix, iy) lies between corners ix .. ix + 1 and iy .. iy + 1 (cell_geometry starts half a step before 0)
    centre = (float(corners_x[1:3, 2:4].mean()), float(corners_y[1:3, 2:4].mean()))
    window.on_click(Click(inaxes=window.map_ax, xdata=centre[0], ydata=centre[1]))
    assert (window.point.ix, window.point.iy) == (1, 2) and window.point.name == "(ix, iy) = (1, 2)"
    assert "picked (ix, iy) = (1, 2)" in map_texts(window)
    lattice, n = window.lattice, window.frames.nx
    orthogonal = (round(centre[0] / lattice.a_bohr * n) % n, round(centre[1] / lattice.b_bohr * n) % n)
    assert orthogonal != (1, 2)
    bridge = window.frames.sites[1]
    window.on_click(click_on(window, bridge.ix, bridge.iy))
    assert window.point.name == "bridge" and "bridge (picked)" in map_texts(window)
    window.on_click(Click(inaxes=window.point_ax, xdata=0.0, ydata=0.0))
    assert window.point.name == "bridge"
    x, y = grid_point_bohr(bridge.ix, bridge.iy, window.frames.nx, window.frames.ny, lattice)
    assert window.point_ax.get_title().startswith(f"rho(Z, r) above bridge, (X, Y) = ({x:.2f}, {y:.2f}) bohr")
    draw_viewer(window)


def test_window_sets_the_map(fixture_runs: Path) -> None:
    """The map shows the probability in the (Z, r) window of the range sliders (both ends included): its drawn
    values times the area equal the probability per cell of that window (compute.frames.xy_probability) at the
    frame; the drawn array repeats the first row and column (the periodic cells), cut off here. A window whose ends
    meet at the top r point holds that one point (it crashed once: an index past the grid)."""
    window = new_viewer(fixture_runs)
    times, z, r = window.frames.series.times_au, window.frames.z_bohr, window.frames.r_bohr
    window.time_slider.set_val(times[30])
    window.z_slider.set_val((z[10], z[59]))
    window.r_slider.set_val((r[5], r[29]))
    area = float(window.info.data["grids"]["xy_volume_element"])
    drawn = np.ma.filled(window.map_ax.collections[0].get_array(), 0.0).reshape(window.frames.nx + 1, -1)[:-1, :-1]
    expected = xy_probability(window.frames, 30, slice(10, 60), slice(5, 30))
    assert np.allclose(drawn * area, expected, rtol=1e-12, atol=0)
    assert f"Z in [{z[10]:.2f}, {z[59]:.2f}], r in [{r[5]:.2f}, {r[29]:.2f}] bohr" in window.map_ax.get_title()
    window.r_slider.set_val((r[-1], r[-1]))
    assert f"r in [{r[-1]:.2f}, {r[-1]:.2f}] bohr" in window.map_ax.get_title()
    draw_viewer(window)


def test_trace_follows_the_point_and_the_window(fixture_runs: Path) -> None:
    """The time trace holds the probability in the window above the picked point at every frame and the window's
    total over the 4 x 4 grid points (recomputed here from prob_zrxy, with einsum), and marks the shown frame;
    exactly the values below the scale's lower end are left out (log scale here)."""
    window = new_viewer(fixture_runs)
    times, z, r = window.frames.series.times_au, window.frames.z_bohr, window.frames.r_bohr
    window.time_slider.set_val(times[30])
    window.on_click(click_on(window, window.frames.sites[1].ix, window.frames.sites[1].iy))
    window.z_slider.set_val((z[10], z[59]))
    window.r_slider.set_val((r[5], r[29]))
    prob = np.asarray(window.frames.series.arrays["prob_zrxy"])[:, 10:60, 5:30]
    bridge = window.frames.sites[1]
    expected = (np.einsum("tzr->t", prob[:, :, :, bridge.ix, bridge.iy]), np.einsum("tzrxy->t", prob) / 16)
    scale = scale_and_limits(list(expected))
    assert scale.kind == "log"
    point_line, mean_line, now = window.trace_ax.get_lines()
    assert point_line.get_label() == "above bridge" and mean_line.get_label() == "mean per grid point"
    for line, values in zip((point_line, mean_line), expected, strict=True):
        drawn = np.asarray(line.get_ydata(), dtype=float)
        shown = values >= scale.lower
        assert np.array_equal(np.isfinite(drawn), shown) and shown.sum() > len(times) // 2
        assert np.allclose(drawn[shown], values[shown], rtol=1e-12, atol=0)
    assert list(now.get_xdata()) == [times[30], times[30]]
    draw_viewer(window)


def contour_sets(window: FrameViewer) -> list[ContourSet]:
    """The contour sets on the cell map."""
    return [c for c in window.map_ax.collections if isinstance(c, ContourSet)]


def test_log_map_and_v_min_contours(fixture_runs: Path) -> None:
    """Map switches the cell map to a log scale over 1e-8 of its top and back. V_min lays contours of V_min(x, y)
    in meV above its lowest over the cell at the window's densest (Z, r) of the frame: that cell holds the most
    probability of the window (checked by definition, einsum); the title names it and the lowest V_min there
    (from potential/); the levels lie within the corrugation recomputed here; off again, the contours go."""
    window = new_viewer(fixture_runs)
    times = window.frames.series.times_au
    window.time_slider.set_val(times[27])  # the packet at the surface (t = 2700 a.u.)
    window.toggle_log()
    norm = window.map_ax.collections[0].norm
    assert isinstance(norm, LogNorm) and norm.vmin == pytest.approx(1e-8 * norm.vmax, rel=1e-12)
    assert window.buttons["Map"].label.get_text() == "Map: log"
    window.toggle_v_min()
    z_slice, r_slice = window.windows()
    iz, ir = densest_zr(window.frames, 27, z_slice, r_slice)
    cells = np.einsum("zrxy->zr", np.asarray(window.frames.series.arrays["prob_zrxy"][27]))
    assert cells[iz, ir] == cells[z_slice, r_slice].max() and z_slice.start <= iz < z_slice.stop
    potential, reference_hartree = window.potential
    assert iz < len(potential.z_bohr)  # below zasymp: contours drawn
    v_ev = (potential.v_min_hartree[iz, ir] - reference_hartree) * window.info.hartree_to_ev
    z, r = window.frames.z_bohr[iz], window.frames.r_bohr[ir]
    expected = f"V_min at (Z, r) = ({z:.2f}, {r:.2f}): meV above {v_ev.min():.3f} eV (vs min V_asy)"
    assert window.map_ax.get_title().split("\n")[-1] == expected
    (contours,) = contour_sets(window)
    corrugation_range = 1000 * (v_ev.max() - v_ev.min())
    assert len(contours.levels) >= 2 and 0 <= min(contours.levels) and max(contours.levels) <= corrugation_range
    draw_viewer(window)
    window.toggle_v_min()
    window.toggle_log()
    assert not contour_sets(window) and not isinstance(window.map_ax.collections[0].norm, LogNorm)


def test_v_min_contours_only_where_they_mean_something(fixture_runs: Path) -> None:
    """No contours, and a title line saying why: at t = 0 (an empty window below zasymp: the packet starts beyond
    it), and for a window beyond zasymp (Z from 8.5 bohr) at the last frame, where V does not depend on (x, y)."""
    window = new_viewer(fixture_runs)
    window.toggle_v_min()
    assert window.map_ax.get_title().endswith("V_min contours: no packet in the window") and not contour_sets(window)
    z, times = window.frames.z_bohr[: window.frames.nz_scattering], window.frames.series.times_au
    window.time_slider.set_val(times[-1])
    window.z_slider.set_val((8.5, z[-1]))
    assert window.map_ax.get_title().endswith("beyond zasymp, flat over the cell") and not contour_sets(window)
    draw_viewer(window)


def test_self_test_steps_run_on_agg(fixture_runs: Path) -> None:
    """Every step of the self-test runs (the window's own methods), the last one at the frame with the most
    probability below zasymp with V_min's contours on."""
    window = new_viewer(fixture_runs)
    for step in self_test_steps(window):
        step()
    assert window.v_min_contours and len(contour_sets(window)) == 1
    draw_viewer(window)


def test_every_panel_and_control_has_its_explanation(fixture_runs: Path) -> None:
    """Hovering over each panel and control shows its text; every text has its panel or control; elsewhere the
    prompt. Explain lays the panels' texts over them and is drawn."""
    window = new_viewer(fixture_runs)
    panels = {"cell_map": window.map_ax, "point_density": window.point_ax, "time_trace": window.trace_ax}
    controls = {"time": window.time_slider.ax, "z_window": window.z_slider.ax, "r_window": window.r_slider.ax}
    controls |= {name: button.ax for name, button in window.buttons.items()}
    assert set(panels) == set(PANELS) and set(controls) == set(CONTROLS)
    for name, ax in panels.items():
        assert window.help_for(ax) == PANELS[name]
    for name, ax in controls.items():
        assert window.help_for(ax) == CONTROLS[name]
    assert window.help_for(None) == HELP_PROMPT
    window.toggle_explanations()
    assert any(t.get_text().startswith("The probability per area") for t in window.map_ax.texts)
    draw_viewer(window)


def test_without_a_window_view_writes_the_movies(
    fixture_runs: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Under Agg (no GUI backend) view.py says so and writes the frame movies of the sites into --out; --self-test
    then refuses (exit 2, naming why), and so does a run it cannot read (no frames/)."""
    assert view.main([str(fixture_runs / "dailight_oblique"), "--out", str(tmp_path)]) == 0
    assert (tmp_path / "frame_movie_site1.gif").stat().st_size > 0
    assert "no window here" in capsys.readouterr().err
    assert view.main([str(fixture_runs / "dailight_oblique"), "--self-test", str(tmp_path / "s")]) == 2
    assert "--self-test needs one" in capsys.readouterr().err
    assert view.main([str(fixture_runs / "potential_only")]) == 2
    assert "error:" in capsys.readouterr().err
