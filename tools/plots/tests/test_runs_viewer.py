"""The viewer (viewer.py, view.py) on the fixture runs, without a window: built on an Agg figure and driven through
the bound methods its widgets call; each state the test checks is drawn for real. Without a GUI backend view.py
writes the frame movies instead."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

import view
from compute.frames import xy_probability
from compute.grids import grid_point_bohr
from explanations import CONTROLS, PANELS
from readers.segments import joined_info, ordered_segments
from viewer import HELP_PROMPT, Click, FrameViewer, click_on, open_window


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
    assert window.point_ax.get_title().startswith("rho(Z, r) above top, (X, Y) = (0.00, 0.00) bohr, t = 0 a.u.\n")
    assert "contours: V_min - min V_asy every" in window.point_ax.get_title()
    assert "| frame t = 0 a.u." in window.info_box.get_text()
    assert window.map_ax.get_title().startswith("t = 0 a.u.: nothing in the window\nZ in [-2.00, 7.92]")
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
    assert window.point_ax.get_title().split("\n")[0].endswith(f"t = {times[5]:g} a.u.")
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


def test_every_panel_and_control_has_its_explanation(fixture_runs: Path) -> None:
    """Hovering over each panel and control shows its text; every text has its panel or control; elsewhere the
    prompt. Explain lays the panels' texts over them and is drawn."""
    window = new_viewer(fixture_runs)
    panels = {"cell_map": window.map_ax, "point_density": window.point_ax}
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
