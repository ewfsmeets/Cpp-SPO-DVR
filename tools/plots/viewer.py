"""The interactive viewer of frames/ (layer 6; view.py opens it): left, the probability over the unit cell summed over
a window of (Z, r), on a linear or log scale, optionally with V_min(x, y)'s contours at the window's densest (Z, r);
a click picks the grid point whose drawn cell holds it; right, rho(Z, r) above that point with V_min's contours
(potential/), Zana, the end of the scattering grid and the absorbers, and below it the probability in the window
above the point against time. A time slider and Play step through the frames, two range sliders set the (Z, r)
window in bohr; hovering shows the explanation of a panel or control (explanations.PANELS, CONTROLS), Explain lays
them over the panels, the key s saves the window as a PNG. Times in a.u. with fs beside them.

The window holds the state and wires the widgets to the panels of panels/frames.py; the numbers come from compute/
(the window picks frames, points, windows and scales). The only module of the tool that uses pyplot.
"""

import textwrap
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.backend_bases import MouseEvent
from matplotlib.backends import BackendFilter, backend_registry
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.widgets import Button, RangeSlider, Slider

from compute.densities import Marks, density_marks
from compute.frames import (
    FrameContours,
    contours_at,
    densest_zr,
    index_window,
    most_at_surface,
    nearest_frame,
    usable_frames,
    window_probability,
    window_probability_at,
    xy_probability,
    xy_probability_records,
    zr_probability_at,
)
from compute.grids import (
    cell_geometry,
    cell_volume_bohr4,
    density_per_volume,
    grid_point_bohr,
    grid_point_of_cell_at,
    other_images,
    periodic,
)
from compute.potential import (
    ENERGY_CAP_EV,
    FLAT_MEV,
    FRAME_CONTOUR_STEP_EV,
    cell_potential,
    corrugation_mev,
    positions_bohr,
)
from explanations import CONTROLS, PANELS
from figures import potential_with_reference
from panels.densities import draw_marks
from panels.frames import (
    add_cell_density_colorbar,
    add_point_density_colorbar,
    draw_cell_density,
    draw_point_density,
    draw_window_trace,
)
from panels.potential import draw_cell_contours
from readers.densities import DensityGrids
from readers.run_info import RunInfo, lattice_of
from style import DENSITY_FLOOR, INK, INK_2, MUTED, SURFACE, apply_style, density_color_scale, info_lines, time_text

PLAY_INTERVAL_MS = 300  # between two frames while playing (explanations.CONTROLS["Play"] says "about three a second")
HELP_WIDTH = 60  # characters per line in the help box
EXPLAIN_WIDTH = 70  # characters per line of an explanation laid over a panel
EMPTY_WINDOW = 1e-8  # a window holding less than this share of the initial norm is shown as empty (heuristic)
HELP_PROMPT = "Hover over a panel, a slider or a button for its explanation; click the cell map to pick a point."
SAVE_KEY = "s"


@dataclass(frozen=True, kw_only=True)
class Point:
    """The picked (x, y) grid point and its name (a site's, else its indices)."""

    ix: int
    iy: int
    name: str


@dataclass(frozen=True, kw_only=True)
class Click:
    """A mouse click as the viewer reads it (the fields of matplotlib's MouseEvent it uses): the self-test's."""

    inaxes: Axes
    xdata: float
    ydata: float


def has_window_backend() -> bool:
    """Whether matplotlib can open a window here (a GUI backend; Agg and the file backends cannot)."""
    return plt.get_backend().lower() not in backend_registry.list_builtin(BackendFilter.NON_INTERACTIVE)


class FrameViewer:
    """The three panels, the sliders, the buttons and the help box on one figure.

    Fixed once: frames (frames/ of the run), info, lattice, cell (the cell map's geometry), site_positions
    ((name, X, Y) of PlotSites) and site_images (their symmetry-equivalent points), potential (potential/ and min
    V_asy, or None), line_marks (the rho(Z, r) panel's lines), initial_norm. The state update() draws: point (the
    picked grid point) with point_probability, point_density and point_scale (the probability and rho at it over
    all frames and rho's log scale, set on each pick), the frame of time_slider, the windows of z_slider and
    r_slider, fixed_scale (one colour scale over all frames for the cell map), log_map (its log scale),
    v_min_contours (V_min's contours over it) and explaining; playing tells whether Play runs. window_totals caches
    the probability in each (Z, r) window over all frames (one pass over frames/ per window, not per redraw)."""

    def __init__(self, figure: Figure, runs: list[Path], info: RunInfo) -> None:
        """Reads frames/ (and potential/ for the contours) of the run's segments, lays out the panels and the
        controls on `figure`, and draws the first frame above the first site (else the origin)."""
        self.figure, self.info = figure, info
        self.frames, self.notes = usable_frames(runs)
        self.lattice = lattice_of(info)
        self.cell = cell_geometry(self.frames.nx, self.frames.ny, self.lattice)
        self.site_positions = [(p.name, x, y) for p, x, y in positions_bohr(self.frames.sites, self.lattice)]
        self.site_images = other_images([(s.fraction_x, s.fraction_y) for s in self.frames.sites], info)
        self.potential = potential_with_reference(runs)
        self.line_marks = self.point_panel_marks()
        self.initial_norm = float(np.asarray(self.frames.series.arrays["prob_zrxy"][0]).sum())
        self.map_ax = figure.add_axes((0.05, 0.42, 0.36, 0.50))
        self.map_colorbar_ax = figure.add_axes((0.42, 0.42, 0.012, 0.50))
        self.point_ax = figure.add_axes((0.55, 0.42, 0.36, 0.50))
        self.point_colorbar_ax = figure.add_axes((0.92, 0.42, 0.012, 0.50))
        self.trace_ax = figure.add_axes((0.64, 0.17, 0.27, 0.15))
        self.add_sliders()
        self.add_buttons()
        self.fixed_scale, self.log_map, self.v_min_contours = False, False, False
        self.explaining, self.playing = False, False
        self.window_cache: dict[tuple[int, int, int, int], np.ndarray] = {}
        self.help = figure.text(0.62, 0.115, self.wrapped(HELP_PROMPT), color=INK, fontsize=8, va="top")
        self.info_box = figure.text(0.01, 0.01, "", ha="left", va="bottom", fontsize=7, color=INK_2, parse_math=False)
        self.timer = figure.canvas.new_timer(interval=PLAY_INTERVAL_MS)
        self.timer.add_callback(self.step)
        figure.canvas.mpl_connect("button_press_event", self.on_click)
        figure.canvas.mpl_connect("motion_notify_event", self.on_hover)
        figure.canvas.mpl_connect("key_press_event", self.on_key)
        self.pick(self.first_point())

    def point_panel_marks(self) -> Marks:
        """The lines of the rho(Z, r) panel: density_marks within the scattering grid's Z range (the specular
        absorber lies beyond it); zasymp is one of them, so the panel does not mark it again."""
        z, r = self.frames.z_bohr, self.frames.r_bohr
        grids = DensityGrids(z, r, self.frames.nz_scattering, float(z[1] - z[0]), float(r[1] - r[0]))
        marks = density_marks(self.info, grids)
        z_end = float(z[self.frames.nz_scattering - 1])
        return Marks([m for m in marks.z_lines if m.position_bohr <= z_end], marks.r_lines)

    def add_sliders(self) -> None:
        """The time slider (on the frames' times) and the Z and r window sliders (on the grid points, in bohr; the
        Z window opens on the points below zasymp)."""
        times = self.frames.series.times_au
        z = self.frames.z_bohr[: self.frames.nz_scattering]
        r = self.frames.r_bohr
        z_below = z[max(self.potential_points() - 1, 0)]
        time_ax = self.figure.add_axes((0.12, 0.29, 0.40, 0.03))
        self.time_slider = Slider(time_ax, "t (a.u.)", times[0], times[-1], valinit=times[0], valstep=times)
        z_ax = self.figure.add_axes((0.12, 0.23, 0.40, 0.03))
        self.z_slider = RangeSlider(z_ax, "Z window (bohr)", z[0], z[-1], valinit=(z[0], z_below), valstep=z)
        r_ax = self.figure.add_axes((0.12, 0.17, 0.40, 0.03))
        self.r_slider = RangeSlider(r_ax, "r window (bohr)", r[0], r[-1], valinit=(r[0], r[-1]), valstep=r)
        for widget in (self.time_slider, self.z_slider, self.r_slider):
            widget.on_changed(self.update)

    def add_buttons(self) -> None:
        """Play, Explain, the cell map's scale (per frame, or one over all frames), linear or log, and V_min's
        contours on it (only with potential/)."""
        v_min_label = "V_min: off" if self.potential is not None else "V_min: no potential/"
        buttons = (
            ("Play", "Play", 0.06, self.toggle_play),
            ("Explain", "Explain", 0.07, self.toggle_explanations),
            ("Scale", "Scale: per frame", 0.12, self.toggle_scale),
            ("Map", "Map: linear", 0.09, self.toggle_log),
            ("V_min", v_min_label, 0.12, self.toggle_v_min),
        )
        self.buttons = {}  # kept: a widget without a reference stops responding
        left = 0.05
        for name, label, width, action in buttons:
            button = Button(self.figure.add_axes((left, 0.09, width, 0.04)), label)
            button.on_clicked(action)
            self.buttons[name] = button
            left += width + 0.01

    def potential_points(self) -> int:
        """The number of Z points below zasymp (nz_potential), where the surface acts."""
        return int(self.info.data["grids"]["nz_potential"])

    def first_point(self) -> Point:
        """The first site of PlotSites, else the cell's origin."""
        if self.frames.sites:
            site = self.frames.sites[0]
            return Point(ix=site.ix, iy=site.iy, name=site.name)
        return Point(ix=0, iy=0, name="(ix, iy) = (0, 0)")

    def point_at(self, ix: int, iy: int) -> Point:
        """The grid point (ix, iy), named after the site there if there is one."""
        name = next((s.name for s in self.frames.sites if (s.ix, s.iy) == (ix, iy)), f"(ix, iy) = ({ix}, {iy})")
        return Point(ix=ix, iy=iy, name=name)

    def pick(self, point: Point) -> None:
        """Picks a grid point: rho(Z, r) at it over all frames and its colour scale (read once per pick), then
        redraws."""
        self.point = point
        volume_bohr4 = cell_volume_bohr4(self.info.data["grids"])
        self.point_probability = zr_probability_at(self.frames, point.ix, point.iy)
        self.point_density = density_per_volume(self.point_probability, volume_bohr4)
        self.point_scale = density_color_scale(self.point_density)
        self.update()

    def record(self) -> int:
        """The frame the time slider shows."""
        return nearest_frame(self.frames.series.times_au, float(self.time_slider.val))

    def windows(self) -> tuple[slice, slice]:
        """The Z and r index windows of the range sliders (both ends included; at least one point)."""
        z_low, z_high = self.z_slider.val
        r_low, r_high = self.r_slider.val
        z = index_window(self.frames.z_bohr[: self.frames.nz_scattering], z_low, z_high)
        return z, index_window(self.frames.r_bohr, r_low, r_high)

    def window_totals(self, z_slice: slice, r_slice: slice) -> np.ndarray:
        """(records,): the probability in the (Z, r) window over all frames, read once per window."""
        key = (z_slice.start, z_slice.stop, r_slice.start, r_slice.stop)
        if key not in self.window_cache:
            self.window_cache[key] = window_probability(self.frames, z_slice, r_slice)
        return self.window_cache[key]

    def update(self, _value: object = None) -> None:
        """Redraws the three panels, the info box and (when on) the explanations for the current state."""
        self.draw_map()
        self.draw_point()
        self.draw_trace()
        if self.explaining:
            self.draw_explanations()
        t_au = float(self.frames.series.times_au[self.record()])
        lines = info_lines(self.info, ["frames/", "run_info.json"], time_au=t_au, notes=self.notes)
        self.info_box.set_text("\n".join(lines))
        self.figure.canvas.draw_idle()

    def map_scale(self, values: np.ndarray, area_bohr2: float, z_slice: slice, r_slice: slice) -> Normalize:
        """The cell map's scale up to this frame's maximum, or with fixed_scale the maximum over all frames of the
        same window (its colours then compare across time): linear from zero, or with log_map log over
        DENSITY_FLOOR of it (linear when nothing is positive)."""
        if self.fixed_scale:
            values = density_per_volume(xy_probability_records(self.frames, z_slice, r_slice), area_bohr2)
        top = float(np.max(values))
        if self.log_map and top > 0.0:
            return density_color_scale(values)
        return Normalize(0.0, max(top, np.finfo(float).tiny))

    def draw_map(self) -> None:
        """The probability over the cell in the (Z, r) window at the current frame, the picked point marked; a
        window that holds less than EMPTY_WINDOW of the initial norm is drawn empty and says so."""
        z_slice, r_slice = self.windows()
        area_bohr2 = float(self.info.data["grids"]["xy_volume_element"])
        values = density_per_volume(xy_probability(self.frames, self.record(), z_slice, r_slice), area_bohr2)
        in_window = float(self.window_totals(z_slice, r_slice)[self.record()])
        empty = in_window < EMPTY_WINDOW * self.initial_norm
        scale = self.map_scale(values, area_bohr2, z_slice, r_slice)
        self.map_ax.clear()
        self.map_colorbar_ax.clear()
        floor = np.inf if empty else DENSITY_FLOOR * scale.vmax  # an empty window: nothing drawn
        title = self.map_title(z_slice, r_slice, in_window, empty)
        labels = self.map_labels()
        draw_cell_density(self.map_ax, self.cell, periodic(values), scale, labels, floor, title, self.site_images)
        add_cell_density_colorbar(self.figure, self.map_ax, scale, cax=self.map_colorbar_ax)
        if self.v_min_contours and self.potential is not None:
            self.draw_v_min_contours(z_slice, r_slice, empty)

    def draw_v_min_contours(self, z_slice: slice, r_slice: slice, empty: bool) -> None:
        """V_min(x, y)'s contours over the map at the (Z, r) of the window holding the most probability in this
        frame (the corrugation the packet meets there), named in a title line; none for an empty window, beyond
        zasymp (V is flat there), or where V_min is flat or not finite."""
        if empty:
            self.add_map_title_line("V_min contours: no packet in the window")
            return
        potential, reference_hartree = self.potential
        iz, ir = densest_zr(self.frames, self.record(), z_slice, r_slice)
        cell_v = cell_potential(potential, reference_hartree, self.info.hartree_to_ev, iz, ir)
        if cell_v is None:
            z, r = self.frames.z_bohr[iz], self.frames.r_bohr[ir]
            self.add_map_title_line(f"V at (Z, r) = ({z:.2f}, {r:.2f}): beyond zasymp, flat over the cell")
            return
        where = f"V_min at (Z, r) = ({cell_v.z_bohr:.2f}, {cell_v.r_bohr:.2f})"
        # the contours span at most the maps' energy cap above the cell's lowest V_min (the walls left out)
        corrugation = corrugation_mev(cell_v.v_ev, ENERGY_CAP_EV)
        if corrugation is None:
            finite = np.isfinite(cell_v.v_ev).any()
            self.add_map_title_line(f"{where}: " + (f"flat (within {FLAT_MEV:g} meV)" if finite else "not finite"))
            return
        draw_cell_contours(self.map_ax, self.cell, periodic(corrugation), "%g meV")
        lowest_ev = float(np.nanmin(cell_v.v_ev))
        self.add_map_title_line(f"{where}: meV above {lowest_ev:.3f} eV (vs min V_asy)")

    def add_map_title_line(self, line: str) -> None:
        """Adds a line to the cell map's title."""
        self.map_ax.set_title(f"{self.map_ax.get_title()}\n{line}")

    def map_title(self, z_slice: slice, r_slice: slice, in_window: float, empty: bool) -> str:
        """The cell map's title: the time, the window in bohr and the probability in it."""
        z, r = self.frames.z_bohr, self.frames.r_bohr
        t_au = float(self.frames.series.times_au[self.record()])
        z_range = f"Z in [{z[z_slice.start]:.2f}, {z[z_slice.stop - 1]:.2f}]"
        r_range = f"r in [{r[r_slice.start]:.2f}, {r[r_slice.stop - 1]:.2f}] bohr"
        what = "nothing in the window" if empty else "probability over the cell"
        when = time_text(t_au, self.info.au_time_to_fs)
        return f"t = {when}: {what}\n{z_range}, {r_range}: P = {in_window:.2e}"

    def map_labels(self) -> list[tuple[str, float, float]]:
        """The sites of the cell map, the picked one labelled so ("bridge (picked)"), or the picked point added."""
        if any(name == self.point.name for name, _, _ in self.site_positions):
            return [
                (f"{name} (picked)" if name == self.point.name else name, x, y) for name, x, y in self.site_positions
            ]
        x, y = grid_point_bohr(self.point.ix, self.point.iy, self.frames.nx, self.frames.ny, self.lattice)
        return [*self.site_positions, (f"picked {self.point.name}", x, y)]

    def draw_point(self) -> None:
        """rho(Z, r) above the picked point at the current frame, on the pick's log scale, with V_min's contours
        (when potential/ exists) and the lines of line_marks."""
        self.point_ax.clear()
        self.point_colorbar_ax.clear()
        x, y = grid_point_bohr(self.point.ix, self.point.iy, self.frames.nx, self.frames.ny, self.lattice)
        t_au = float(self.frames.series.times_au[self.record()])
        when = time_text(t_au, self.info.au_time_to_fs)
        title = f"rho(Z, r) above {self.point.name}, (X, Y) = ({x:.2f}, {y:.2f}) bohr, t = {when}"
        if self.potential is not None:
            title += f"\ncontours: V_min - min V_asy every {FRAME_CONTOUR_STEP_EV:g} eV (0 eV labelled)"
        z = self.frames.z_bohr[: self.frames.nz_scattering]
        density = self.point_density[self.record()]
        draw_point_density(
            self.point_ax, z, self.frames.r_bohr, density, self.point_scale, self.contours(), title, mark_zasymp=False
        )
        draw_marks(self.point_ax, self.line_marks)
        add_point_density_colorbar(self.figure, self.point_scale, self.point_colorbar_ax)

    def draw_trace(self) -> None:
        """The probability in the (Z, r) window above the picked point against time, the cell's mean per grid point
        beside it, the shown frame marked."""
        self.trace_ax.clear()
        z_slice, r_slice = self.windows()
        at_point = window_probability_at(self.point_probability, z_slice, r_slice)
        cell_mean = self.window_totals(z_slice, r_slice) / (self.frames.nx * self.frames.ny)
        times = self.frames.series.times_au
        t_au = float(times[self.record()])
        title = (
            f"In the window above {self.point.name} against time (1000 a.u. = {1000 * self.info.au_time_to_fs:.3g} fs)"
        )
        draw_window_trace(self.trace_ax, times, (at_point, cell_mean), self.point.name, t_au, title)

    def contours(self) -> FrameContours | None:
        """V_min's contours at the picked point, or None without potential/."""
        if self.potential is None:
            return None
        potential, reference_hartree = self.potential
        z_asymptotic_bohr = self.info.data["z_asymptotic"]
        return contours_at(
            potential, reference_hartree, self.info.hartree_to_ev, self.point.ix, self.point.iy, z_asymptotic_bohr
        )

    def on_click(self, event: MouseEvent | Click) -> None:
        """A click on the cell map picks the grid point whose cell holds it (not while the toolbar zooms or pans)."""
        toolbar = getattr(self.figure.canvas, "toolbar", None)
        if event.inaxes is not self.map_ax or event.xdata is None or (toolbar is not None and toolbar.mode):
            return
        ix, iy = grid_point_of_cell_at(event.xdata, event.ydata, self.frames.nx, self.frames.ny, self.lattice)
        self.pick(self.point_at(ix, iy))

    def step(self) -> None:
        """One frame on (after the last, the first again)."""
        times = self.frames.series.times_au
        self.time_slider.set_val(times[(self.record() + 1) % len(times)])

    def toggle_play(self, _event: object = None) -> None:
        """Starts or stops stepping through the frames."""
        self.playing = not self.playing
        self.buttons["Play"].label.set_text("Pause" if self.playing else "Play")
        if self.playing:
            self.timer.start()
        else:
            self.timer.stop()
        self.figure.canvas.draw_idle()

    def toggle_scale(self, _event: object = None) -> None:
        """Switches the cell map between a scale per frame and one over all frames."""
        self.fixed_scale = not self.fixed_scale
        self.buttons["Scale"].label.set_text("Scale: all frames" if self.fixed_scale else "Scale: per frame")
        self.update()

    def toggle_log(self, _event: object = None) -> None:
        """Switches the cell map between a linear and a log colour scale."""
        self.log_map = not self.log_map
        self.buttons["Map"].label.set_text("Map: log" if self.log_map else "Map: linear")
        self.update()

    def toggle_v_min(self, _event: object = None) -> None:
        """Shows or hides V_min's contours over the cell map (nothing without potential/)."""
        if self.potential is None:
            return
        self.v_min_contours = not self.v_min_contours
        self.buttons["V_min"].label.set_text("V_min: on" if self.v_min_contours else "V_min: off")
        self.update()

    def toggle_explanations(self, _event: object = None) -> None:
        """Shows or hides the explanations over the panels."""
        self.explaining = not self.explaining
        self.update()

    def draw_explanations(self) -> None:
        """Each panel's explanation over it, in a box (the panels are redrawn without it)."""
        box = {"facecolor": SURFACE, "alpha": 0.93, "edgecolor": MUTED, "pad": 4}
        for ax, name in ((self.map_ax, "cell_map"), (self.point_ax, "point_density"), (self.trace_ax, "time_trace")):
            text = textwrap.fill(PANELS[name], EXPLAIN_WIDTH)
            ax.text(0.01, 0.98, text, transform=ax.transAxes, va="top", fontsize=7.5, color=INK, zorder=10, bbox=box)

    def help_for(self, ax: Axes | None) -> str:
        """The explanation of the panel or control that `ax` belongs to; elsewhere the prompt."""
        if ax is self.map_ax:
            return PANELS["cell_map"]
        if ax is self.point_ax:
            return PANELS["point_density"]
        if ax is self.trace_ax:
            return PANELS["time_trace"]
        for name, widget in (("time", self.time_slider), ("z_window", self.z_slider), ("r_window", self.r_slider)):
            if ax is widget.ax:
                return CONTROLS[name]
        for name, button in self.buttons.items():
            if ax is button.ax:
                return CONTROLS[name]
        return HELP_PROMPT

    def wrapped(self, text: str) -> str:
        """The text wrapped to the help box."""
        return "\n".join(textwrap.wrap(text, HELP_WIDTH))

    def on_hover(self, event: MouseEvent) -> None:
        """Shows the explanation of what the mouse is over in the help box."""
        text = self.wrapped(self.help_for(event.inaxes))
        if text != self.help.get_text():
            self.help.set_text(text)
            self.figure.canvas.draw_idle()

    def on_key(self, event: object) -> None:
        """s saves the window as viewer_t<time>.png in the working directory."""
        if getattr(event, "key", None) == SAVE_KEY:
            self.save(Path(f"viewer_t{self.frames.series.times_au[self.record()]:g}.png"))

    def save(self, path: Path) -> None:
        """Saves the window as it is now."""
        self.figure.savefig(path)


def open_window(runs: list[Path], info: RunInfo) -> tuple[Figure, FrameViewer]:
    """A window figure (14 x 9 in) with the viewer on it; the key s becomes the viewer's (ctrl+s still opens the
    toolbar's save dialog)."""
    matplotlib.rcParams["keymap.save"] = [key for key in matplotlib.rcParams["keymap.save"] if key != SAVE_KEY]
    apply_style()
    figure = plt.figure(figsize=(14.0, 9.0))
    figure.canvas.manager.set_window_title(f"spodvr frames: {runs[-1]}")
    return figure, FrameViewer(figure, runs, info)


def show() -> None:
    """Hands the open windows to the GUI until they are closed."""
    plt.show()


def second_point(window: FrameViewer) -> tuple[int, int]:
    """The grid point of the second site of PlotSites, else (1, 1)."""
    sites = window.frames.sites
    return (sites[1].ix, sites[1].iy) if len(sites) > 1 else (1, 1)


def click_on(window: FrameViewer, ix: int, iy: int) -> Click:
    """A click on the cell map exactly at grid point (ix, iy)."""
    x, y = grid_point_bohr(ix, iy, window.frames.nx, window.frames.ny, window.lattice)
    return Click(inaxes=window.map_ax, xdata=x, ydata=y)


SELF_TEST_PNGS = (
    "viewer_1_opened.png",
    "viewer_2_time.png",
    "viewer_3_clicked.png",
    "viewer_4_played.png",
    "viewer_5_v_min.png",
)


def self_test_steps(window: FrameViewer) -> list[Callable[[], None]]:
    """The self-test's steps, one per PNG: nothing (the window as opened), the time slider to the middle frame, a
    click on the second site's grid point (else (1, 1)), one Play step, V_min's contours over the cell map at the
    frame with the packet at the surface."""
    times = window.frames.series.times_au

    def middle_frame() -> None:
        """Moves the time slider to the middle frame."""
        window.time_slider.set_val(times[len(times) // 2])

    def click_second_point() -> None:
        """Clicks the cell map at the second site's grid point."""
        window.on_click(click_on(window, *second_point(window)))

    def as_opened() -> None:
        """Leaves the window as it opened."""

    def v_min_at_the_surface() -> None:
        """Moves to the frame with the most probability below zasymp and switches V_min's contours on."""
        window.time_slider.set_val(times[most_at_surface(window.frames, window.potential_points())])
        window.toggle_v_min()

    return [as_opened, middle_frame, click_second_point, window.step, v_min_at_the_surface]


def run_self_test(runs: list[Path], info: RunInfo, out: Path) -> list[Path]:
    """Opens the window, then (from a timer, once the GUI runs) drives it through the bound methods its widgets
    call (self_test_steps), saves a PNG of the window after each step into `out` (SELF_TEST_PNGS) and closes it,
    also when a step fails. Returns the PNGs written: all of them only when every step ran."""
    figure, window = open_window(runs, info)
    out.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    def drive() -> None:
        """Runs the steps, a PNG after each; closes the window whatever happens."""
        try:
            for name, step in zip(SELF_TEST_PNGS, self_test_steps(window), strict=True):
                step()
                window.save(out / name)
                saved.append(out / name)
        finally:
            plt.close(figure)

    timer = figure.canvas.new_timer(interval=500)
    timer.single_shot = True
    timer.add_callback(drive)
    timer.start()
    plt.show()
    return saved
