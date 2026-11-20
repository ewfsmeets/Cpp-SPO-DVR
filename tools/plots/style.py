"""The one look of every figure (layer 3): colours and their meanings, rcParams, axis labels, the scale rule,
the info box at the bottom of a figure, the density colour scale, and saving (a figure as PNG + PDF, a movie as
GIF + a PNG of its still).

Colours: the wave-packet designer's (tools/wave_packet_designer/plots.py; tests/test_style.py checks that they
stay equal) with the same meanings in grid and absorber panels: the scattering grid blue, the specular grid
aqua, r orange. Elsewhere series differ by colour (SERIES_COLORS, colour-blind safe), line style, marker and
direct labels, never by colour alone. The
only colour maps: viridis and cividis (sequential), PuOr (diverging). This module holds every colour of the
tool: tests/test_style.py fails on a hex colour or a colour-map name in any other module.
"""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.animation import PillowWriter
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
from matplotlib.text import Text

from readers.run_info import RunInfo

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#ffffff"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GRID_COLOR = {"scattering": BLUE, "specular": AQUA, "r": ORANGE}
SEQUENTIAL = "viridis"
SEQUENTIAL_GREY_SAFE = "cividis"
DIVERGING = "PuOr"
# The colours of a panel's series, in order (Okabe & Ito's colour-blind-safe set without its yellow, black first),
# beside a line style and marker each, so identity never rests on colour alone. Its blue, vermillion and green are
# close in hue to BLUE, ORANGE and AQUA: where a panel draws both (the gas-phase |psi| over the r absorber's faint
# band), the band carries a label and the lines differ by style.
SERIES_COLORS = [INK, "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9"]
LINE_STYLES = ["-", "--", ":", "-."]
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]

LOG_FLOOR = 1e-12  # lower limit of a log axis, relative to the largest value
DENSITY_FLOOR = 1e-8  # floor of colour maps and of ratios, relative to the largest value
NOT_DRAWN_HATCH = "////"  # the background of maps where values beyond their scale are not drawn
LOG_RANGE = 1e3  # a series spanning more than this (max / smallest positive) gets a log axis
LOG_HEADROOM = 1.5  # a log axis ends this factor above the largest value
LINEAR_HEADROOM = 1.05  # a linear axis ends this factor above the largest value
DPI = 110
REPORT_WIDTH_IN = 8.0  # 880 px: the width of GitHub's content column
INFO_LINES = 4  # title and system; energies, incidence, times; grids; files and program version (info_lines)
INFO_LINE_IN = 0.16  # height of one info line (inches)
INFO_MARGIN_IN = 0.2  # space between the info box and the axes above it


def apply_style() -> None:
    """The rcParams of every figure (the designer's: small fonts, recessive axes and grid, frameless legends)."""
    matplotlib.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "lines.linewidth": 2.0,
            "legend.frameon": False,
            "axes.prop_cycle": matplotlib.cycler(color=SERIES_COLORS),
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
            "image.cmap": SEQUENTIAL,
        }
    )


def axis_label(quantity: str, unit: str = "") -> str:
    """ "Collision energy (eV)"; a dimensionless quantity has no unit."""
    return f"{quantity} ({unit})" if unit else quantity


def legend_over_data(ax: Axes) -> None:
    """A legend placed where matplotlib finds room, on a light frame so that it stays readable where it covers
    curves (panels with many series)."""
    ax.legend(loc="best", fontsize=7, frameon=True, framealpha=0.85, facecolor=SURFACE, edgecolor=GRID)


def masked(values: np.ndarray, floor: float) -> np.ndarray:
    """The values with those below `floor` set to NaN, so that nothing lies along an axis that looks like data."""
    values = np.asarray(values, dtype=float)
    return np.where(values < floor, np.nan, values)


@dataclass
class Scale:
    """The y (or colour) scale of a panel: "log" or "linear", its limits, and how many values a log axis
    cannot show (values <= 0, counted for the info box's notes)."""

    kind: str
    lower: float
    upper: float
    not_shown: int


def scale_and_limits(arrays: list[np.ndarray], signed: bool = False) -> Scale:
    """The scale rule of every panel, on all series together: log when the positive finite values span more
    than LOG_RANGE (max / smallest positive), from max(smallest positive, max * LOG_FLOOR); otherwise linear.
    Signed quantities are always linear, over all finite values."""
    values = np.concatenate([np.ravel(a) for a in arrays]) if arrays else np.array([])
    finite = values[np.isfinite(values)]
    if signed or finite.size == 0:
        lower = float(np.min(finite)) if finite.size else 0.0
        upper = float(np.max(finite)) if finite.size else 1.0
        room = (LINEAR_HEADROOM - 1.0) * (upper - lower)  # so that the extremes do not touch the frame
        return Scale("linear", lower - room, upper + room, 0)
    positive = finite[finite > 0]
    if positive.size == 0:
        return Scale("linear", 0.0, 1.0, int(np.sum(finite <= 0)))
    largest, smallest = float(np.max(positive)), float(np.min(positive))
    if largest / smallest > LOG_RANGE:
        return Scale("log", max(smallest, largest * LOG_FLOOR), largest * LOG_HEADROOM, int(np.sum(finite <= 0)))
    return Scale("linear", min(0.0, float(np.min(finite))), largest * LINEAR_HEADROOM, 0)


def density_color_scale(values: np.ndarray) -> LogNorm:
    """A log colour scale over all values (for a movie: all its frames, fixed before the first): from
    DENSITY_FLOOR of the largest to the largest."""
    top = float(np.nanmax(values))
    return LogNorm(DENSITY_FLOOR * top, top)


def linear_color_scale(values: np.ndarray) -> Normalize:
    """A linear colour scale from zero to the largest value (to 1 when nothing is positive)."""
    top = float(np.max(values)) if np.size(values) else 0.0
    return Normalize(0.0, top if top > 0.0 else 1.0)


def hatch_not_drawn(ax: Axes) -> None:
    """A hatched background for a map whose values beyond its colour scale are not drawn, so that a gap cannot be
    taken for a value at the middle of a diverging scale."""
    ax.set_facecolor(SURFACE)
    ax.patch.set_hatch(NOT_DRAWN_HATCH)
    ax.patch.set_edgecolor(GRID)


def apply_scale(ax: Axes, scale: Scale) -> None:
    """Sets the y scale and limits of an Axes."""
    ax.set_yscale(scale.kind)
    if scale.upper > scale.lower:
        ax.set_ylim(scale.lower, scale.upper)


def new_figure(n_notes: int = 0, width_in: float = REPORT_WIDTH_IN, height_in: float = 4.5) -> Figure:
    """A figure (no pyplot) with the style applied, `height_in` for the axes and a strip below them for the info
    box (its INFO_LINES lines and `n_notes` note lines); the constrained layout keeps the axes above the strip."""
    apply_style()
    strip_in = (INFO_LINES + n_notes) * INFO_LINE_IN + INFO_MARGIN_IN
    total_in = height_in + strip_in
    fig = Figure(figsize=(width_in, total_in), layout="constrained")
    strip = strip_in / total_in  # as a fraction of the figure height
    fig.get_layout_engine().set(rect=(0.0, strip, 1.0, 1.0 - strip))
    return fig


def info_lines(
    info: RunInfo,
    sources: list[str],
    time_au: float | None = None,
    notes: tuple[str, ...] = (),
) -> list[str]:
    """The INFO_LINES lines of the info box (title, system, initial state; energies, incidence, times; grids;
    files read, program version), then one line per note."""
    d = info.data
    g = d["grids"]
    s = d["initial_state"]
    energies = d["analysis_energies"]["collision_ev"]
    energy_range = f"E_coll {energies[0]:.3g}-{energies[-1]:.3g} eV" if energies else "no analysis energies"
    packet = d["wave_packet"]
    incidence = (
        f"E_par {packet['e_parallel_ev']:g} eV at phi {packet['phi_parallel_deg']:g} deg (E_n = E_coll - E_par)"
        if packet["e_parallel_ev"] > 0
        else "normal incidence"
    )
    lines = [
        f"{d['title']} | {d['system_id']} | (v, j, mJ) = ({s['v']}, {s['j']}, {s['mJ']})",
        f"{energy_range}, {incidence} | t {info.t_start_au:g}-{info.t_end_au:g} a.u."
        + (f" | frame t = {time_au:g} a.u." if time_au is not None else ""),
        f"grids nZ {g['nz']}/{g['nz_specular']} x nR {g['nr']} x {g['nx']} x {g['ny']}, JMax {g['j_max']}, "
        f"MJMax {g['mj_max']}",
        f"files: {', '.join(sources)} | spodvr {d['program']['version']} (git {d['program']['git']})",
    ]
    return lines + [f"note: {n}" for n in notes]


def add_info_box(
    fig: Figure, info: RunInfo, sources: list[str], time_au: float | None = None, notes: tuple[str, ...] = ()
) -> Text:
    """Writes the info box into the figure's bottom strip (never over an axes) and returns its Text."""
    lines = info_lines(info, sources, time_au, notes)
    return fig.text(
        0.01, 0.01, "\n".join(lines), ha="left", va="bottom", fontsize=7, color=INK_2, parse_math=False, wrap=True
    )


MOVIE_FPS = 5  # frames per second of the GIF movies


@dataclass(frozen=True, kw_only=True)
class Movie:
    """A figure that changes with time: update(k) redraws it for frame k (0 <= k < n_frames).

    update must keep what the layout depends on (the axes, tick labels, colour bars): save_movie fixes the
    layout after the first frame. The colour scale is the caller's, fixed for all frames before the first."""

    fig: Figure
    update: Callable[[int], None]
    n_frames: int
    still: int | None = None  # the frame saved as the PNG (the movie's picture in the report); None: the last


def file_names(stem: str, fig: Figure | Movie) -> list[str]:
    """The files a figure is saved as: stem.png and stem.pdf; a movie: stem.gif and stem.png (its still: Movie.still).
    The first is the one a report shows."""
    return [f"{stem}.gif", f"{stem}.png"] if isinstance(fig, Movie) else [f"{stem}.png", f"{stem}.pdf"]


def save_movie(movie: Movie, stem: Path) -> list[Path]:
    """The movie's files (file_names: all frames as a GIF, which GitHub plays inline, and its still as a PNG);
    returns their paths. The layout is computed for the first frame and then fixed (the figure's layout engine
    stays "none" afterwards), so the axes stay put while the frames change."""
    gif, png = (stem.parent / name for name in file_names(stem.name, movie))
    writer = PillowWriter(fps=MOVIE_FPS)
    with writer.saving(movie.fig, gif, dpi=DPI):
        for k in range(movie.n_frames):
            movie.update(k)
            if k == 0:
                movie.fig.draw_without_rendering()  # the constrained layout of the first frame ...
                movie.fig.set_layout_engine("none")  # ... kept for the others
            writer.grab_frame()
    if movie.still is not None:
        movie.update(movie.still)
    movie.fig.savefig(png, metadata={"Software": None})
    return [gif, png]


def save_figure(fig: Figure, stem: Path) -> list[Path]:
    """The figure's files (file_names: PNG without metadata, so that the bytes are reproducible, and PDF);
    returns their paths."""
    png, pdf = (stem.parent / name for name in file_names(stem.name, fig))
    fig.savefig(png, metadata={"Software": None})
    fig.savefig(pdf, metadata={"Creator": None, "Producer": None, "CreationDate": None})
    return [png, pdf]
