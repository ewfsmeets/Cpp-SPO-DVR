"""Panels of frames/: the density over (Z, r) at one (x, y) grid point (a movie frame), the density over the unit
cell in a window of (Z, r), the populations of the rotational states over Z, the J populations and the basis edges
against time."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure

from compute.frames import EdgeSeries, FrameContours
from compute.grids import CellGeometry
from panels.potential import draw_cell_outline
from panels.time_series import mark_segments
from style import (
    LINE_STYLES,
    MARKERS,
    MUTED,
    SEQUENTIAL,
    apply_scale,
    axis_label,
    masked,
    scale_and_limits,
)

MAX_STATE_LABELS = 20  # the state map labels every state up to this many


def draw_point_density(
    ax: Axes,
    z_bohr: np.ndarray,
    r_bohr: np.ndarray,
    density: np.ndarray,
    color_scale: LogNorm,
    contours: FrameContours | None,
    title: str,
    mark_zasymp: bool = True,
) -> None:
    """rho(Z, r) at one (x, y) point and time (density (nz, nr), bohr^-4; x: Z, y: r) on the log colour scale,
    values below its lower end not drawn; with contours, V_min's lines over the map (every level of
    FrameContours, the 0 eV line labelled) and, unless mark_zasymp is false (the caller marks it), zasymp, where
    they end."""
    shown = masked(density, color_scale.vmin)
    ax.pcolormesh(z_bohr, r_bohr, shown.T, shading="nearest", cmap=SEQUENTIAL, norm=color_scale, rasterized=True)
    if contours is not None:
        lines = ax.contour(
            contours.z_bohr, contours.r_bohr, contours.v_ev.T, levels=contours.levels_ev, colors=MUTED, linewidths=0.6
        )
        ax.clabel(lines, levels=[0.0], fmt="0 eV", fontsize=6)
    if contours is not None and mark_zasymp:
        ax.axvline(contours.z_asymptotic_bohr, color=MUTED, linewidth=0.8, linestyle=":")
        ax.annotate(
            " zasymp", (contours.z_asymptotic_bohr, 1.0), xycoords=("data", "axes fraction"), va="top", fontsize=6
        )
    ax.grid(False)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("r", "bohr"))
    ax.set_title(title)


def add_point_density_colorbar(fig: Figure, color_scale: LogNorm, cax: Axes) -> None:
    """The colour bar of rho(Z, r) at a point, in its own axes."""
    fig.colorbar(ScalarMappable(norm=color_scale, cmap=SEQUENTIAL), cax=cax, label="rho (bohr^-4)")


def draw_cell_density(
    ax: Axes,
    cell: CellGeometry,
    density: np.ndarray,
    color_scale: Normalize,
    marks: list[tuple[str, float, float]],
    floor: float,
    title: str,
    images: list[list[tuple[float, float]]] | None = None,
) -> None:
    """A density over the unit cell (bohr^-2): density (nx + 1, ny + 1) at the grid points with the first row and
    column repeated (periodic), each drawn as the parallelogram around it, on a linear colour scale from zero;
    values below `floor` (bohr^-2) not drawn (noise); the cell's outline, the marked points and their symmetry
    images (draw_cell_outline)."""
    shown = masked(density, floor)
    ax.pcolormesh(cell.corners_x_bohr, cell.corners_y_bohr, shown, cmap=SEQUENTIAL, norm=color_scale, rasterized=True)
    draw_cell_outline(ax, cell, marks, title, images)


def add_cell_density_colorbar(fig: Figure, ax: Axes, color_scale: Normalize, cax: Axes | None = None) -> None:
    """The colour bar of a density over the unit cell (beside ax, or in cax when given: the viewer's fixed place)."""
    mappable = ScalarMappable(norm=color_scale, cmap=SEQUENTIAL)
    label = "probability per area (bohr^-2)"
    if cax is None:
        fig.colorbar(mappable, ax=ax, label=label)
    else:
        fig.colorbar(mappable, cax=cax, label=label)


def draw_state_populations(
    ax: Axes, z_bohr: np.ndarray, populations: np.ndarray, states: np.ndarray, color_scale: LogNorm, title: str
) -> None:
    """The probability per (Z, rotational state) at one time (populations (nz, nrot) in the basis's order, which
    is mJ by mJ; x: Z, y: the states sorted by J, then mJ) on the log colour scale, values below its lower end not
    drawn; every state labelled (J, mJ) when there are at most MAX_STATE_LABELS, else the first of each J."""
    order = np.lexsort((states[:, 1], states[:, 0]))  # by J, then mJ
    sorted_states, sorted_populations = states[order], populations[:, order]
    shown = masked(sorted_populations, color_scale.vmin).T
    rows = np.arange(len(states))
    ax.pcolormesh(z_bohr, rows, shown, shading="nearest", cmap=SEQUENTIAL, norm=color_scale, rasterized=True)
    if len(states) <= MAX_STATE_LABELS:
        ax.set_yticks(rows, [f"J {j}, mJ {m}" for j, m in sorted_states])
    else:
        first_of_j = [k for k in rows if k == 0 or sorted_states[k, 0] != sorted_states[k - 1, 0]]
        ax.set_yticks(first_of_j, [f"J = {sorted_states[k, 0]}" for k in first_of_j])
    ax.grid(False)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel("rotational state (by J, then mJ)")
    ax.set_title(title)


def draw_j_populations(
    ax: Axes, times_au: np.ndarray, j_values: np.ndarray, per_j: np.ndarray, starts_au: list[float]
) -> None:
    """The probability of each J on the grids against time (per_j (records, n_J); one line style and marker per J,
    labelled): what the absorbers took is no longer counted."""
    for k, j in enumerate(j_values):
        style = LINE_STYLES[k % len(LINE_STYLES)]
        marker = MARKERS[k % len(MARKERS)]
        ax.plot(times_au, per_j[:, k], linestyle=style, marker=marker, markersize=3, label=f"J = {j}")
    apply_scale(ax, scale_and_limits([per_j]))
    ax.locator_params(axis="x", nbins=5)  # long times: fewer ticks, so that their labels do not collide
    mark_segments(ax, starts_au)
    ax.legend(loc="best", ncols=2)
    ax.set_xlabel(axis_label("t", "a.u."))
    ax.set_ylabel(axis_label("Probability on the grids"))
    ax.set_title("J populations against time\n(what the absorbers took is gone)")


def draw_edge_series(
    ax: Axes,
    times_au: np.ndarray,
    series: list[EdgeSeries],
    limit: float | None,
    starts_au: list[float],
    y_label: str,
    title: str,
) -> None:
    """Edges of the basis against time relative to the initial norm (log axis when they span more than LOG_RANGE;
    one line style and marker each, labelled), with the verdict's limit (dashed, in the legend) when given."""
    for k, edge in enumerate(series):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        ax.plot(times_au, edge.values, linestyle=style, marker=marker, markersize=3, label=edge.label)
    limits = [np.array([limit])] if limit is not None else []
    if limit is not None:
        ax.axhline(limit, color=MUTED, linestyle="--", linewidth=1.0, label=f"limit {limit:g} (heuristic)")
    apply_scale(ax, scale_and_limits([edge.values for edge in series] + limits))
    ax.locator_params(axis="x", nbins=5)
    mark_segments(ax, starts_au)
    ax.legend(loc="best")
    ax.set_xlabel(axis_label("t", "a.u."))
    ax.set_ylabel(y_label)
    ax.set_title(title)
