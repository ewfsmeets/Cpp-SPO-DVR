"""Panels of the run's potential (potential/): a map over the unit cell (V_min at one (Z, r), or the barrier of each
(x, y) point), an elbow V(Z, r) with the minimax path and its barrier, and V over the DVR nodes at one point.
Energies in eV relative to the gas-phase minimum, on one colour scale per figure (potential_color_scale): sequential
unless a well goes deeper than a few per cent of the top, then diverging around zero; values above the scale's top
(at most the cap) are not drawn (the hatched background; the walls reach keV)."""

import numpy as np
from matplotlib import patheffects
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from compute.grids import CellGeometry
from compute.potential import AT_R_C, AT_START, UNREACHED, Barrier, PointOfInterest
from style import DIVERGING, INK, LINE_STYLES, MARKERS, MUTED, SEQUENTIAL, SURFACE, axis_label, hatch_not_drawn

LABEL_BOX = {"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 1.0}  # readable on any cell
# The minimax path: white with a dark outline, visible on the dark valleys of the sequential map and the light zero
# of the diverging one alike.
ON_ANY_COLOUR = {"color": SURFACE, "path_effects": [patheffects.withStroke(linewidth=2.6, foreground=INK)]}


SHALLOW_WELL = 0.05  # heuristic: a well shallower than this fraction of the scale's top keeps the scale sequential


def potential_color_scale(cap_ev: float, values_ev: list[np.ndarray]) -> Normalize:
    """One colour scale for all the potential panels of a figure (values_ev: what they draw, eV), its top the largest
    value but never beyond the cap (the walls reach keV; the cap when nothing finite is above zero). Sequential from
    the lowest value (at most zero) when no well is deeper than SHALLOW_WELL of the top (the usual case, e.g. a
    physisorption dip of 10 meV: the whole colour range for the barriers), else symmetric around zero (wells and
    barriers apart)."""
    finite = np.concatenate([np.ravel(v) for v in values_ev])
    finite = finite[np.isfinite(finite)]
    largest = float(np.max(finite)) if finite.size else 0.0
    top = min(cap_ev, largest) if largest > 0 else cap_ev
    lowest = min(float(np.min(finite)), 0.0) if finite.size else 0.0
    if lowest >= -SHALLOW_WELL * top:
        return Normalize(lowest, top)
    return Normalize(-top, top)


def is_diverging(color_scale: Normalize) -> bool:
    """Whether a scale of potential_color_scale is the one symmetric around zero."""
    return color_scale.vmin < 0 and color_scale.vmin == -color_scale.vmax


def potential_colormap(color_scale: Normalize) -> str:
    """Diverging (zero in the middle) for a scale symmetric around zero, else sequential."""
    return DIVERGING if is_diverging(color_scale) else SEQUENTIAL


def add_potential_colorbar(fig: Figure, axes: list[Axes], color_scale: Normalize) -> Colorbar:
    """One colour bar for the potential panels `axes`, made once; a diverging one open at the bottom (values below
    get the lowest colour); values above the top are not drawn."""
    mappable = ScalarMappable(norm=color_scale, cmap=potential_colormap(color_scale))
    extend = "min" if is_diverging(color_scale) else "neither"
    return fig.colorbar(mappable, ax=axes, label="V - min V_asy (eV)", extend=extend)


def not_above(values_ev: np.ndarray, color_scale: Normalize) -> np.ndarray:
    """The values with those above the colour scale's top set to NaN (not drawn)."""
    return np.where(np.asarray(values_ev) <= color_scale.vmax, values_ev, np.nan)


def point_label(point: PointOfInterest) -> str:
    """A point's name, or for the lowest barrier's grid point its fractions along the lattice vectors."""
    if point.name is not None:
        return point.name
    return f"lowest barrier (x {point.fraction_x:g}, y {point.fraction_y:g})"


def barrier_title(
    label: str,
    barrier: Barrier,
    z_bohr: np.ndarray,
    r_bohr: np.ndarray,
    height_ev: float,
    found_on: str = "at the run's nodes",
) -> str:
    """A title naming the point and what the minimax pass found (found_on: on which V, e.g. the run's nodes): a
    barrier (height in eV above the gas-phase minimum and place), or its highest point at r_c (no saddle before the
    dissociated side), or none (V falls from the gas-phase side), or that r_c is not reached; and when the barrier
    is the VMax clamp."""
    if barrier.kind == UNREACHED:
        return f"{label}: r >= r_c not reached"
    if barrier.kind == AT_START:
        return f"{label}: no barrier\n(V falls from the gas-phase side)"
    if barrier.kind == AT_R_C:
        return f"{label}: no saddle before r_c\n(highest point {height_ev:.2f} eV at r_c)"
    clamp = " = VMax clamp" if barrier.at_v_max else ""
    place = f"Z {z_bohr[barrier.iz]:.2f}, r {r_bohr[barrier.ir]:.2f} bohr"
    return f"{label}: barrier {height_ev:.2f} eV{clamp}\n({found_on}; {place})"


def draw_cell_map(
    ax: Axes,
    cell: CellGeometry,
    values_ev: np.ndarray,
    color_scale: Normalize,
    marks: list[tuple[str, float, float]],
    title: str,
    images: list[list[tuple[float, float]]] | None = None,
) -> None:
    """A quantity of the (x, y) grid over the unit cell in Cartesian bohr (x: X, y: Y): values_ev (nx + 1, ny + 1)
    are the grid points with the first row and column repeated (periodic), each drawn as the parallelogram around
    it, so the drawn region covers the outlined unit cell; values above the scale's top not drawn; the marked points
    (label, X, Y) labelled, and their symmetry images (draw_cell_outline)."""
    ax.pcolormesh(
        cell.corners_x_bohr,
        cell.corners_y_bohr,
        not_above(values_ev, color_scale),
        cmap=potential_colormap(color_scale),
        norm=color_scale,
        rasterized=True,
    )
    hatch_not_drawn(ax)
    draw_cell_outline(ax, cell, marks, title, images)


def draw_cell_outline(
    ax: Axes,
    cell: CellGeometry,
    marks: list[tuple[str, float, float]],
    title: str,
    images: list[list[tuple[float, float]]] | None = None,
) -> None:
    """What every map over the unit cell has around its values: the cell's outline, the marked points (label, X, Y)
    labelled on a light box (readable on any colour), equal axes in Cartesian bohr (x: X, y: Y) and the title. With
    images (one list of (X, Y) per mark: compute.grids.site_images), each mark's symmetry-equivalent points get its
    marker, open and unlabelled."""
    ax.plot(cell.outline_x_bohr, cell.outline_y_bohr, color=INK, linewidth=0.8)
    for k, (label, x, y) in enumerate(marks):
        marker = MARKERS[k % len(MARKERS)]
        for x_image, y_image in images[k] if images and k < len(images) else []:
            ax.plot(x_image, y_image, marker=marker, color=INK, markerfacecolor="none", markersize=5, linestyle="none")
        ax.plot(x, y, marker=marker, color=INK, markersize=5, linestyle="none")
        ax.annotate(label, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=6, bbox=LABEL_BOX)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel(axis_label("X", "bohr"))
    ax.set_ylabel(axis_label("Y", "bohr"))
    ax.set_title(title)


CELL_CONTOUR_LEVELS = 5  # about this many contour lines over the cell (matplotlib's MaxNLocator)


def draw_cell_contours(ax: Axes, cell: CellGeometry, values: np.ndarray, label_format: str) -> None:
    """Contour lines of a quantity over the unit cell (values (nx + 1, ny + 1) at cell's periodic grid points; NaN
    not drawn) at about CELL_CONTOUR_LEVELS round values within the data's range, labelled with label_format (e.g. "%g meV"), thin white
    lines with a dark outline (readable over any colour of a density map)."""
    low, high = float(np.nanmin(values)), float(np.nanmax(values))
    levels = MaxNLocator(CELL_CONTOUR_LEVELS).tick_values(low, high)
    levels = levels[(levels >= low) & (levels <= high)]  # the locator rounds past the data's ends
    lines = ax.contour(cell.points_x_bohr, cell.points_y_bohr, values, levels=levels, colors=SURFACE, linewidths=0.7)
    lines.set_path_effects([patheffects.withStroke(linewidth=1.8, foreground=INK)])
    ax.clabel(lines, fmt=label_format, fontsize=6, colors=INK)


def draw_elbow(
    ax: Axes,
    z_bohr: np.ndarray,
    r_bohr: np.ndarray,
    values_ev: np.ndarray,
    color_scale: Normalize,
    contours: tuple[float, int],
    barrier: Barrier | None,
    title: str,
) -> None:
    """V(Z, r) (nz, nr; eV) as filled contours every contours[0] eV up to the scale's top, with a line every
    contours[1] steps (x: Z, y: r); values above the top not drawn, below the bottom the lowest colour; the minimax
    path when given, and its highest point (x) when that is a barrier."""
    step_ev, line_every = contours
    first = step_ev * np.floor(color_scale.vmin / step_ev)  # levels at whole steps from zero
    levels = np.arange(first, color_scale.vmax + 0.5 * step_ev, step_ev)
    shown = not_above(values_ev, color_scale).T
    ax.contourf(
        z_bohr, r_bohr, shown, levels=levels, cmap=potential_colormap(color_scale), norm=color_scale, extend="min"
    )
    ax.contour(z_bohr, r_bohr, shown, levels=levels[::line_every], colors=MUTED, linewidths=0.4)
    if barrier is not None and barrier.path:
        path = np.array(barrier.path)
        ax.plot(z_bohr[path[:, 0]], r_bohr[path[:, 1]], linewidth=1.2, linestyle=LINE_STYLES[1], **ON_ANY_COLOUR)
        if barrier.saddle:
            ax.plot(z_bohr[barrier.iz], r_bohr[barrier.ir], marker="x", markersize=8, linestyle="none", **ON_ANY_COLOUR)
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("r", "bohr"))
    ax.set_title(title)


def draw_orientation(
    ax: Axes,
    edges_deg: tuple[np.ndarray, np.ndarray],
    nodes_deg: tuple[np.ndarray, np.ndarray],
    values_ev: np.ndarray,
    color_scale: Normalize,
    title: str,
) -> None:
    """V at the DVR nodes (ntheta, nphi; eV) at one point as a map (x: phi, y: theta; each node's cell between the
    edges (theta_edges, phi_edges)), the nodes (theta, phi) marked: the only orientations the run has; values above
    the scale's top not drawn."""
    theta_edges, phi_edges = edges_deg
    theta, phi = nodes_deg
    ax.pcolormesh(
        phi_edges,
        theta_edges,
        not_above(values_ev, color_scale),
        cmap=potential_colormap(color_scale),
        norm=color_scale,
    )
    phi_grid, theta_grid = np.meshgrid(phi, theta)
    ax.plot(phi_grid.ravel(), theta_grid.ravel(), marker=".", color=INK, markersize=3, linestyle="none")
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("phi, from the first lattice vector", "deg"))
    ax.set_ylabel(axis_label("theta, from the surface normal", "deg"))
    ax.set_title(title)
