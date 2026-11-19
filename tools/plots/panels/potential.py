"""Panels of the run's potential (potential/): a map over the unit cell (V_min at one (Z, r), or the barrier of each
(x, y) point), an elbow V(Z, r) with the minimax path and its barrier, and V over the DVR nodes at one point.
Energies in eV relative to the gas-phase minimum, on one diverging colour scale from -cap to +cap: values above the
cap are not drawn (the hatched background; the walls reach keV), values below it get the lowest colour."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.grids import CellGeometry
from compute.potential import AT_R_C, AT_START, UNREACHED, Barrier, PointOfInterest
from style import DIVERGING, INK, LINE_STYLES, MARKERS, MUTED, SURFACE, axis_label, hatch_not_drawn

POTENTIAL_COLORS = DIVERGING  # zero (the gas-phase minimum) in the middle: wells and barriers apart
LABEL_BOX = {"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 1.0}  # readable on any cell


def potential_color_scale(cap_ev: float) -> Normalize:
    """The colour scale of the potential panels: from -cap to +cap eV around the gas-phase minimum."""
    return Normalize(-cap_ev, cap_ev)


def add_potential_colorbar(fig: Figure, axes: list[Axes], color_scale: Normalize) -> Colorbar:
    """One colour bar for the potential panels `axes`, made once; open at the bottom only (values below get the
    lowest colour, values above are not drawn)."""
    mappable = ScalarMappable(norm=color_scale, cmap=POTENTIAL_COLORS)
    return fig.colorbar(mappable, ax=axes, label="V - min V_asy (eV)", extend="min")


def not_above(values_ev: np.ndarray, color_scale: Normalize) -> np.ndarray:
    """The values with those above the colour scale's top set to NaN (not drawn)."""
    return np.where(np.asarray(values_ev) <= color_scale.vmax, values_ev, np.nan)


def point_label(point: PointOfInterest) -> str:
    """A point's name, or for the lowest barrier's grid point its fractions along the lattice vectors."""
    if point.name is not None:
        return point.name
    return f"lowest barrier (x {point.fraction_x:g}, y {point.fraction_y:g})"


def barrier_title(label: str, barrier: Barrier, z_bohr: np.ndarray, r_bohr: np.ndarray, height_ev: float) -> str:
    """A title naming the point and what the minimax pass found on the run's nodes: a barrier (height in eV above
    the gas-phase minimum and place), or its highest point at r_c (no saddle before the dissociated side), or none
    (V falls from the gas-phase side), or that r_c is not reached; and when the barrier is the VMax clamp."""
    if barrier.kind == UNREACHED:
        return f"{label}: r >= r_c not reached"
    if barrier.kind == AT_START:
        return f"{label}: no barrier\n(V falls from the gas-phase side)"
    if barrier.kind == AT_R_C:
        return f"{label}: no saddle before r_c\n(highest point {height_ev:.2f} eV at r_c)"
    clamp = " = VMax clamp" if barrier.at_v_max else ""
    place = f"Z {z_bohr[barrier.iz]:.2f}, r {r_bohr[barrier.ir]:.2f} bohr"
    return f"{label}: barrier {height_ev:.2f} eV{clamp}\n(at the run's nodes; {place})"


def draw_cell_map(
    ax: Axes,
    cell: CellGeometry,
    values_ev: np.ndarray,
    color_scale: Normalize,
    marks: list[tuple[str, float, float]],
    title: str,
) -> None:
    """A quantity of the (x, y) grid over the unit cell in Cartesian bohr (x: X, y: Y): values_ev (nx + 1, ny + 1)
    are the grid points with the first row and column repeated (periodic), each drawn as the parallelogram around
    it, so the drawn region covers the outlined unit cell; values above the scale's top not drawn; the marked points
    (label, X, Y) labelled."""
    ax.pcolormesh(
        cell.corners_x_bohr,
        cell.corners_y_bohr,
        not_above(values_ev, color_scale),
        cmap=POTENTIAL_COLORS,
        norm=color_scale,
        rasterized=True,
    )
    ax.plot(cell.outline_x_bohr, cell.outline_y_bohr, color=INK, linewidth=0.8)
    for k, (label, x, y) in enumerate(marks):
        ax.plot(x, y, marker=MARKERS[k % len(MARKERS)], color=INK, markersize=5, linestyle="none")
        ax.annotate(label, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=6, bbox=LABEL_BOX)
    hatch_not_drawn(ax)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel(axis_label("X", "bohr"))
    ax.set_ylabel(axis_label("Y", "bohr"))
    ax.set_title(title)


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
    levels = np.arange(color_scale.vmin, color_scale.vmax + 0.5 * step_ev, step_ev)
    shown = not_above(values_ev, color_scale).T
    ax.contourf(z_bohr, r_bohr, shown, levels=levels, cmap=POTENTIAL_COLORS, norm=color_scale, extend="min")
    ax.contour(z_bohr, r_bohr, shown, levels=levels[::line_every], colors=MUTED, linewidths=0.4)
    if barrier is not None and barrier.path:
        path = np.array(barrier.path)
        ax.plot(z_bohr[path[:, 0]], r_bohr[path[:, 1]], color=INK, linewidth=1.2, linestyle=LINE_STYLES[1])
        if barrier.saddle:
            ax.plot(z_bohr[barrier.iz], r_bohr[barrier.ir], marker="x", color=INK, markersize=8, linestyle="none")
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
    ax.pcolormesh(phi_edges, theta_edges, not_above(values_ev, color_scale), cmap=POTENTIAL_COLORS, norm=color_scale)
    phi_grid, theta_grid = np.meshgrid(phi, theta)
    ax.plot(phi_grid.ravel(), theta_grid.ravel(), marker=".", color=INK, markersize=3, linestyle="none")
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("phi, from the first lattice vector", "deg"))
    ax.set_ylabel(axis_label("theta, from the surface normal", "deg"))
    ax.set_title(title)
