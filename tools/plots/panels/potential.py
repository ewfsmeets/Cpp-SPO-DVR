"""Panels of the run's potential (potential/): V_min over the unit cell at one (Z, r), an elbow V(Z, r) with the
minimax path and its barrier, and V over the DVR orientations at one point. Energies in eV relative to the
gas-phase minimum, on one diverging colour scale from -cap to +cap: values above the cap are not drawn (the hatched
background; the walls reach keV), values below it get the lowest colour."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.potential import Barrier, CellGeometry
from style import DIVERGING, INK, LINE_STYLES, MARKERS, MUTED, SURFACE, axis_label, hatch_not_drawn

POTENTIAL_COLORS = DIVERGING  # zero (the gas-phase minimum) in the middle: wells and barriers apart
LABEL_BOX = {"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 1.0}  # readable on any cell


def potential_color_scale(cap_ev: float) -> Normalize:
    """The colour scale of the potential panels: from -cap to +cap eV around the gas-phase minimum."""
    return Normalize(-cap_ev, cap_ev)


def add_potential_colorbar(fig: Figure, axes: list[Axes], color_scale: Normalize) -> None:
    """One colour bar for the potential panels `axes`, made once; open at the bottom only (values below get the
    lowest colour, values above are not drawn)."""
    fig.colorbar(
        ScalarMappable(norm=color_scale, cmap=POTENTIAL_COLORS), ax=axes, label="V - min V_asy (eV)", extend="min"
    )


def not_above(values_ev: np.ndarray, color_scale: Normalize) -> np.ndarray:
    """The values with those above the colour scale's top set to NaN (not drawn)."""
    return np.where(np.asarray(values_ev) <= color_scale.vmax, values_ev, np.nan)


def barrier_title(name: str, barrier: Barrier, z_bohr: np.ndarray, r_bohr: np.ndarray, height_ev: float) -> str:
    """A title naming the point and its barrier: height (eV above the gas-phase minimum) and place, or that the
    dissociated side is not reached, or that the barrier is the VMax clamp."""
    if not np.isfinite(height_ev):
        return f"{name}: r >= r_c not reached"
    clamp = " = VMax clamp" if barrier.at_v_max else ""
    return f"{name}: barrier {height_ev:.2f} eV{clamp}\n(Z {z_bohr[barrier.iz]:.2f}, r {r_bohr[barrier.ir]:.2f} bohr)"


def draw_vmin_map(
    ax: Axes,
    cell: CellGeometry,
    values_ev: np.ndarray,
    color_scale: Normalize,
    marks: list[tuple[str, float, float]],
    title: str,
) -> None:
    """V_min over the unit cell in Cartesian bohr (x: X, y: Y): values_ev (nx + 1, ny + 1) are the grid points
    with the first row and column repeated (periodic), each drawn as the parallelogram around it, so the drawn
    region covers the outlined unit cell; values above the scale's top not drawn; the marked points labelled."""
    ax.pcolormesh(
        cell.corners_x_bohr,
        cell.corners_y_bohr,
        not_above(values_ev, color_scale),
        cmap=POTENTIAL_COLORS,
        norm=color_scale,
        rasterized=True,
    )
    ax.plot(cell.outline_x_bohr, cell.outline_y_bohr, color=INK, linewidth=0.8)
    for k, (name, x, y) in enumerate(marks):
        ax.plot(x, y, marker=MARKERS[k % len(MARKERS)], color=INK, markersize=5, linestyle="none")
        ax.annotate(name, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7, bbox=LABEL_BOX)
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
    step_ev: float,
    barrier: Barrier | None,
    title: str,
) -> None:
    """V(Z, r) (nz, nr; eV) as filled contours every step_ev up to the scale's top, with lines every 5 steps (x: Z,
    y: r); values above the top not drawn, below the bottom the lowest colour; the minimax path and its barrier
    marked when given."""
    levels = np.arange(color_scale.vmin, color_scale.vmax + 0.5 * step_ev, step_ev)
    shown = not_above(values_ev, color_scale).T
    ax.contourf(z_bohr, r_bohr, shown, levels=levels, cmap=POTENTIAL_COLORS, norm=color_scale, extend="min")
    ax.contour(z_bohr, r_bohr, shown, levels=levels[::5], colors=MUTED, linewidths=0.4)
    if barrier is not None and barrier.path:
        path = np.array(barrier.path)
        ax.plot(z_bohr[path[:, 0]], r_bohr[path[:, 1]], color=INK, linewidth=1.2, linestyle=LINE_STYLES[1])
        ax.plot(z_bohr[barrier.iz], r_bohr[barrier.ir], marker="x", color=INK, markersize=8, linestyle="none")
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("r", "bohr"))
    ax.set_title(title)


def draw_orientation(
    ax: Axes, theta_deg: np.ndarray, phi_deg: np.ndarray, values_ev: np.ndarray, color_scale: Normalize, title: str
) -> None:
    """V at the DVR nodes (ntheta, nphi; eV) at one point as a map (x: phi, y: theta; each node's cell around it),
    the nodes marked (the only orientations the run has); values above the scale's top not drawn."""
    ax.pcolormesh(
        phi_deg,
        theta_deg,
        not_above(values_ev, color_scale),
        shading="nearest",
        cmap=POTENTIAL_COLORS,
        norm=color_scale,
    )
    phi, theta = np.meshgrid(phi_deg, theta_deg)
    ax.plot(phi.ravel(), theta.ravel(), marker=".", color=INK, markersize=3, linestyle="none")
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("phi", "deg"))
    ax.set_ylabel(axis_label("theta", "deg"))
    ax.set_title(title)
