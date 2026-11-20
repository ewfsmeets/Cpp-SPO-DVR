"""Panels of the FrameExtras datasets of frames/: the distribution of the molecular axis per steradian, probability
profiles over Z split into labelled parts (diffraction orders, vibrational levels), a rotational moment over Z at
several times, and a momentum distribution at several times with the edge of its grid."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.densities import Marks
from compute.frame_extras import ISOTROPIC_PER_SR, MomentumDistribution, StericDistribution, ZProfiles
from panels.densities import draw_marks
from style import (
    DENSITY_FLOOR,
    LINE_STYLES,
    MARKERS,
    MUTED,
    SEQUENTIAL,
    SERIES_COLORS,
    apply_scale,
    axis_label,
    legend_over_data,
    masked,
    scale_and_limits,
)

MARKS_EVERY = 12  # a line over many grid points carries a marker at every this many points (offset per line)


def plot_series(ax: Axes, x: np.ndarray, y: np.ndarray, k: int, label: str) -> None:
    """Line k of a panel: its own colour, line style and marker (identity never by colour alone), the markers of
    different lines at different points."""
    style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
    every = (k % MARKS_EVERY, MARKS_EVERY)
    color = SERIES_COLORS[k % len(SERIES_COLORS)]
    ax.plot(x, y, linestyle=style, marker=marker, markevery=every, markersize=3, color=color, label=label)


SQUARE_FRACTION = 0.6  # a node's square spans this fraction of the smallest node spacing: neighbours never touch


def node_squares(theta_deg: np.ndarray, phi_deg: np.ndarray) -> tuple[np.ndarray, float, float]:
    """The squares drawn at the nodes: their corners (ntheta * nphi, 4, 2) as (phi, theta) in degrees, theta fastest
    in phi order (as per_sr.ravel() of (ntheta, nphi)), and their half width in phi and half height in theta, both a
    SQUARE_FRACTION of the smallest spacing, so the squares shrink as the grid grows."""
    half_phi = 0.5 * SQUARE_FRACTION * 360.0 / len(phi_deg)
    spacing = np.diff(theta_deg)
    half_theta = 0.5 * SQUARE_FRACTION * (float(np.min(spacing)) if spacing.size else 90.0)
    phi, theta = (a.ravel() for a in np.meshgrid(phi_deg, theta_deg))
    dx, dy = np.array([-1, 1, 1, -1]) * half_phi, np.array([-1, -1, 1, 1]) * half_theta
    corners = np.stack([phi[:, None] + dx, theta[:, None] + dy], axis=-1)
    return corners, half_phi, half_theta


def draw_steric_distribution(ax: Axes, steric: StericDistribution, color_scale: Normalize, title: str) -> None:
    """The distribution of the axis per steradian at the DVR nodes only (x: phi, y: theta from the surface normal),
    each node a coloured square (node_squares) on a linear scale from zero: the run knows nothing between them, so
    no cells are drawn; theta from 0 to past the last node (homonuclear: 90 degrees, with a note on the mirror)."""
    corners, half_phi, half_theta = node_squares(steric.theta_deg, steric.phi_deg)
    squares = PolyCollection(corners, cmap=SEQUENTIAL, norm=color_scale, edgecolors=MUTED, linewidths=0.5, zorder=3)
    squares.set_array(steric.per_sr.ravel())
    ax.add_collection(squares)
    ax.set_xlim(steric.phi_deg[0] - 2.0 * half_phi, steric.phi_deg[-1] + 2.0 * half_phi)
    ax.set_ylim(max(steric.theta_deg[-1] + 2.0 * half_theta, 90.0 if steric.homonuclear else 180.0), 0.0)
    ax.set_xlabel(axis_label("phi, from the first lattice vector", "deg"))
    ax.set_ylabel(axis_label("theta, from the surface normal", "deg"))
    if steric.homonuclear:
        ax.annotate(
            "theta > 90: the mirror (180 - theta, phi + 180)",
            (0.5, 1.0),
            xycoords="axes fraction",
            va="bottom",
            ha="center",
            fontsize=7,
            color=MUTED,
        )
    ax.set_title(title, pad=14)


def add_steric_colorbar(fig: Figure, axes: list[Axes], color_scale: Normalize) -> None:
    """The colour bar of the distribution maps, naming the isotropic value."""
    label = f"probability per steradian (isotropic {ISOTROPIC_PER_SR:.3f})"
    fig.colorbar(ScalarMappable(norm=color_scale, cmap=SEQUENTIAL), ax=axes, label=label)


def draw_z_profiles(ax: Axes, z_bohr: np.ndarray, profiles: ZProfiles, marks: Marks, title: str) -> None:
    """Probability per Z point of each labelled part (one line each; log axis when they span more than LOG_RANGE,
    values below its lower end not drawn), with the marked Z lines (equal positions labelled together)."""
    scale = scale_and_limits(list(profiles.values))
    for k, (label, values) in enumerate(zip(profiles.labels, profiles.values, strict=True)):
        plot_series(ax, z_bohr, masked(values, scale.lower) if scale.kind == "log" else values, k, label)
    apply_scale(ax, scale)
    draw_marks(ax, marks)
    legend_over_data(ax)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("Probability per Z point"))
    ax.set_title(title)


def draw_moment(
    ax: Axes,
    z_bohr: np.ndarray,
    series: list[tuple[str, np.ndarray, np.ndarray]],
    references: list[tuple[str, float]],
    y_label: str,
    title: str,
) -> None:
    """A rotational moment against Z at several times: series of (label, moment (nz,), probability (nz,)), the
    moment not drawn where the probability at that Z is below DENSITY_FLOOR of its record's largest (a ratio of
    noise); horizontal reference levels (label, value), dashed."""
    for k, (label, moment, probability) in enumerate(series):
        plot_series(ax, z_bohr, np.where(probability >= DENSITY_FLOOR * np.max(probability), moment, np.nan), k, label)
    for label, value in references:
        ax.axhline(value, color=MUTED, linestyle="--", linewidth=0.8)
        ax.annotate(f" {label}", (0.0, value), xycoords=("axes fraction", "data"), va="bottom", fontsize=6)
    legend_over_data(ax)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(y_label)
    ax.set_title(title)


def draw_momentum_distribution(
    ax: Axes, series: list[tuple[str, MomentumDistribution]], quantity: str, title: str
) -> None:
    """Probability per momentum point at several times (one line each) against the momentum (bohr^-1), on the scale
    rule; dashed lines at +-the edge of the grid, beyond which the verdict counts the packet at the edge."""
    scale = scale_and_limits([d.probability for _, d in series])
    for k, (label, distribution) in enumerate(series):
        values = masked(distribution.probability, scale.lower) if scale.kind == "log" else distribution.probability
        plot_series(ax, distribution.k_per_bohr, values, k, label)
    edge = series[0][1].edge_per_bohr
    for side in (-edge, edge):
        ax.axvline(side, color=MUTED, linestyle="--", linewidth=0.8)
    apply_scale(ax, scale)
    legend_over_data(ax)
    ax.set_xlabel(axis_label(quantity, "bohr^-1"))
    ax.set_ylabel(axis_label("Probability per momentum point"))
    ax.set_title(title)
