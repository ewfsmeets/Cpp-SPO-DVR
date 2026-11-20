"""Panels of the energy-resolved stationary states (energy_resolved/): A rho_E(Z, r) at one site as a map, and its
r integral against Z at several sites, the valid region marked (outside it the state is not the stationary one)."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.figure import Figure

from compute.energy_resolved import ValidRegion
from panels.frame_extras import plot_series
from style import MUTED, SEQUENTIAL, axis_label, hatch_not_drawn, legend_over_data


def mark_valid_region(ax: Axes, region: ValidRegion, horizontal_r: bool) -> None:
    """Dashed lines at the valid region's Z ends (and its r end on a map), labelled once."""
    ax.axvline(region.z_low_bohr, color=MUTED, linestyle="--", linewidth=0.8)
    ax.axvline(region.z_high_bohr, color=MUTED, linestyle="--", linewidth=0.8, label="valid region")
    if horizontal_r and np.isfinite(region.r_high_bohr):
        ax.axhline(region.r_high_bohr, color=MUTED, linestyle="--", linewidth=0.8)


def draw_rho_e_map(
    ax: Axes, z_bohr: np.ndarray, r_bohr: np.ndarray, density: np.ndarray, color_scale: Normalize, title: str
) -> None:
    """A rho_E (nz, nr; per bohr of r) as a map (x: Z, y: r) on a linear scale from zero; NaN (outside the valid
    region) not drawn, the hatching shows there."""
    ax.pcolormesh(z_bohr, r_bohr, density.T, shading="nearest", cmap=SEQUENTIAL, norm=color_scale, rasterized=True)
    hatch_not_drawn(ax)
    ax.grid(False)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("r", "bohr"))
    ax.set_title(title)


def add_rho_e_colorbar(fig: Figure, ax: Axes, color_scale: Normalize) -> None:
    """The colour bar of a density map."""
    fig.colorbar(ScalarMappable(norm=color_scale, cmap=SEQUENTIAL), ax=ax, label="A rho_E (bohr^-1)")


def draw_profiles(
    ax: Axes, z_bohr: np.ndarray, profiles: list[tuple[str, np.ndarray]], region: ValidRegion, title: str
) -> None:
    """The r integral of A rho_E against Z, one line per site (labelled), only in the valid region; 1 is the
    incident wave alone, the oscillations its interference with the reflected waves."""
    for k, (label, values) in enumerate(profiles):
        plot_series(ax, z_bohr, values, k, label)
    mark_valid_region(ax, region, horizontal_r=False)
    ax.axhline(1.0, color=MUTED, linewidth=0.8)
    legend_over_data(ax)
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel("A rho_E integrated over r")
    ax.set_title(title)
