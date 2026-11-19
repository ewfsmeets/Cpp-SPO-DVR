"""Panels of the density over (Z, r): a map at one time (a movie frame), carpets of the Z and r profiles
against time, the norm per region and the specular line against time."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure

from compute.densities import DensityGrids, Mark, Marks, RegionNorms
from panels.time_series import mark_segments
from style import (
    GRID_COLOR,
    LINE_STYLES,
    MARKERS,
    MUTED,
    SEQUENTIAL,
    SURFACE,
    apply_scale,
    axis_label,
    density_color_scale,
    masked,
    scale_and_limits,
)


def add_density_colorbar(fig: Figure, ax: Axes, color_scale: LogNorm, label: str) -> None:
    """A colour bar for `color_scale`, made once (a movie's frames redraw the axes, not the bar)."""
    fig.colorbar(ScalarMappable(norm=color_scale, cmap=SEQUENTIAL), ax=ax, label=label, extend="min")


def grouped(lines: list[Mark]) -> dict[float, list[Mark]]:
    """The marks at each position (lines at the same place share one label: "zasymp = Zana")."""
    groups: dict[float, list[Mark]] = {}
    for mark in lines:
        groups.setdefault(mark.position_bohr, []).append(mark)
    return groups


def mark_color(marks: list[Mark]) -> str:
    """An absorber's grid colour (the r absorber in the r colour), else muted, for the marks at one place."""
    grid = next((mark.grid for mark in marks if mark.grid is not None), None)
    if grid is None:
        return MUTED
    return GRID_COLOR["r" if grid == "both" else grid]


def draw_marks(ax: Axes, marks: Marks) -> None:
    """The Z lines (vertical) and r lines (horizontal) of a map, absorbers in their grid's colour, labelled at
    the top or right edge."""
    for k, (z, at_z) in enumerate(grouped(marks.z_lines).items()):
        ax.axvline(z, color=mark_color(at_z), linewidth=0.8, linestyle=LINE_STYLES[k % len(LINE_STYLES)])
        label = " = ".join(mark.label for mark in at_z)
        ax.annotate(f" {label}", (z, 1.0), xycoords=("data", "axes fraction"), rotation=90, va="top", fontsize=6)
    for r, at_r in grouped(marks.r_lines).items():
        ax.axhline(r, color=mark_color(at_r), linewidth=0.8)
        label = " = ".join(mark.label for mark in at_r)
        ax.annotate(f"{label} ", (1.0, r), xycoords=("axes fraction", "data"), ha="right", va="bottom", fontsize=6)


def draw_density_map(
    ax: Axes,
    grids: DensityGrids,
    density_per_bohr2: np.ndarray,
    color_scale: LogNorm,
    marks: Marks,
    time_au: float,
    specular_line_bohr: float,
) -> None:
    """rho(Z, r) of both grids at one time (bohr^-2; the scattering and the specular wave summed; x: Z, y: r),
    on the log colour scale `color_scale`, values below its lower end not drawn; the marks and the specular
    analysis line at this time."""
    shown = masked(density_per_bohr2, color_scale.vmin)
    ax.pcolormesh(
        grids.z_bohr, grids.r_bohr, shown.T, shading="nearest", cmap=SEQUENTIAL, norm=color_scale, rasterized=True
    )
    ax.grid(False)
    draw_marks(ax, marks)
    ax.axvline(specular_line_bohr, color=GRID_COLOR["specular"], linewidth=1.2, linestyle=":")
    ax.annotate(
        " specular analysis line",
        (specular_line_bohr, 0.0),
        xycoords=("data", "axes fraction"),
        rotation=90,
        va="bottom",
        fontsize=6,
    )
    ax.set_xlabel(axis_label("Z", "bohr"))
    ax.set_ylabel(axis_label("r", "bohr"))
    ax.set_title(f"Where is the wave packet? rho(Z, r) at t = {time_au:g} a.u.")


def draw_carpet(
    ax: Axes,
    times_au: np.ndarray,
    coordinate_bohr: np.ndarray,
    per_bohr: np.ndarray,
    name: str,
    lines: list[Mark],
    starts_au: list[float],
) -> None:
    """A profile (n_records, n_coordinate; per bohr) against time (x) and the coordinate (y) on a log colour
    scale (style.density_color_scale: values below DENSITY_FLOOR of the largest not drawn), with horizontal
    lines at the marks (equal ones share a label) and a vertical line where each restart segment starts."""
    color_scale = density_color_scale(per_bohr)
    mesh = ax.pcolormesh(
        times_au,
        coordinate_bohr,
        masked(per_bohr, color_scale.vmin).T,
        shading="nearest",
        cmap=SEQUENTIAL,
        norm=color_scale,
        rasterized=True,
    )
    ax.figure.colorbar(mesh, ax=ax, label=f"probability per bohr of {name}", extend="min")
    for position, at in grouped(lines).items():
        ax.axhline(position, color=MUTED, linewidth=0.8, linestyle=":")
        label = " = ".join(mark.label for mark in at)
        ax.annotate(f" {label}", (0.0, position), xycoords=("axes fraction", "data"), va="bottom", fontsize=6)
    mark_segments(ax, starts_au)
    ax.grid(False)
    ax.set_xlabel(axis_label("Time", "a.u."))
    ax.set_ylabel(axis_label(name, "bohr"))
    ax.set_title(f"The packet in {name} against time")


def draw_region_norms(
    ax: Axes, legend_ax: Axes, times_au: np.ndarray, regions: RegionNorms, starts_au: list[float]
) -> None:
    """The norm in each region (together the norm) and 1 - norm against time; y log when they span more than
    1e3, values below its lower limit not drawn. The legend goes into legend_ax (its axis turned off), since
    the curves fill the panel."""
    series = {**regions.regions, "1 - norm (absorbed)": regions.absorbed}
    scale = scale_and_limits(list(series.values()))
    for k, (label, values) in enumerate(series.items()):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        ax.plot(
            times_au,
            masked(values, scale.lower),
            linestyle=style,
            marker=marker,
            markevery=max(1, len(times_au) // 10),
            markersize=3,
            label=label,
        )
    apply_scale(ax, scale)
    mark_segments(ax, starts_au)
    ax.set_xlabel(axis_label("Time", "a.u."))
    ax.set_ylabel("Norm")
    ax.set_title("Where is the norm?")
    legend_ax.axis("off")
    legend_ax.legend(*ax.get_legend_handles_labels(), loc="center left", fontsize=7)


def draw_specular_line(ax: Axes, times_au: np.ndarray, line_bohr: np.ndarray) -> None:
    """The Z of the specular analysis line against time (x: time, y: Z), over a Z carpet: it follows the incoming
    packet's tail outward until it is fixed."""
    ax.step(times_au, line_bohr, where="post", color=GRID_COLOR["specular"], linewidth=1.2)
    ax.annotate(
        "specular analysis line",
        (times_au[-1], line_bohr[-1]),
        ha="right",
        va="bottom",
        fontsize=6,
        bbox={"facecolor": SURFACE, "edgecolor": "none", "alpha": 0.8, "pad": 1.0},
    )
