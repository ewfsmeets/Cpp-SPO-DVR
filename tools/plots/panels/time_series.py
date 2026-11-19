"""Panels against time: the averages and the Balint-Kurti coefficients."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

from compute.time_series import ScatteringInTime, convergence_deviation
from style import (
    DENSITY_FLOOR,
    LINE_STYLES,
    MARKERS,
    MUTED,
    SEQUENTIAL,
    apply_scale,
    axis_label,
    masked,
    scale_and_limits,
)

DEVIATION_RANGE = (1e-4, 1.0)  # the convergence map's colour scale; values beyond get the end colours (arrows)


def mark_segments(ax: Axes, starts_au: list[float]) -> None:
    """A thin vertical line where each restart segment after the first starts."""
    for t in starts_au[1:]:
        ax.axvline(t, color=MUTED, linewidth=0.8)
        ax.annotate(" restart", (t, 1.0), xycoords=("data", "axes fraction"), va="top", fontsize=7, color=MUTED)


def draw_norm(ax: Axes, times_au: np.ndarray, norm: np.ndarray, rest_norm: np.ndarray, starts_au: list[float]) -> None:
    """Total norm, 1 - norm (absorbed) and the rest norm against time; y log when they span more than 1e3
    (style.scale_and_limits), values below its lower limit not drawn."""
    series = {
        "norm": norm,
        "1 - norm (absorbed)": 1.0 - norm,
        "rest norm (specular wave between the analysis lines)": rest_norm,
    }
    scale = scale_and_limits(list(series.values()))
    for (label, values), style in zip(series.items(), LINE_STYLES, strict=False):  # 3 series, 4 styles
        ax.plot(times_au, masked(values, scale.lower) if scale.kind == "log" else values, linestyle=style, label=label)
    apply_scale(ax, scale)
    mark_segments(ax, starts_au)
    ax.set_xlabel(axis_label("Time", "a.u."))
    ax.set_ylabel("Norm")
    ax.set_title("Where is the wave packet? Norms")
    ax.legend(loc="best")


def draw_mean_positions(
    ax: Axes, times_au: np.ndarray, z_mean_bohr: np.ndarray, r_mean_bohr: np.ndarray, starts_au: list[float]
) -> None:
    """<Z> and <r> against time (bohr)."""
    ax.plot(times_au, z_mean_bohr, linestyle="-", label="<Z>")
    ax.plot(times_au, r_mean_bohr, linestyle="--", label="<r>")
    mark_segments(ax, starts_au)
    ax.set_xlabel(axis_label("Time", "a.u."))
    ax.set_ylabel(axis_label("Mean position", "bohr"))
    ax.set_title("Mean positions")
    ax.legend(loc="best")


def draw_channel_densities(ax: Axes, times_au: np.ndarray, densities: np.ndarray, labels: list[str]) -> None:
    """|c_f(t)|^2 of the given channels at the analysis line (bohr^-1) against time; y log when they span more
    than 1e3, values below its lower limit (1e-12 of the largest, or zero) not drawn. Each channel has its own
    line style and marker."""
    scale = scale_and_limits([densities])
    every = max(1, len(times_au) // 12)  # about 12 markers per line
    for k, label in enumerate(labels):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        values = masked(densities[:, k], scale.lower)
        ax.plot(
            times_au, values, linestyle=style, marker=marker, markevery=every, markersize=4, linewidth=1.4, label=label
        )
    apply_scale(ax, scale)
    ax.set_xlabel(axis_label("Time", "a.u."))
    ax.set_ylabel(axis_label("|c(t)|^2 at the analysis line", "bohr^-1"))
    ax.set_title("When does each channel pass the analysis line?")
    ax.legend(loc="best", fontsize=7)


def draw_convergence(ax: Axes, s: ScatteringInTime) -> None:
    """|P_scattered(E; t) / P_scattered(E; tEnd) - 1| over time and energy on a log colour scale (1e-4 .. 1):
    how far the time sums stopped at t are from their final value; energies where P_scattered(tEnd) is below
    1e-8 of its maximum are not drawn. A colour still changing near tEnd means not converged."""
    final = s.p_scattered[-1]
    deviation = np.where((final > DENSITY_FLOOR * np.max(final))[:, None], convergence_deviation(s), np.nan)
    norm = LogNorm(*DEVIATION_RANGE)
    mesh = ax.pcolormesh(
        s.times_au, s.collision_energy_ev, deviation, shading="nearest", cmap=SEQUENTIAL, norm=norm, rasterized=True
    )
    ax.grid(False)
    ax.figure.colorbar(mesh, ax=ax, extend="both", label="|P(E; t) / P(E; tEnd) - 1|")
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.set_xlabel(axis_label("Time the sums stop", "a.u."))
    ax.set_ylabel(axis_label("Collision energy", "eV"))
    ax.set_title("Converged in time? (dark: close to the final value)")


def draw_energy_series(ax: Axes, energies_ev: np.ndarray, values: np.ndarray, labels: list[str], signed: bool) -> None:
    """One line per column of values (n_energies, n_series) against the collision energy, each with its own
    line style and marker; y by style.scale_and_limits (always linear when signed), NaN not drawn."""
    scale = scale_and_limits([values], signed=signed)
    for k, label in enumerate(labels):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        shown = values[:, k] if signed else masked(values[:, k], scale.lower)
        ax.plot(energies_ev, shown, linestyle=style, marker=marker, markersize=3, linewidth=1.4, label=label)
    apply_scale(ax, scale)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.legend(loc="best", fontsize=7)


def draw_time_delays(ax: Axes, energies_ev: np.ndarray, tau_au: np.ndarray, labels: list[str]) -> None:
    """tau_f(E) of the given channels (a.u.), NaN (little weight, or closed) not drawn: negative values mean
    the molecule turned before reaching Z = 0."""
    draw_energy_series(ax, energies_ev, tau_au, labels, signed=True)
    ax.set_ylabel(axis_label("Time delay relative to reflection at Z = 0", "a.u."))
    ax.set_title("Time delay d arg S / dE (S referred to its analysis line)")


def draw_channel_probabilities(ax: Axes, energies_ev: np.ndarray, probability: np.ndarray, labels: list[str]) -> None:
    """|S_f|^2 of the given channels against the collision energy (y log when they span more than 1e3)."""
    draw_energy_series(ax, energies_ev, probability, labels, signed=False)
    ax.set_ylabel("Probability")
    ax.set_title("Their probabilities")
