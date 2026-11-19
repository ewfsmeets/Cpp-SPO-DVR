"""Panels of the final rovibrational states: P(v', j') against energy and at one energy, P(mJ') within each
j', and the alignment A0(2)."""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LogNorm, Normalize

from compute.states import StateToState, helicopter_limit
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

HATCHES = ["", "//", "..", "xx", "\\\\", "oo"]  # bars of different j' differ by hatch, not by colour


def level_label(v: int, j: int) -> str:
    """ "v'=0 j'=2"."""
    return f"v'={v} j'={j}"


def draw_level_probabilities(ax: Axes, s: StateToState, chosen: np.ndarray) -> None:
    """P(v', j') of the chosen levels (indices into s.levels) against the collision energy; y log when they
    span more than 1e3 (style.scale_and_limits), values below its lower limit not drawn."""
    scale = scale_and_limits([s.p_vj[:, chosen]])
    for k, level in enumerate(chosen):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        label = level_label(*s.levels[level])
        ax.plot(
            s.collision_energy_ev,
            masked(s.p_vj[:, level], scale.lower),
            linestyle=style,
            marker=marker,
            markersize=3,
            label=label,
        )
    apply_scale(ax, scale)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Probability")
    ax.set_title("Final rovibrational levels (summed over mJ' and diffraction)")
    ax.legend(loc="best", fontsize=7)


def draw_level_map(
    ax: Axes, v_values: list[int], j_values: list[int], grid: np.ndarray, energy_ev: float, j_top: int
) -> None:
    """P(v', j') at one energy over j' (x) and v' (y); colours log when the values span more than 1e3, cells
    below 1e-8 of the largest (or without probability) left blank; the top j' of the basis is marked."""
    top = float(np.nanmax(grid, initial=0.0)) if np.isfinite(grid).any() else 0.0
    if not top > 0:  # nothing scattered at this energy (or NaN): no map
        ax.set_title(f"P(v', j') at E_coll = {energy_ev:.3g} eV: nothing scattered")
        ax.set_xlabel("j'")
        ax.set_ylabel("v'")
        return
    shown = masked(grid, DENSITY_FLOOR * top)
    scale = scale_and_limits([shown])
    norm = LogNorm(scale.lower, top) if scale.kind == "log" else Normalize(0.0, top)
    mesh = ax.pcolormesh(j_values, v_values, shown, shading="nearest", cmap=SEQUENTIAL, norm=norm, rasterized=True)
    ax.figure.colorbar(mesh, ax=ax, label="Probability")
    ax.set_xticks(j_values, [f"{j}\n(JMax)" if j == j_top else str(j) for j in j_values])
    ax.set_yticks(v_values)
    ax.grid(False)
    ax.set_xlabel("j'")
    ax.set_ylabel("v'")
    ax.set_title(f"P(v', j') at E_coll = {energy_ev:.3g} eV (blank: none)")


def draw_mj_distribution(ax: Axes, fractions: dict[int, np.ndarray], energy_ev: float) -> None:
    """P(mJ') / P(j') of each j' at one energy (compute.states.mj_fractions), bars side by side, told apart
    by hatch."""
    width = 0.8 / max(len(fractions), 1)
    for k, (j, fraction) in enumerate(sorted(fractions.items())):
        offset = (k - (len(fractions) - 1) / 2) * width
        ax.bar(
            np.arange(-j, j + 1) + offset,
            fraction,
            width,
            fill=False,
            hatch=HATCHES[k % len(HATCHES)],
            label=f"j' = {j}",
        )
    ax.set_xlabel("mJ' (quantisation axis: the surface normal)")
    ax.set_ylabel("P(mJ') / P(j')")
    ax.set_title(f"mJ' distribution at E_coll = {energy_ev:.3g} eV")
    if fractions:
        ax.legend(loc="best", fontsize=7)


def draw_alignment(ax: Axes, s: StateToState) -> None:
    """A0(2)(j'; E) for j' >= 1, between -1 (cartwheel) and its helicopter limit (2j' - 1) / (j' + 1) (dotted);
    energies where P(j') is noise and j' above the run's MJMax (truncated basis) are left out."""
    drawn = [j for j in sorted(s.alignment) if not np.isnan(s.alignment[j]).all()]
    for k, j in enumerate(drawn):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        ax.plot(s.collision_energy_ev, s.alignment[j], linestyle=style, marker=marker, markersize=3, label=f"j' = {j}")
        ax.axhline(helicopter_limit(j), linestyle=":", linewidth=0.8, color=MUTED)
    ax.set_ylim(-1.05, max([helicopter_limit(j) for j in drawn], default=1.0) + 0.1)
    truncated = [j for j in sorted(s.alignment) if j > s.mj_max]
    note = f"; j' > MJMax = {s.mj_max} left out" if truncated else ""
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("A0(2)")
    ax.set_title(f"Final-state alignment (dotted: helicopter limit{note})")
    if drawn:
        ax.legend(loc="best", fontsize=7)
