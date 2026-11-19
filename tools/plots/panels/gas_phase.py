"""Panels of the gas phase: the asymptotic potential with its levels, and the wave functions."""

import numpy as np
from matplotlib.axes import Axes

from compute.gas_phase import Levels
from readers.tables import AsymptoticPotential, GasPhaseStates
from style import GRID_COLOR, LINE_STYLES, MUTED, apply_scale, axis_label, masked, scale_and_limits

LEVELS_SHOWN_EV = 3.0  # levels up to this far above the potential minimum are drawn
MARGIN_EV = 0.2  # room below the potential minimum


def draw_potential_and_levels(
    ax: Axes, potential: AsymptoticPotential, levels: Levels, r_opt_start_bohr: float, j: int
) -> None:
    """V_asy(r) (eV) with the levels E(v, j) of one j as horizontal lines (dashed: not bound on the r grid), the
    r absorber shaded from ROptStart (orange: r, as the designer)."""
    ax.plot(potential.r_bohr, potential.v_ev, linestyle="-", label="V_asy(r)")
    bottom = np.min(potential.v_ev)
    top = bottom + LEVELS_SHOWN_EV
    shown = (levels.j == j) & (levels.energy_ev < top)
    for v, energy_ev, bound in zip(levels.v[shown], levels.energy_ev[shown], levels.bound[shown], strict=True):
        ax.axhline(energy_ev, linestyle="-" if bound else "--", linewidth=0.8, color=MUTED)
        ax.annotate(f" v={v}", (potential.r_bohr[-1], energy_ev), fontsize=7, va="bottom", ha="right")
    ax.axvspan(r_opt_start_bohr, potential.r_bohr[-1], color=GRID_COLOR["r"], alpha=0.15, label="r absorber")
    ax.set_ylim(bottom - MARGIN_EV, top)
    ax.set_xlabel(axis_label("r", "bohr"))
    ax.set_ylabel(axis_label("Energy", "eV"))
    ax.set_title(f"Gas phase: V_asy(r) and the levels of j = {j}")
    ax.legend(loc="lower right", fontsize=7)


def draw_wave_functions(ax: Axes, states: GasPhaseStates, j: int, count: int, r_opt_start_bohr: float) -> None:
    """|psi_v,j(r)| of the lowest `count` levels of one j (bohr^-1/2) on a log axis (style.scale_and_limits;
    values below its lower limit not drawn), so that the tails at the grid ends show; where the r absorber
    starts."""
    shown = np.flatnonzero(states.j == j)[:count]
    magnitudes = [np.abs(states.wave_functions[i]) for i in shown]
    scale = scale_and_limits(magnitudes)
    for k, (i, magnitude) in enumerate(zip(shown, magnitudes, strict=True)):
        style = LINE_STYLES[k % len(LINE_STYLES)]
        ax.plot(states.r_bohr, masked(magnitude, scale.lower), linestyle=style, label=f"v = {states.v[i]}")
    apply_scale(ax, scale)
    ax.axvspan(r_opt_start_bohr, states.r_bohr[-1], color=GRID_COLOR["r"], alpha=0.15, label="r absorber")
    ax.set_xlabel(axis_label("r", "bohr"))
    ax.set_ylabel(axis_label("|psi(r)|", "bohr^-1/2"))
    ax.set_title(f"|psi|, j = {j}: tails in the r absorber?")
    ax.legend(loc="best", fontsize=7)
