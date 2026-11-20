"""Panels against energy: the probabilities and the mean final energies."""

import numpy as np
from matplotlib.axes import Axes

from compute.results import Results
from style import GRID, LINE_STYLES, MARKERS, MUTED, apply_scale, axis_label, masked, scale_and_limits

MEAN_ENERGY_AXIS_EV = 0.01  # the mean-energy axis spans at least 0-10 meV, so rounding noise stays flat at 0


def draw_reaction_log(ax: Axes, results: Results) -> None:
    """P_reaction = 1 - P_scattered on a log axis, the field's usual view: the points above the unitarity error
    filled, those inside it (not distinguishable from 0) open, values <= 0 not drawn."""
    energies, p = results.collision_energy_ev, results.p_reaction
    trusted = p > results.unitarity_error
    ax.plot(energies, np.where(p > 0, p, np.nan), linestyle="-", color=MUTED, linewidth=1.0)
    if np.any(trusted):
        ax.plot(
            energies[trusted], p[trusted], linestyle="none", marker="o", markersize=3, label="above the unitarity error"
        )
    inside = (p > 0) & ~trusted
    if np.any(inside):
        ax.plot(
            energies[inside],
            p[inside],
            linestyle="none",
            marker="o",
            markersize=3,
            markerfacecolor="none",
            label="inside the unitarity error",
        )
    ax.set_yscale("log")
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("P_reaction = 1 - P_scattered")
    ax.set_title("Reaction probability (log)\n(includes flux not yet analysed and basis loss)")
    ax.legend(loc="best", fontsize=7)


def draw_packet_spectrum(ax: Axes, energies_ev: np.ndarray, spectrum_per_ev: np.ndarray) -> None:
    """The packet's energy spectrum |a(E)|^2 (per eV) at the analysis energies, on a linear axis from zero."""
    ax.plot(energies_ev, spectrum_per_ev, linestyle="-", marker="o", markersize=3)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel(axis_label("Packet spectrum |a(E)|^2", "eV^-1"))
    ax.set_title("The packet's energy spectrum\n(the S matrix divides by it)")


def draw_reaction(ax: Axes, results: Results) -> None:
    """P_reaction = 1 - P_scattered against the collision energy on its own (signed, linear) axis, with the
    band -e..e of the unitarity error e = -min P_reaction: values inside it are not distinguishable from 0."""
    ax.plot(
        results.collision_energy_ev, results.p_reaction, linestyle="-", marker="o", markersize=3, label="P_reaction"
    )
    if results.unitarity_error > 0:
        label = f"+- unitarity error {results.unitarity_error:.1e}"
        ax.axhspan(-results.unitarity_error, results.unitarity_error, color=GRID, linewidth=0, label=label)
    ax.axhline(0.0, linewidth=0.8, color=MUTED)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("P_reaction = 1 - P_scattered")
    ax.set_title("Reaction probability (linear)\n(includes flux not yet analysed and basis loss)")
    ax.legend(loc="best", fontsize=7)


def draw_vibrational_probabilities(ax: Axes, results: Results) -> None:
    """P(v') of each level scattered into against the collision energy; y log when the values span more than
    1e3 (style.scale_and_limits), values below its lower limit not drawn."""
    scale = scale_and_limits(list(results.p_vibration.values()))
    for k, (v, p) in enumerate(results.p_vibration.items()):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        values = masked(p, scale.lower)
        ax.plot(results.collision_energy_ev, values, linestyle=style, marker=marker, markersize=3, label=f"P(v' = {v})")
    apply_scale(ax, scale)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Probability")
    ax.set_title("Final vibrational state")
    ax.legend(loc="best", fontsize=7)


MEAN_ENERGY_LABELS = {
    "vibration": "<E_vib> (above v = 0)",
    "rotation": "<E_rot>",
    "diffraction": "<E_dif> (parallel, incl. the initial E_parallel)",
}


def draw_mean_energies(ax: Axes, results: Results, e_parallel_ev: float) -> None:
    """The mean final vibrational, rotational and diffraction (parallel) energy of the scattered molecules
    (eV), not drawn where P_scattered < 1e-6, and the initial parallel energy for comparison; the axis starts
    at 0 and spans at least 10 meV."""
    for (name, values), style in zip(results.mean_energy_ev.items(), LINE_STYLES, strict=False):  # 3 series, 4 styles
        ax.plot(results.collision_energy_ev, values, linestyle=style, label=MEAN_ENERGY_LABELS[name])
    if e_parallel_ev > 0:
        ax.axhline(e_parallel_ev, linewidth=0.8, color=MUTED, label=f"initial E_parallel {e_parallel_ev:g} eV")
    top = max(float(np.nanmax(np.concatenate(list(results.mean_energy_ev.values())), initial=0.0)), e_parallel_ev)
    ax.set_ylim(0.0, max(1.05 * top, MEAN_ENERGY_AXIS_EV))
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel(axis_label("Mean final energy", "eV"))
    ax.set_title("Energy disposal of the scattered molecules")
    ax.legend(loc="best", fontsize=7)


def draw_fine_comparison(ax: Axes, on_the_fly: Results, fine: Results, spectrum: Results) -> None:
    """P_scattered against the collision energy: on the fly (line), the fine analysis on its spectral grid
    (dots) and interpolated to the analysis energies (open squares)."""
    ax.plot(on_the_fly.collision_energy_ev, on_the_fly.p_scattered, linestyle="-", label="on the fly")
    ax.plot(spectrum.collision_energy_ev, spectrum.p_scattered, linestyle="none", marker=".", label="fine, its grid")
    ax.plot(
        fine.collision_energy_ev,
        fine.p_scattered,
        linestyle="none",
        marker="s",
        fillstyle="none",
        label="fine, interpolated",
    )
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Probability")
    ax.set_title("P_scattered: on-the-fly and fine analysis")
    ax.legend(loc="best", fontsize=7)


def draw_fine_difference(ax: Axes, energies_ev: np.ndarray, difference: np.ndarray) -> None:
    """Fine minus on-the-fly P_scattered at the analysis energies both have (linear: signed)."""
    ax.plot(energies_ev, difference, linestyle="-", marker="o", markersize=3)
    ax.axhline(0.0, linewidth=0.8, color=MUTED)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("P_scattered: fine - on the fly")
    ax.set_title("Do the two analyses agree?")
