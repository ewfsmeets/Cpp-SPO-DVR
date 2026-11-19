"""Panels against energy: the probabilities and the mean final energies."""

import numpy as np
from matplotlib.axes import Axes

from compute.results import Results
from style import LINE_STYLES, MARKERS, apply_scale, axis_label, masked, scale_and_limits

MEAN_ENERGY_AXIS_EV = 0.01  # the mean-energy axis spans at least 0-10 meV, so rounding noise stays flat at 0


def draw_probabilities(ax: Axes, results: Results) -> None:
    """P_scattered and P_reaction against the collision energy (linear: P_reaction is signed); the
    unitarity-error floor when P_reaction < 0."""
    ax.plot(results.collision_energy_ev, results.p_scattered, linestyle="-", label="P_scattered")
    ax.plot(results.collision_energy_ev, results.p_reaction, linestyle="--", label="P_reaction = 1 - P_scattered")
    if results.unitarity_error > 0:
        label = f"unitarity error >= {results.unitarity_error:.1e}"
        ax.axhline(-results.unitarity_error, linestyle=LINE_STYLES[3], linewidth=1, label=label)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Probability")
    ax.set_title("Scattering and reaction")
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


def draw_mean_energies(ax: Axes, results: Results) -> None:
    """The mean final vibrational, rotational and diffraction (parallel) energy of the scattered molecules
    (eV), not drawn where P_scattered < 1e-6; the axis starts at 0 and spans at least 10 meV."""
    for (name, values), style in zip(results.mean_energy_ev.items(), LINE_STYLES, strict=False):  # 3 series, 4 styles
        ax.plot(results.collision_energy_ev, values, linestyle=style, label=MEAN_ENERGY_LABELS[name])
    top = max(float(np.nanmax(np.concatenate(list(results.mean_energy_ev.values())), initial=0.0)), 0.0)
    ax.set_ylim(0.0, max(1.05 * top, MEAN_ENERGY_AXIS_EV))
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel(axis_label("Mean final energy", "eV"))
    ax.set_title("Energy disposal of the scattered molecules")
    ax.legend(loc="best", fontsize=7)
