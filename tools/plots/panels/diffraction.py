"""Panels of diffraction: the spot map at one energy, P per order and the out-of-plane fraction against
energy, and the scattering angles of the final states."""

import numpy as np
from matplotlib.axes import Axes

from compute.diffraction import Diffraction, FinalAngles
from style import INK, LINE_STYLES, MARKERS, MUTED, apply_scale, axis_label, masked, scale_and_limits

SPOT_AREA = 250.0  # marker area (points^2) of the most probable spot or final state
SPOT_FLOOR = 1e-6  # spots and final states below this fraction of the largest are not drawn
LABEL_FLOOR = 1e-2  # spots above this fraction of the largest are labelled (n,m)
MARGIN_PER_BOHR = 0.8  # room around the channels in the spot map (bohr^-1)


def draw_spot_map(ax: Axes, d: Diffraction, energy_index: int, incident: np.ndarray) -> None:
    """The analysed diffraction channels at their Cartesian parallel momentum K (open circles), the probability
    at one energy as filled circles of area proportional to P (the strongest labelled (n,m)), the incident
    parallel momentum (arrow) and the plane of incidence (dotted)."""
    k, p = d.k_cart_per_bohr, d.p_channel[energy_index]
    ax.scatter(k[:, 0], k[:, 1], s=12, facecolors="none", edgecolors=MUTED, linewidths=0.6, label="analysed channels")
    largest = float(np.max(p))
    if largest > 0:  # else nothing is scattered (or NaN): the analysed channels only
        shown = p >= SPOT_FLOOR * largest
        ax.scatter(k[shown, 0], k[shown, 1], s=SPOT_AREA * p[shown] / largest, color=INK, zorder=3, label="P (area)")
        for (n, m), (kx, ky), prob in zip(d.channels, k, p, strict=True):
            if prob >= LABEL_FLOOR * largest:
                ax.annotate(f"({n},{m})", (kx, ky), xytext=(6, 6), textcoords="offset points", fontsize=6)
    ax.annotate("", xy=incident, xytext=(0.0, 0.0), arrowprops={"arrowstyle": "->", "color": INK})
    corners = np.vstack([k, [[0.0, 0.0]]])  # the channels and the origin of K
    ax.set_xlim(np.min(corners[:, 0]) - MARGIN_PER_BOHR, np.max(corners[:, 0]) + MARGIN_PER_BOHR)
    ax.set_ylim(np.min(corners[:, 1]) - MARGIN_PER_BOHR, np.max(corners[:, 1]) + MARGIN_PER_BOHR)
    ax.axline((0.0, 0.0), tuple(d.plane_direction), linestyle=":", linewidth=1, color=MUTED, label=d.plane)
    ax.set_aspect("equal")
    ax.set_xlabel(axis_label("K_x", "bohr^-1"))
    ax.set_ylabel(axis_label("K_y", "bohr^-1"))
    incidence = "arrow: incidence" if np.hypot(*incident) > 0 else "normal incidence"
    ax.set_title(f"Diffraction at E_coll = {d.collision_energy_ev[energy_index]:.3g} eV ({incidence})")
    ax.legend(loc="lower left", fontsize=7)


def draw_orders(ax: Axes, energies_ev: np.ndarray, p_order: dict[int, np.ndarray]) -> None:
    """P summed over each diffraction order against the collision energy; y log when the orders span more
    than 1e3, values below its lower limit not drawn."""
    scale = scale_and_limits(list(p_order.values()))
    for k, (order, p) in enumerate(sorted(p_order.items())):
        style, marker = LINE_STYLES[k % len(LINE_STYLES)], MARKERS[k % len(MARKERS)]
        ax.plot(
            energies_ev, masked(p, scale.lower), linestyle=style, marker=marker, markersize=3, label=f"order {order}"
        )
    apply_scale(ax, scale)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Probability")
    ax.set_title("Diffraction orders")
    ax.legend(loc="best", fontsize=7)


def draw_out_of_plane(ax: Axes, energies_ev: np.ndarray, fraction: np.ndarray, plane: str, only_specular: bool) -> None:
    """The fraction of P_scattered out of the plane against the collision energy (0-1); only_specular: the
    plane holds no channel but (0, 0), so the fraction is 1 - P(0, 0) / P_scattered."""
    ax.plot(energies_ev, fraction, linestyle="-", marker="o", markersize=3)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(axis_label("Collision energy", "eV"))
    ax.set_ylabel("Fraction of P_scattered")
    note = " (only (0,0) is in it)" if only_specular else ""
    ax.set_title(f"Out of {plane}?{note}")


def draw_scattering_angles(ax: Axes, angles: FinalAngles, energy_ev: float, specular: tuple[float, float]) -> None:
    """On polar axes: every open final state (v', j', n, m) at one energy, summed over mJ', at azimuth phi_f and
    radius theta_f (degrees from the surface normal), marker area proportional to P; the specular direction
    (phi, theta) of this energy as an open circle. States below 1e-6 of the largest are not drawn."""
    largest = float(np.max(angles.probability, initial=0.0))
    shown = angles.probability >= SPOT_FLOOR * largest if largest > 0 else np.zeros(len(angles.probability), bool)
    theta_deg = np.degrees(angles.theta_rad[shown])
    area = SPOT_AREA * angles.probability[shown] / largest if largest > 0 else []
    ax.scatter(angles.phi_rad[shown], theta_deg, s=area, color=INK, alpha=0.6, zorder=3, label="final states (area: P)")
    specular_deg = np.degrees(specular[1])
    ax.scatter([specular[0]], [specular_deg], s=120, facecolors="none", edgecolors=INK, zorder=4, label="specular")
    outermost_deg = max(float(np.max(theta_deg, initial=0.0)), specular_deg)
    ax.set_ylim(0.0, min(90.0, 10.0 * np.ceil((outermost_deg + 5.0) / 10.0)))
    if shown.any():
        ax.set_rlabel_position(np.degrees(np.median(angles.phi_rad[shown])) + 180.0)  # radius labels opposite the data
    ax.set_xlabel("azimuth phi_f (deg, from the first lattice vector, the x axis of the spot map)")
    ax.set_ylabel("theta_f (deg)", labelpad=28)
    ax.set_title(f"Scattering angles at E_coll = {energy_ev:.3g} eV (summed over mJ')")
    ax.legend(loc="lower right", fontsize=7)
