"""Diffraction (from s_matrix.tsv): the probability of each diffraction channel (n, m) and of each order, in
and out of the plane of incidence, and the scattering angles of the final states. K denotes Cartesian
parallel momenta here (k_cart_*), not the program's components along the skewed reciprocal axes."""

from dataclasses import dataclass

import numpy as np

from compute.analysis import channel_energy_of_rows, energy_of_rows
from compute.grids import incident_parallel_momentum, parallel_momentum_of
from readers.binary import BkCoefficients
from readers.run_info import RunInfo
from readers.tables import Table

# heuristic: |K x d| <= 1e-9 |K| counts as parallel to d. K is exact to ~1e-15 |K|, and a channel off the
# plane is at least ~1e-3 |K| away from it on any lattice the program supports, so any value in between works.
IN_PLANE_TOLERANCE = 1e-9


def in_plane(k_cart_x: np.ndarray, k_cart_y: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Whether Cartesian parallel momenta lie along the unit vector `direction` (either sense; K = 0 too)."""
    cross = np.abs(k_cart_x * direction[1] - k_cart_y * direction[0])
    return cross <= IN_PLANE_TOLERANCE * np.hypot(k_cart_x, k_cart_y)


def plane_of_incidence(info: RunInfo, azimuth_rad: float | None) -> tuple[np.ndarray, str]:
    """The in-plane direction: along the incident parallel momentum; at normal incidence (no such momentum)
    along `azimuth_rad` from the first lattice vector (default 0; ignored at oblique incidence). Returns the
    unit vector and a label."""
    incident = incident_parallel_momentum(info)
    length = float(np.hypot(*incident))
    if length > 0:
        return incident / length, "the plane of incidence"
    angle = azimuth_rad or 0.0
    return np.array([np.cos(angle), np.sin(angle)]), f"the plane at azimuth {np.degrees(angle):.3g} deg"


@dataclass
class Diffraction:
    """The diffraction channels and their probabilities, summed over v', j' and mJ' (n_energies first)."""

    collision_energy_ev: np.ndarray  # (n_energies,)
    channels: list[tuple[int, int]]  # (n, m), sorted
    order: np.ndarray  # (n_channels,) diffraction order
    k_cart_per_bohr: np.ndarray  # (n_channels, 2): Cartesian K
    p_channel: np.ndarray  # (n_energies, n_channels)
    in_plane: np.ndarray  # (n_channels,) bool
    plane_direction: np.ndarray  # (2,) unit vector of the in-plane direction
    plane: str  # which plane "in plane" refers to


def order_of_channel(table: Table, channel: tuple[int, int]) -> int:
    """The diffraction order of channel (n, m), as s_matrix.tsv gives it (its first row of that channel)."""
    n, m = table["n"].astype(int), table["m"].astype(int)
    return int(table["order"][(n == channel[0]) & (m == channel[1])][0])


def diffraction(table: Table, info: RunInfo, collision_energy_ev: np.ndarray, azimuth_rad: float | None) -> Diffraction:
    """Per diffraction channel: order, Cartesian K, probability per energy, and whether it is in plane."""
    energy_of_row = energy_of_rows(table, collision_energy_ev)
    n, m = table["n"].astype(int), table["m"].astype(int)
    channels = sorted(set(zip(n.tolist(), m.tolist(), strict=True)))
    k_cart_x, k_cart_y = parallel_momentum_of(
        info, np.array([c[0] for c in channels]), np.array([c[1] for c in channels])
    )
    p_channel = np.zeros((len(collision_energy_ev), len(channels)))
    for k, (a, b) in enumerate(channels):
        here = (n == a) & (m == b)
        p_channel[:, k] = np.bincount(energy_of_row[here], weights=table["P"][here], minlength=len(collision_energy_ev))
    direction, plane = plane_of_incidence(info, azimuth_rad)
    orders = np.array([order_of_channel(table, c) for c in channels])
    inside = in_plane(k_cart_x, k_cart_y, direction)
    k_cart = np.column_stack([k_cart_x, k_cart_y])
    return Diffraction(collision_energy_ev, channels, orders, k_cart, p_channel, inside, direction, plane)


def p_per_order(d: Diffraction) -> dict[int, np.ndarray]:
    """P summed over the channels of each diffraction order, per energy."""
    return {int(o): np.sum(d.p_channel[:, d.order == o], axis=1) for o in np.unique(d.order)}


def out_of_plane_fraction(d: Diffraction) -> np.ndarray:
    """The fraction of P_scattered out of the plane, per energy (NaN where nothing is scattered). At an
    azimuth off every mirror line only the specular channel is in plane, so it is 1 - P(0, 0) / P_scattered."""
    total = np.sum(d.p_channel, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0, np.sum(d.p_channel[:, ~d.in_plane], axis=1) / total, np.nan)


@dataclass
class FinalAngles:
    """The open final states (v', j', n, m) of one energy, summed over mJ' (they leave at the same angles),
    and their scattering angles."""

    v: np.ndarray
    j: np.ndarray
    n: np.ndarray
    m: np.ndarray
    theta_rad: np.ndarray  # polar angle from the surface normal, atan2(|K|, k_z)
    phi_rad: np.ndarray  # azimuth of K from the first lattice vector
    probability: np.ndarray  # summed over mJ'


def scattering_angles(table: Table, bk: BkCoefficients, info: RunInfo, energy_index: int) -> FinalAngles:
    """theta_f = atan2(|K|, k_z,f) and phi_f = atan2(K_y, K_x) of every open final state (v', j', n, m) at one
    analysis energy, with k_z,f^2 = 2 M (E_total - E_channel) (the program's channel momentum)."""
    rows = np.flatnonzero(
        energy_of_rows(table, np.array(info.data["analysis_energies"]["collision_ev"])) == energy_index
    )
    keys = np.column_stack([table[c][rows].astype(int) for c in ("v", "j", "n", "m")])
    states, first, state_of_row = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    probability = np.bincount(state_of_row.ravel(), weights=table["P"][rows], minlength=len(states))
    first_row = rows[first]  # one row of each state: all its mJ' share the channel energy
    mass_me = info.data["molecule"]["total_mass"]
    channel_hartree = channel_energy_of_rows(table, bk)[first_row]
    kz_sq_per_bohr2 = 2.0 * mass_me * (table["E_total"][first_row] / info.hartree_to_ev - channel_hartree)
    k_cart_x, k_cart_y = parallel_momentum_of(info, states[:, 2], states[:, 3])
    theta = np.arctan2(np.hypot(k_cart_x, k_cart_y), np.sqrt(np.clip(kz_sq_per_bohr2, 0.0, None)))
    v, j, n, m = states.T
    return FinalAngles(v, j, n, m, theta, np.arctan2(k_cart_y, k_cart_x), probability)
