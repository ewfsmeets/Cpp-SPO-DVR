"""Derived grid quantities: the Cartesian parallel momentum K of a diffraction channel on the skewed surface
lattice (gamma the skew angle, x along the first lattice vector).

The program keeps kx, ky and the reciprocal steps delta_kx, delta_ky as components along the skewed
reciprocal axes (src/input/process.cpp, miscellaneous_setup):
    kx = K_x / sin(gamma),  ky = (cos(gamma) K_x + sin(gamma) K_y) / sin(gamma).
So channel (n, m), with kx -> kx + n delta_kx and ky -> ky + m delta_ky, has
    K = (kx sin(gamma), ky - kx cos(gamma)),  |K|^2 = kx^2 - 2 cos(gamma) kx ky + ky^2,
and |K|^2 / 2M is its parallel kinetic energy (src/analysis/analysis_basis.cpp).
"""

import numpy as np

from readers.run_info import RunInfo


def parallel_momentum_cartesian(
    n: np.ndarray,
    m: np.ndarray,
    *,
    kx_per_bohr: float,
    ky_per_bohr: float,
    delta_kx_per_bohr: float,
    delta_ky_per_bohr: float,
    skew_angle_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """(K_x, K_y) in bohr^-1 of diffraction channels (n, m): the incident parallel momentum plus n and m
    reciprocal steps, in Cartesian components (x along the first lattice vector)."""
    kx = kx_per_bohr + np.asarray(n) * delta_kx_per_bohr
    ky = ky_per_bohr + np.asarray(m) * delta_ky_per_bohr
    return kx * np.sin(skew_angle_rad), ky - kx * np.cos(skew_angle_rad)


def incident_parallel_momentum(info: RunInfo) -> np.ndarray:
    """The incident packet's Cartesian parallel momentum (2,), bohr^-1: K of channel (0, 0)."""
    k_cart_x, k_cart_y = parallel_momentum_of(info, np.array([0]), np.array([0]))
    return np.array([k_cart_x[0], k_cart_y[0]])


def parallel_momentum_of(info: RunInfo, n: np.ndarray, m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """parallel_momentum_cartesian with the run's packet and lattice (run_info.json)."""
    packet, lattice = info.data["wave_packet"], info.data["lattice"]
    return parallel_momentum_cartesian(
        n,
        m,
        kx_per_bohr=packet["kx"],
        ky_per_bohr=packet["ky"],
        delta_kx_per_bohr=lattice["delta_kx"],
        delta_ky_per_bohr=lattice["delta_ky"],
        skew_angle_rad=lattice["skew_angle_rad"],
    )
