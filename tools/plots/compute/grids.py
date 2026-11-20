"""Derived grid quantities: the Cartesian parallel momentum K of a diffraction channel on the skewed surface
lattice (gamma the skew angle, x along the first lattice vector), Cartesian positions and the periodic cells of the
(x, y) grid, the r absorber's first point, and densities per volume.

The program keeps kx, ky and the reciprocal steps delta_kx, delta_ky as components along the skewed
reciprocal axes (src/input/process.cpp, miscellaneous_setup):
    kx = K_x / sin(gamma),  ky = (cos(gamma) K_x + sin(gamma) K_y) / sin(gamma).
So channel (n, m), with kx -> kx + n delta_kx and ky -> ky + m delta_ky, has
    K = (kx sin(gamma), ky - kx cos(gamma)),  |K|^2 = kx^2 - 2 cos(gamma) kx ky + ky^2,
and |K|^2 / 2M is its parallel kinetic energy (src/analysis/analysis_basis.cpp).
"""

from dataclasses import dataclass

import numpy as np

from readers.run_info import Lattice, RunInfo


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


def density_per_volume(probability: np.ndarray, cell_volume_bohr_n: float) -> np.ndarray:
    """A probability per grid cell divided by the cell's volume: the density whose integral over the cells is
    the probability (per bohr^2 for a (Z, r) cell dZ dr, per bohr for a cell dZ of a Z profile). Every
    probability per cell becomes a density here."""
    return np.asarray(probability) / cell_volume_bohr_n


def r_absorber(info: RunInfo) -> dict | None:
    """The r absorber (ROpt) of run_info.json's absorber list, or None."""
    return next((a for a in info.data["absorbers"]["list"] if a["coordinate"] == "r"), None)


def r_absorber_first_index(info: RunInfo, nr: int) -> int:
    """The first r point of the r absorber (its points in run_info.json, the program's own index): r_c, where the
    dissociated side starts; nr (no such side) without one."""
    absorber = r_absorber(info)
    if absorber is None or absorber["points"] is None:
        return nr
    return int(absorber["points"]["first"])


def cartesian_bohr(fx: np.ndarray, fy: np.ndarray, lattice: Lattice) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian (X, Y) in bohr of fractional coordinates (fx, fy) along the lattice vectors.

    Grid point (ix, iy) lies at x = ix a / nx along the first lattice vector and y = iy b / ny along the second
    (src/grids/spatial_grids.cpp: the x grid starts at 0; src/input/process.cpp: dx = a / nx), which makes the skew
    angle gamma with the first: R = (x + y cos(gamma), y sin(gamma)). The program's PES uses the inverse map: it
    splits a bond (cos(phi), sin(phi)) into (sin(gamma - phi), sin(phi)) / sin(gamma) along the two vectors
    (src/potential/dai_light_cu111_potential.cpp, value). The reciprocal axes of parallel_momentum_cartesian,
    (sin(gamma), -cos(gamma)) and (0, 1), are perpendicular to the second and the first lattice vector,
    (cos(gamma), sin(gamma)) and (1, 0), respectively.
    """
    x = np.asarray(fx) * lattice.a_bohr
    y = np.asarray(fy) * lattice.b_bohr
    return x + y * np.cos(lattice.skew_angle_rad), y * np.sin(lattice.skew_angle_rad)


@dataclass
class CellGeometry:
    """The (x, y) grid drawn in Cartesian bohr, periodically: the corners of the cells around the grid points
    ix = 0 .. nx and iy = 0 .. ny (the last row and column repeat the first), and the unit cell's outline."""

    corners_x_bohr: np.ndarray  # (nx + 2, ny + 2): cell (ix, iy) spans fractions (ix -+ 1/2) / nx, (iy -+ 1/2) / ny
    corners_y_bohr: np.ndarray
    outline_x_bohr: np.ndarray  # (5,): the unit cell, closed
    outline_y_bohr: np.ndarray


def cell_geometry(nx: int, ny: int, lattice: Lattice) -> CellGeometry:
    """The periodic cells around the (x, y) grid points (half a step before and after each, one more row and column
    so that they cover the unit cell) and the unit cell's outline."""
    fx, fy = np.meshgrid((np.arange(nx + 2) - 0.5) / nx, (np.arange(ny + 2) - 0.5) / ny, indexing="ij")
    corners_x, corners_y = cartesian_bohr(fx, fy, lattice)
    outline_x, outline_y = cartesian_bohr(np.array([0, 1, 1, 0, 0]), np.array([0, 0, 1, 1, 0]), lattice)
    return CellGeometry(corners_x, corners_y, outline_x, outline_y)


def periodic(values: np.ndarray) -> np.ndarray:
    """(nx, ny) values of the grid points with the first row and column repeated at the end: (nx + 1, ny + 1), the
    values of the cells of cell_geometry."""
    nx, ny = values.shape
    return np.asarray(values)[np.arange(nx + 1) % nx][:, np.arange(ny + 1) % ny]


def cell_volume_bohr4(grids: dict) -> float:
    """dZ dr dx dy sin(skew) of a (Z, r, x, y) cell (bohr^4), from run_info.json's grids."""
    return float(grids["dz"]) * float(grids["dr"]) * float(grids["xy_volume_element"])
