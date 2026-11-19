"""Diffraction geometry on a hand-made lattice: the channels in the plane of incidence on a hexagonal cell."""

import numpy as np

from compute.diffraction import in_plane
from compute.grids import parallel_momentum_cartesian

GAMMA = np.pi / 3  # the 60 degree cell of a (111) surface, a = b
DELTA_K = 2 * np.pi / (4.84 * np.sin(GAMMA))
K_PARALLEL = 1.3  # bohr^-1, along x (phi = 0)
CHANNELS = [(n, m) for n in range(-4, 5) for m in range(-4, 5)]


def in_plane_channels(k_x: np.ndarray, k_y: np.ndarray) -> set[tuple[int, int]]:
    """The channels whose K lies along x (the incident direction)."""
    return {c for c, inside in zip(CHANNELS, in_plane(k_x, k_y, np.array([1.0, 0.0])), strict=True) if inside}


def skewed_incidence() -> tuple[float, float]:
    """kx, ky along the skewed reciprocal axes of a parallel momentum K_PARALLEL along x (process.cpp)."""
    kx = K_PARALLEL / np.sin(GAMMA)
    return kx, np.cos(GAMMA) * kx


def test_in_plane_set_on_the_hexagonal_cell() -> None:
    """phi = 0 on the 60 degree cell: K_y = m dk - n dk cos 60 = 0 for m = n / 2, so (0, 0), +-(2, 1), +-(4, 2)."""
    kx, ky = skewed_incidence()
    n, m = np.array(CHANNELS).T
    k_x, k_y = parallel_momentum_cartesian(
        n, m, kx_per_bohr=kx, ky_per_bohr=ky, delta_kx_per_bohr=DELTA_K, delta_ky_per_bohr=DELTA_K, skew_angle_rad=GAMMA
    )
    assert in_plane_channels(k_x, k_y) == {(0, 0), (2, 1), (-2, -1), (4, 2), (-4, -2)}
    assert abs(np.hypot(k_x[CHANNELS.index((0, 0))], k_y[CHANNELS.index((0, 0))]) - K_PARALLEL) < 1e-15


def test_skewed_components_taken_as_cartesian_give_another_set() -> None:
    """Negative test: K = (kx + n dk, ky + m dk), the skewed components used as Cartesian (a realistic slip),
    finds another set of channels along x."""
    kx, ky = skewed_incidence()
    n, m = np.array(CHANNELS).T
    wrong = in_plane_channels(kx + n * DELTA_K, ky + m * DELTA_K)
    assert wrong != {(0, 0), (2, 1), (-2, -1), (4, 2), (-4, -2)}
