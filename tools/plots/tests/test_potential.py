"""compute/potential.py on synthetic surfaces: the minimax barrier, its path, the VMax flag, Cartesian positions on
a skewed cell, the periodic cells, and the choice of points on a grid with nx != ny (an (ix, iy) swap shows)."""

import numpy as np

from compute.potential import (
    barrier_on,
    cartesian_bohr,
    cell_geometry,
    elbow_points,
    lowest_barrier,
    map_point,
    minimax_path,
    orientation_at,
    periodic,
)
from readers.potential import PotentialData, PotentialLimits, Site

NO_LIMITS = PotentialLimits(False, 6.0, 8.0, False, -1.0, False, 10.0)


def ridge_surface() -> np.ndarray:
    """A 5 x 5 V(Z, r): 0 everywhere but a ridge at r index 2 of height 2.0 with a gap of 1.5 at Z index 0."""
    v = np.zeros((5, 5))
    v[:, 2] = 2.0
    v[0, 2] = 1.5
    return v


def test_minimax_goes_through_the_lowest_gap() -> None:
    """From (4, 0) to r index >= 3 the lowest highest point is the gap (1.5), not the straight crossing at the start's
    Z (2.0); the path starts at the start, crosses at the gap and ends on the far side."""
    v = ridge_surface()
    is_product = np.zeros(v.shape, dtype=bool)
    is_product[:, 3:] = True
    height, path = minimax_path(v, (4, 0), is_product)
    assert height == 1.5
    assert path[0] == (4, 0) and (0, 2) in path and path[-1][1] == 3
    straight = max(v[4, 0], v[4, 1], v[4, 2], v[4, 3])  # a slip: the straight line at the start's Z
    assert straight == 2.0


def test_unreachable_and_nan_points() -> None:
    """A NaN gap cannot be crossed (the barrier becomes the ridge); a side that cannot be reached gives inf; an
    all-NaN last row gives an infinite barrier without a path."""
    v = ridge_surface()
    v[0, 2] = np.nan
    is_product = np.zeros(v.shape, dtype=bool)
    is_product[:, 3:] = True
    assert minimax_path(v, (4, 0), is_product)[0] == 2.0
    v[:, 2] = np.nan
    assert minimax_path(v, (4, 0), is_product) == (float("inf"), [])
    v[4, :] = np.nan
    barrier = barrier_on(v, 3, NO_LIMITS)
    assert barrier.height_hartree == float("inf") and barrier.path == []


def test_barrier_starts_at_the_lowest_point_of_the_last_row() -> None:
    """barrier_on starts at the lowest V of the last Z row and stops at the first point with r index >= ir_c; its
    highest point; the VMax flag only when the barrier equals VMax (to 1e-12 relative)."""
    v = ridge_surface()
    v[4, 0] = 0.5
    barrier = barrier_on(v, 3, NO_LIMITS)
    assert barrier.path[0] == (4, 1) and barrier.path[-1][1] == 3
    assert barrier.height_hartree == 1.5 and (barrier.iz, barrier.ir) == (0, 2) and not barrier.at_v_max
    clamp = PotentialLimits(False, 6.0, 8.0, False, -1.0, True, 1.5)
    assert barrier_on(v, 3, clamp).at_v_max
    above = PotentialLimits(False, 6.0, 8.0, False, -1.0, True, 1.5 * (1 + 2e-12))
    assert not barrier_on(v, 3, above).at_v_max


def test_cartesian_positions_and_periodic_cells() -> None:
    """On a hexagonal cell (a = b = 2, gamma = 60 deg) the fraction (0, 0.5) lies at (0.5, sqrt(3)/2) (treating the
    skewed coordinates as Cartesian would give (0, 1)); the cells' corners lie half a step around the points, one
    more row and column cover the cell, whose values repeat the first; the outline closes."""
    lattice = {"a": 2.0, "b": 2.0, "skew_angle_rad": np.pi / 3}
    x, y = cartesian_bohr(np.array([0.5, 0.0]), np.array([0.0, 0.5]), lattice)
    assert np.allclose(x, [1.0, 0.5], rtol=0, atol=1e-15) and np.allclose(y, [0.0, np.sqrt(3) / 2], rtol=0, atol=1e-15)
    cell = cell_geometry(2, 4, lattice)
    assert cell.corners_x_bohr.shape == (4, 6)
    corner_x, corner_y = cartesian_bohr(np.array(-0.25), np.array(-0.125), lattice)
    assert np.isclose(cell.corners_x_bohr[0, 0], corner_x) and np.isclose(cell.corners_y_bohr[0, 0], corner_y)
    assert periodic(np.arange(8.0).reshape(2, 4)).tolist() == [[0, 1, 2, 3, 0], [4, 5, 6, 7, 4], [0, 1, 2, 3, 0]]
    assert (cell.outline_x_bohr[0], cell.outline_y_bohr[0]) == (cell.outline_x_bohr[-1], cell.outline_y_bohr[-1])


def synthetic_potential(sites: list[Site]) -> PotentialData:
    """A 3 x 2 (x, y) grid (nx != ny) of the ridge surface with a gap that is lower at (ix, iy) = (2, 1) (0.5) than
    anywhere else (1.5); the orientation node of the lowest V is (itheta, iphi) = (ix, iy) at every point."""
    v_min = np.repeat(ridge_surface()[:, :, None, None], 3, axis=2).repeat(2, axis=3)
    v_min[0, 2, 2, 1] = 0.5
    node = np.zeros((*v_min.shape, 2), dtype=np.int32)
    node[..., 0] = np.arange(3)[None, None, :, None]
    node[..., 1] = np.arange(2)[None, None, None, :]
    z, r = np.arange(5.0), np.arange(5.0)
    return PotentialData(v_min, node, None, z, r, 3, 2, np.zeros(3), np.zeros(2), sites, NO_LIMITS, 27.2114)


def test_points_are_chosen_on_a_grid_with_nx_not_ny() -> None:
    """The lowest barrier lies at (ix, iy) = (2, 1) (a swap would look for (1, 2), outside the 3 x 2 grid); without
    sites the elbow is that point, with a site its own; the maps' point follows --z, --r or the barrier; the
    orientation is read at (ix, iy), not (iy, ix)."""
    potential = synthetic_potential([])
    lowest = lowest_barrier(potential, 3)
    assert (lowest.ix, lowest.iy) == (2, 1) and lowest.barrier.height_hartree == 0.5 and lowest.isite is None
    assert "x 0.666667, y 0.5" in lowest.name
    assert [(p.ix, p.iy) for p in elbow_points(potential, 3)] == [(2, 1)]
    with_site = synthetic_potential([Site("corner", 2 / 3, 0.0, 2, 0)])
    (point,) = elbow_points(with_site, 3)
    assert (point.ix, point.iy, point.isite) == (2, 0, 0) and point.barrier.height_hartree == 1.5
    assert map_point(potential, lowest, None, None) == (0, 2) and map_point(potential, lowest, 3.2, 0.9) == (3, 1)
    assert orientation_at(potential, lowest, 0, 2) == (2, 1)
