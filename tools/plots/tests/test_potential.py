"""compute/potential.py on synthetic surfaces (the minimax pass and its kinds, the VMax flag, the choice of points on a
grid with nx != ny, where an (ix, iy) swap shows, the nodes' cells), compute/grids.py's Cartesian positions and
periodic cells, and the potential colour bar."""

import numpy as np

from compute.grids import cartesian_bohr, cell_geometry, periodic
from compute.potential import (
    AT_R_C,
    AT_START,
    SADDLE,
    barrier_on,
    barriers_on_grid,
    lowest_barrier,
    map_point,
    minimax_path,
    node_edges_deg,
    orientation_at,
    site_points,
)
from panels.potential import add_potential_colorbar, potential_color_scale, potential_colormap
from readers.potential import PotentialData, PotentialLimits
from readers.run_info import Lattice
from readers.sites import Site
from style import DIVERGING, SEQUENTIAL, new_figure

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
    highest point is a saddle (inside the path); the VMax flag only when the barrier equals VMax (to 1e-12
    relative)."""
    v = ridge_surface()
    v[4, 0] = 0.5
    barrier = barrier_on(v, 3, NO_LIMITS)
    assert barrier.path[0] == (4, 1) and barrier.path[-1][1] == 3
    assert barrier.height_hartree == 1.5 and (barrier.iz, barrier.ir) == (0, 2) and barrier.saddle
    assert not barrier.at_v_max
    clamp = PotentialLimits(False, 6.0, 8.0, False, -1.0, True, 1.5)
    assert barrier_on(v, 3, clamp).at_v_max
    above = PotentialLimits(False, 6.0, 8.0, False, -1.0, True, 1.5 * (1 + 2e-12))
    assert not barrier_on(v, 3, above).at_v_max


def test_kinds_of_the_highest_point() -> None:
    """V rising with r all the way to r_c: the highest point is the path's end at r_c (no saddle); V falling from
    the gas-phase side: the highest point is the start (no barrier); a hump before r_c: a saddle. Only the last is
    a barrier."""
    rising = np.tile(np.arange(5.0), (5, 1))  # V = r index everywhere
    barrier = barrier_on(rising, 3, NO_LIMITS)
    assert (barrier.kind, barrier.height_hartree, barrier.ir) == (AT_R_C, 3.0, 3) and not barrier.saddle
    falling = -rising
    barrier = barrier_on(falling, 3, NO_LIMITS)
    assert barrier.kind == AT_START and barrier.path[0] == (barrier.iz, barrier.ir) and not barrier.saddle
    hump = rising.copy()
    hump[:, 2] = 4.0
    assert barrier_on(hump, 3, NO_LIMITS).kind == SADDLE


def test_cartesian_positions_and_periodic_cells() -> None:
    """On a hexagonal cell (a = b = 2, gamma = 60 deg) the fraction (0, 0.5) lies at (0.5, sqrt(3)/2) (treating the
    skewed coordinates as Cartesian would give (0, 1)); the cells' corners lie half a step around the points, one
    more row and column cover the cell, whose values repeat the first; the outline closes."""
    lattice = Lattice(2.0, 2.0, np.pi / 3)
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
    return PotentialData(v_min, node, None, z, r, 3, 2, np.zeros(3), np.zeros(2), True, sites, NO_LIMITS)


def test_points_are_chosen_on_a_grid_with_nx_not_ny() -> None:
    """The lowest barrier lies at (ix, iy) = (2, 1) (a swap would look for (1, 2), outside the 3 x 2 grid), with its
    fractions; a site's barrier is its own; the map's point follows --z, --r or the barrier; the orientation is read
    at (ix, iy), not (iy, ix)."""
    potential = synthetic_potential([])
    lowest = lowest_barrier(potential, barriers_on_grid(potential, 3))
    assert (lowest.ix, lowest.iy) == (2, 1) and lowest.barrier.height_hartree == 0.5 and lowest.isite is None
    assert (lowest.fraction_x, lowest.fraction_y) == (2 / 3, 0.5) and lowest.name is None
    with_site = synthetic_potential([Site("corner", 2 / 3, 0.0, 2, 0)])
    (point,) = site_points(with_site, 3)
    assert (point.ix, point.iy, point.isite, point.name) == (2, 0, 0, "corner") and point.barrier.height_hartree == 1.5
    assert map_point(potential, lowest, None, None) == (0, 2) and map_point(potential, lowest, 3.2, 0.9) == (3, 1)
    assert orientation_at(potential, lowest, 0, 2) == (2, 1)


def test_node_edges() -> None:
    """The theta cells reach 0 and 90 deg (homonuclear) or 180, midway between nodes; the phi cells half a step
    around each node."""
    theta, phi = np.radians([39.2, 90.0]), np.radians(np.arange(0.0, 360.0, 60.0))
    theta_edges, phi_edges = node_edges_deg(theta, phi, homonuclear=True)
    assert np.allclose(theta_edges, [0.0, 64.6, 90.0]) and np.allclose(phi_edges, np.arange(-30.0, 360.0, 60.0))
    assert node_edges_deg(theta, phi, homonuclear=False)[0][-1] == 180.0


def test_colour_scale_is_sequential_unless_a_deep_well() -> None:
    """Barriers only (0.2 to 1.5 eV, walls far above the cap): sequential from zero to the largest value below the cap
    (3 eV), the colour bar closed (nothing below its bottom). A 10 meV dip (below 5 % of the top) keeps it sequential
    from the dip. A 0.5 eV well: diverging, -3 to 3 eV, the bar open at the bottom. NaN is ignored; only NaN gives
    the cap."""
    barriers = np.array([0.2, 1.5, np.nan])
    scale = potential_color_scale(3.0, [barriers, np.array([1.0])])
    assert (scale.vmin, scale.vmax) == (0.0, 1.5) and potential_colormap(scale) == SEQUENTIAL
    assert potential_color_scale(3.0, [np.array([0.2, 2e3])]).vmax == 3.0
    dip = potential_color_scale(3.0, [np.array([-0.01, 0.2, 2e3])])
    assert (dip.vmin, dip.vmax) == (-0.01, 3.0) and potential_colormap(dip) == SEQUENTIAL
    well = potential_color_scale(3.0, [np.array([-0.5, 0.2, 2e3])])
    assert (well.vmin, well.vmax) == (-3.0, 3.0) and potential_colormap(well) == DIVERGING
    assert (potential_color_scale(3.0, [np.array([np.nan])]).vmax) == 3.0
    fig = new_figure()
    assert add_potential_colorbar(fig, [fig.add_subplot()], scale).extend == "neither"
    assert add_potential_colorbar(fig, [fig.add_subplot()], well).extend == "min"
