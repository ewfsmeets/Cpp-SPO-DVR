"""The viewer's pure helpers on hand-made grids: the grid point whose drawn cell holds a click on the skewed cell
(compute.grids.grid_point_of_cell_at), its inverse (grid_point_bohr), the index window of a range in bohr
(compute.frames.index_window) and the frame nearest a time (compute.frames.nearest_frame)."""

import numpy as np

from compute.frames import index_window, nearest_frame
from compute.grids import cartesian_bohr, fractions_of, grid_point_bohr, grid_point_of_cell_at, site_images
from copies.point_group import GENERATORS, apply_to_cartesian, group_closure, rotation
from readers.run_info import Lattice

HEXAGONAL = Lattice(4.84, 4.84, np.radians(60.0))  # the Cu(111) cell of the Dai-Light fixtures
OBLONG = Lattice(5.0, 3.0, np.radians(70.0))  # a != b: a swapped a and b shows


def test_a_click_in_each_cell_picks_its_grid_point() -> None:
    """On a 6 x 4 grid of the 60-degree cell and of a 70-degree cell with a != b: every grid point's Cartesian
    position (grid_point_bohr), and a point a third of a cell step off it, give that point; a point just inside the
    cell's far corner wraps to (0, 0). Reading the Cartesian position as if the axes were orthogonal (a slip: fx =
    X / a, fy = Y / b) picks another point for 18 of the 24 points of the first, 12 of the second."""
    nx, ny = 6, 4
    for lattice, misses in ((HEXAGONAL, 18), (OBLONG, 12)):
        slipped = 0
        for ix in range(nx):
            for iy in range(ny):
                x, y = grid_point_bohr(ix, iy, nx, ny, lattice)
                assert grid_point_of_cell_at(x, y, nx, ny, lattice) == (ix, iy)
                x_off, y_off = cartesian_bohr(np.array((ix + 0.33) / nx), np.array((iy - 0.33) / ny), lattice)
                assert grid_point_of_cell_at(float(x_off), float(y_off), nx, ny, lattice) == (ix, iy)
                wrong = (round(x / lattice.a_bohr * nx) % nx, round(y / lattice.b_bohr * ny) % ny)
                slipped += wrong != (ix, iy)
        assert slipped == misses
        corner_x, corner_y = cartesian_bohr(np.array(1.0 - 0.01), np.array(1.0 - 0.01), lattice)
        assert grid_point_of_cell_at(float(corner_x), float(corner_y), nx, ny, lattice) == (0, 0)


def test_index_window() -> None:
    """The points nearest the two ends, both included, in either order; ends that meet give one point."""
    points = np.arange(6) * 0.5
    assert index_window(points, 0.4, 1.6) == slice(1, 4)
    assert index_window(points, 1.6, 0.4) == slice(1, 4)
    assert index_window(points, 2.5, 2.5) == slice(5, 6)


def test_nearest_frame() -> None:
    """The record nearest a time, the first of two equally near."""
    times = np.array([0.0, 100.0, 200.0])
    assert nearest_frame(times, 140.0) == 1
    assert nearest_frame(times, 150.0) == 1
    assert nearest_frame(times, 1e9) == 2


def image_fractions(fraction_x: float, fraction_y: float, lattice: Lattice, group: str) -> set[tuple[float, float]]:
    """site_images as fractions rounded to 6 digits."""
    return {
        tuple(round(v, 6) for v in fractions_of(x, y, lattice))
        for x, y in site_images(fraction_x, fraction_y, lattice, group)
    }


def test_site_images_by_hand() -> None:
    """The points of the cell (its outline included) equivalent to a site: on the 60-degree cell (C6v) top at the
    four corners, a bridge at the three edge midpoints and the long diagonal's middle (five with the copies on the
    far edges), the two hollows; on a square cell (C4v) the point (0.25, 0) at its four images, all on edges, eight with their copies; C3v (p3m1:
    a mirror at 90 degrees) splits the hollows; "none" keeps the site alone. The site's own position comes first.
    Rotations alone for C6v (a slip: the mirrors left out) give (0.2, 0.1) 6 images instead of 12."""
    third, two_thirds = round(1 / 3, 6), round(2 / 3, 6)
    assert image_fractions(0.0, 0.0, HEXAGONAL, "C6v") == {(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)}
    bridges = {(0.5, 0.0), (0.5, 1.0), (0.0, 0.5), (1.0, 0.5), (0.5, 0.5)}
    assert image_fractions(0.5, 0.0, HEXAGONAL, "C6v") == bridges
    assert image_fractions(1 / 3, 1 / 3, HEXAGONAL, "C6v") == {(third, third), (two_thirds, two_thirds)}
    assert image_fractions(1 / 3, 1 / 3, HEXAGONAL, "C3v") == {(third, third)}
    square = Lattice(4.0, 4.0, np.radians(90.0))
    on_edges = {(0.25, 0.0), (0.75, 0.0), (0.0, 0.25), (0.0, 0.75)}  # each also on the opposite edge
    assert image_fractions(0.25, 0.0, square, "C4v") == on_edges | {(0.25, 1.0), (0.75, 1.0), (1.0, 0.25), (1.0, 0.75)}
    assert image_fractions(0.25, 0.1, HEXAGONAL, "none") == {(0.25, 0.1)}
    first = site_images(0.5, 0.0, HEXAGONAL, "C6v")[0]
    assert tuple(round(v, 6) for v in fractions_of(*first, HEXAGONAL)) == (0.5, 0.0)
    assert len(site_images(0.2, 0.1, HEXAGONAL, "C6v")) == 12
    rotations_only = {
        tuple(
            round(v, 6)
            for v in fractions_of(
                *apply_to_cartesian(g, *cartesian_bohr(np.array(0.2), np.array(0.1), HEXAGONAL)), HEXAGONAL
            )
        )
        for g in group_closure([rotation(2)])
    }
    assert len(rotations_only) == 6 and set(GENERATORS) == {"C2v", "C4v", "C6v", "C3v", "C3", "none"}
