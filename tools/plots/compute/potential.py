"""The run's potential (potential/) for plots: energies relative to the gas-phase minimum, the (x, y) grid and the
unit cell in Cartesian coordinates, and the minimax (mountain-pass) barrier on V(Z, r) at one (x, y) point.

Cartesian positions: grid point (ix, iy) lies at x = ix a / nx along the first lattice vector and y = iy b / ny
along the second (src/grids/spatial_grids.cpp: the x grid starts at 0; src/input/process.cpp: dx = a / nx), which
makes the skew angle gamma with the first: R = (x + y cos(gamma), y sin(gamma)). The program's PES uses the inverse
of this map: it splits a bond (cos(phi), sin(phi)) into the skewed components (sin(gamma - phi), sin(phi)) /
sin(gamma) (src/potential/dai_light_cu111_potential.cpp, value). The reciprocal axes of compute/grids.py,
(sin(gamma), -cos(gamma)) and (0, 1), are perpendicular to the lattice vectors (cos(gamma), sin(gamma)) and (1, 0),
as reciprocal axes must be.
"""

import heapq
from dataclasses import dataclass

import numpy as np

from readers.potential import PotentialData, PotentialLimits, Site
from readers.tables import AsymptoticPotential

ENERGY_CAP_EV = 3.0  # heuristic: maps and elbows show V up to this above the gas-phase minimum (walls reach keV)
CONTOUR_STEP_EV = 0.1  # the elbows' contour spacing
V_MAX_CLAMP_RTOL = 1e-12  # relative: a barrier this close to VMax is the clamp, not the surface


def relative_ev(v_hartree: np.ndarray, reference_hartree: float, hartree_to_ev: float) -> np.ndarray:
    """V - V_ref in eV (the reference: the gas-phase minimum, min V_asy)."""
    return (np.asarray(v_hartree) - reference_hartree) * hartree_to_ev


def gas_phase_minimum_hartree(asymptotic: AsymptoticPotential) -> float:
    """min V_asy(r): the reference of the potential plots (hartree)."""
    return float(np.min(asymptotic.v_hartree))


def nearest_index(grid_bohr: np.ndarray, value_bohr: float) -> int:
    """The index of the grid point nearest to a value."""
    return int(np.argmin(np.abs(np.asarray(grid_bohr) - value_bohr)))


def bond_length_index(asymptotic: AsymptoticPotential, r_bohr: np.ndarray) -> int:
    """The r grid index nearest the gas-phase bond length (where V_asy is lowest)."""
    return nearest_index(r_bohr, float(asymptotic.r_bohr[int(np.argmin(asymptotic.v_hartree))]))


def cartesian_bohr(fx: np.ndarray, fy: np.ndarray, lattice: dict) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian (X, Y) in bohr of fractional coordinates (fx, fy) along the lattice vectors (run_info.json's
    lattice: a, b, skew_angle_rad)."""
    x = np.asarray(fx) * lattice["a"]
    y = np.asarray(fy) * lattice["b"]
    gamma = lattice["skew_angle_rad"]
    return x + y * np.cos(gamma), y * np.sin(gamma)


@dataclass
class CellGeometry:
    """The (x, y) grid drawn in Cartesian bohr, periodically: the corners of the cells around the grid points
    ix = 0 .. nx and iy = 0 .. ny (the last row and column repeat the first), and the unit cell's outline."""

    corners_x_bohr: np.ndarray  # (nx + 2, ny + 2): cell (ix, iy) spans fractions (ix -+ 1/2) / nx, (iy -+ 1/2) / ny
    corners_y_bohr: np.ndarray
    outline_x_bohr: np.ndarray  # (5,): the unit cell, closed
    outline_y_bohr: np.ndarray


def cell_geometry(nx: int, ny: int, lattice: dict) -> CellGeometry:
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


def site_positions_bohr(sites: list[Site], lattice: dict) -> list[tuple[str, float, float]]:
    """(name, X, Y) of each site in Cartesian bohr."""
    positions = []
    for site in sites:
        x, y = cartesian_bohr(np.array(site.x), np.array(site.y), lattice)
        positions.append((site.name, float(x), float(y)))
    return positions


@dataclass
class Barrier:
    """The minimax barrier of V(Z, r): the lowest possible highest point of a path over the grid from the gas-phase
    side to the dissociated side, where it lies, and one path that reaches it (not a minimum-energy path: below the
    barrier any route will do)."""

    height_hartree: float  # inf when the dissociated side cannot be reached
    iz: int  # the barrier's point (the path's start when there is no path)
    ir: int
    path: list[tuple[int, int]]  # (iz, ir) from the start to the first point with r >= r_c the search reaches
    at_v_max: bool  # the barrier is the VMax clamp


def minimax_path(v: np.ndarray, start: tuple[int, int], is_product: np.ndarray) -> tuple[float, list[tuple[int, int]]]:
    """The lowest possible highest V of a path over the grid from `start` to a point where is_product holds, and one
    such path; (inf, []) when none is reachable. A bottleneck Dijkstra: points are taken in the order of the highest
    V on the way to them; steps go to the 4 neighbours (no diagonal steps, so the barrier is an upper bound on one
    with them); among equal ways the lower point goes first; a NaN point cannot be crossed."""
    passable = np.where(np.isnan(v), np.inf, v).tolist()  # plain lists: fast scalar access
    n_z, n_r = v.shape
    best = [[np.inf] * n_r for _ in range(n_z)]
    best[start[0]][start[1]] = passable[start[0]][start[1]]
    previous: dict[tuple[int, int], tuple[int, int]] = {}
    heap = [(best[start[0]][start[1]], passable[start[0]][start[1]], start)]
    while heap:
        highest, _, (iz, ir) = heapq.heappop(heap)
        if highest > best[iz][ir]:
            continue
        if is_product[iz, ir]:
            return float(highest), path_to((iz, ir), previous)
        for jz, jr in ((iz + 1, ir), (iz - 1, ir), (iz, ir + 1), (iz, ir - 1)):
            if not (0 <= jz < n_z and 0 <= jr < n_r):
                continue
            through = max(highest, passable[jz][jr])
            if through < best[jz][jr]:
                best[jz][jr] = through
                previous[(jz, jr)] = (iz, ir)
                heapq.heappush(heap, (through, passable[jz][jr], (jz, jr)))
    return float("inf"), []


def path_to(cell: tuple[int, int], previous: dict[tuple[int, int], tuple[int, int]]) -> list[tuple[int, int]]:
    """The cells from the search's start to `cell`, following `previous` back."""
    path = [cell]
    while path[-1] in previous:
        path.append(previous[path[-1]])
    return path[::-1]


def barrier_on(v_zr_hartree: np.ndarray, ir_c: int, limits: PotentialLimits) -> Barrier:
    """The minimax barrier of V(Z, r) (nz_potential, nr): from the lowest point of the last Z row (the gas-phase
    side, just below zasymp) to any point with r index >= ir_c (the dissociated side: the r absorber's first point,
    the program's own index). An all-NaN last row has no start: an infinite barrier at its first point."""
    v = np.asarray(v_zr_hartree, dtype=float)
    last = v.shape[0] - 1
    if np.all(np.isnan(v[last])):
        return Barrier(float("inf"), last, 0, [], False)
    start = (last, int(np.nanargmin(v[last])))
    is_product = np.zeros(v.shape, dtype=bool)
    is_product[:, ir_c:] = True
    height, path = minimax_path(v, start, is_product)
    if not path:
        return Barrier(height, start[0], start[1], [], False)
    heights = [v[cell] for cell in path]
    iz, ir = path[int(np.argmax(heights))]
    at_v_max = limits.impose_v_max and height >= limits.v_max_clamp_hartree * (1 - V_MAX_CLAMP_RTOL)
    return Barrier(height, iz, ir, path, bool(at_v_max))


@dataclass
class ElbowPoint:
    """An (x, y) point whose elbow is drawn: a site (isite its index) or, without sites, the lowest barrier's point
    (isite None)."""

    name: str
    ix: int
    iy: int
    barrier: Barrier
    isite: int | None


def lowest_barrier(potential: PotentialData, ir_c: int) -> ElbowPoint:
    """The (x, y) grid point whose V_min(Z, r) has the lowest minimax barrier (the first in (ix, iy) order), named
    by its fractional coordinates."""
    lowest = None
    for ix in range(potential.nx):
        for iy in range(potential.ny):
            barrier = barrier_on(potential.v_min_hartree[:, :, ix, iy], ir_c, potential.limits)
            if lowest is None or barrier.height_hartree < lowest.barrier.height_hartree:
                name = f"lowest barrier (x {ix / potential.nx:g}, y {iy / potential.ny:g})"
                lowest = ElbowPoint(name, ix, iy, barrier, None)
    return lowest


def elbow_points(potential: PotentialData, ir_c: int) -> list[ElbowPoint]:
    """The sites with their barriers on V_min(Z, r); without sites, the grid point with the lowest barrier."""
    if not potential.sites:
        return [lowest_barrier(potential, ir_c)]
    points = []
    for k, site in enumerate(potential.sites):
        barrier = barrier_on(potential.v_min_hartree[:, :, site.ix, site.iy], ir_c, potential.limits)
        points.append(ElbowPoint(site.name, site.ix, site.iy, barrier, k))
    return points


def lowest_site(potential: PotentialData, ir_c: int) -> ElbowPoint:
    """The site with the lowest barrier on V_min(Z, r) (the first of equal ones)."""
    points = elbow_points(potential, ir_c)
    return points[int(np.argmin([p.barrier.height_hartree for p in points]))]


def map_point(
    potential: PotentialData, lowest: ElbowPoint, z_bohr: float | None, r_bohr: float | None
) -> tuple[int, int]:
    """(iz, ir) of the left potential map: the grid points nearest to the asked Z and r, by default the lowest
    barrier's."""
    iz = lowest.barrier.iz if z_bohr is None else nearest_index(potential.z_bohr, z_bohr)
    ir = lowest.barrier.ir if r_bohr is None else nearest_index(potential.r_bohr, r_bohr)
    return iz, ir


def orientation_at(potential: PotentialData, point: ElbowPoint, iz: int, ir: int) -> tuple[int, int]:
    """(itheta, iphi) of the lowest V over the orientations at an (x, y) point's (Z, r) point (v_min_node)."""
    node = potential.v_min_node[iz, ir, point.ix, point.iy]
    return int(node[0]), int(node[1])


def point_as_site(potential: PotentialData, point: ElbowPoint) -> Site:
    """An elbow point as a site to mark on the maps (its fractional coordinates from its grid point)."""
    return Site(point.name, point.ix / potential.nx, point.iy / potential.ny, point.ix, point.iy)
