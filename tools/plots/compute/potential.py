"""The run's potential (potential/) for plots: energies relative to the gas-phase minimum, the minimax
(mountain-pass) barrier on V(Z, r) at one (x, y) point and at every one, the points the figures show (and their
Cartesian positions), and the cell edges of the DVR nodes. The Cartesian map and the cell geometry are in
compute/grids.py.

V_min, the minimum over the DVR nodes at fixed (x, y), is the lower envelope over the run's orientations: each point
has its own orientation, its contours kink where the node changes, and it holds no rotational or zero-point energy.
"""

import heapq
from dataclasses import dataclass

import numpy as np

from compute.grids import cartesian_bohr
from readers.potential import PotentialData, PotentialLimits
from readers.run_info import Lattice
from readers.sites import Site
from readers.tables import AsymptoticPotential

ENERGY_CAP_EV = 3.0  # heuristic: maps and elbows show V up to this above the gas-phase minimum (walls reach keV)
CONTOUR_STEP_EV = 0.1  # the elbows' filled-contour spacing
FRAME_CONTOUR_STEP_EV = 0.5  # the spacing of V_min's contour lines over the frame movies (a divisor of the cap)
CONTOUR_LINE_EVERY = 5  # a contour line every this many filled steps
MAX_ELBOW_SITES = 4  # the elbows of the first this many sites of PlotSites
V_MAX_CLAMP_RTOL = 1e-12  # relative: a barrier this close to VMax is the clamp, not the surface
# A write2DPES cut: the PES itself, no switching and no clamp.
RAW_PES = PotentialLimits(
    use_switching=False,
    z_min_switch_bohr=0.0,
    z_max_switch_bohr=0.0,
    impose_v_min=False,
    v_min_clamp_hartree=0.0,
    impose_v_max=False,
    v_max_clamp_hartree=0.0,
)


def relative_ev(v_hartree: np.ndarray, reference_hartree: float, hartree_to_ev: float) -> np.ndarray:
    """V - V_ref in eV (the reference: the gas-phase minimum, min V_asy)."""
    return (np.asarray(v_hartree) - reference_hartree) * hartree_to_ev


def gas_phase_minimum_hartree(asymptotic: AsymptoticPotential) -> float:
    """The reference of the potential plots (hartree): the lowest value on the r grid of asymptotic_potential.dat,
    V(r) at Z = zasymp (with switching: ZmaxSwitch), theta = pi/2, at the program's x = y = 0."""
    return float(np.min(asymptotic.v_hartree))


def nearest_index(grid_bohr: np.ndarray, value_bohr: float) -> int:
    """The index of the grid point nearest to a value."""
    return int(np.argmin(np.abs(np.asarray(grid_bohr) - value_bohr)))


def bond_length_index(asymptotic: AsymptoticPotential, r_bohr: np.ndarray) -> int:
    """The r grid index nearest the gas-phase bond length (where V_asy is lowest)."""
    return nearest_index(r_bohr, float(asymptotic.r_bohr[int(np.argmin(asymptotic.v_hartree))]))


SADDLE, AT_R_C, AT_START, UNREACHED = "saddle", "at r_c", "at the start", "unreached"  # where the highest point is


@dataclass
class Barrier:
    """The minimax pass of V(Z, r): the lowest possible highest point of a path over the grid from the gas-phase
    side to the dissociated side, where it lies, and one path that reaches it (not a minimum-energy path: below it
    any route will do). Only when that highest point lies inside the path is it a barrier (a saddle); it can also be
    the path's end at r_c (no saddle before the dissociated side; a lower hump may lie before it) or its start (V
    falls from the gas-phase side: no barrier)."""

    height_hartree: float  # inf when the dissociated side cannot be reached
    iz: int  # the highest point (the path's start when there is no path)
    ir: int
    path: list[tuple[int, int]]  # (iz, ir) from the start to the first point with r >= r_c the search reaches
    kind: str  # SADDLE, AT_R_C, AT_START or UNREACHED
    saddle: bool  # kind == SADDLE: the highest point is a barrier (inside the path)
    at_v_max: bool  # the highest point is the VMax clamp


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
    """The minimax barrier of V(Z, r) (nz, nr: the run's potential below zasymp, or a write2DPES cut over the whole
    Z grid): from the lowest point of the last Z row (the gas-phase side) to any point with r index >= ir_c (the
    dissociated side: the r absorber's first point, the program's own index). An all-NaN last row has no start: an
    infinite barrier at its first point."""
    v = np.asarray(v_zr_hartree, dtype=float)
    last = v.shape[0] - 1
    if np.all(np.isnan(v[last])):
        return Barrier(float("inf"), last, 0, [], UNREACHED, False, False)
    start = (last, int(np.nanargmin(v[last])))
    is_product = np.zeros(v.shape, dtype=bool)
    is_product[:, ir_c:] = True
    height, path = minimax_path(v, start, is_product)
    if not path:
        return Barrier(height, start[0], start[1], [], UNREACHED, False, False)
    heights = [v[cell] for cell in path]
    highest = int(np.argmax(heights))  # the first of equal ones
    iz, ir = path[highest]
    kind = AT_START if highest == 0 else AT_R_C if highest == len(path) - 1 else SADDLE  # only the end has r >= r_c
    at_v_max = limits.impose_v_max and height >= limits.v_max_clamp_hartree * (1 - V_MAX_CLAMP_RTOL)
    return Barrier(height, iz, ir, path, kind, kind == SADDLE, bool(at_v_max))


@dataclass
class PointOfInterest:
    """An (x, y) grid point a figure shows: a site of PlotSites (name and isite set) or the grid point with the
    lowest barrier (name and isite None), with its fractions and its barrier on V_min(Z, r)."""

    name: str | None
    ix: int
    iy: int
    fraction_x: float
    fraction_y: float
    barrier: Barrier
    isite: int | None


def barriers_on_grid(potential: PotentialData, ir_c: int) -> list[list[Barrier]]:
    """The minimax barrier on V_min(Z, r) at every (x, y) grid point: [ix][iy]."""
    return [
        [barrier_on(potential.v_min_hartree[:, :, ix, iy], ir_c, potential.limits) for iy in range(potential.ny)]
        for ix in range(potential.nx)
    ]


def lowest_barrier(potential: PotentialData, barriers: list[list[Barrier]]) -> PointOfInterest:
    """The (x, y) grid point of barriers_on_grid with the lowest saddle (the first in (ix, iy) order); without any
    saddle, the one with the lowest pass."""
    ranks = [(not b.saddle, b.height_hartree, ix, iy) for ix, row in enumerate(barriers) for iy, b in enumerate(row)]
    _, _, ix, iy = min(ranks)
    return PointOfInterest(None, ix, iy, ix / potential.nx, iy / potential.ny, barriers[ix][iy], None)


def site_points(potential: PotentialData, ir_c: int) -> list[PointOfInterest]:
    """The sites of PlotSites with their barriers on V_min(Z, r)."""
    points = []
    for k, site in enumerate(potential.sites):
        barrier = barrier_on(potential.v_min_hartree[:, :, site.ix, site.iy], ir_c, potential.limits)
        points.append(PointOfInterest(site.name, site.ix, site.iy, site.fraction_x, site.fraction_y, barrier, k))
    return points


def lowest_site(points: list[PointOfInterest]) -> PointOfInterest:
    """The site with the lowest saddle (the first of equal ones); without any saddle, the one with the lowest pass."""
    return min(points, key=lambda p: (not p.barrier.saddle, p.barrier.height_hartree))


def positions_bohr(
    points: list[Site | PointOfInterest], lattice: Lattice
) -> list[tuple[Site | PointOfInterest, float, float]]:
    """Each site or point with its Cartesian (X, Y) in bohr."""
    positions = []
    for point in points:
        x, y = cartesian_bohr(np.array(point.fraction_x), np.array(point.fraction_y), lattice)
        positions.append((point, float(x), float(y)))
    return positions


def map_point(
    potential: PotentialData, lowest: PointOfInterest, z_bohr: float | None, r_bohr: float | None
) -> tuple[int, int]:
    """(iz, ir) of the V_min map: the grid points nearest to the asked Z and r, by default the lowest barrier's."""
    iz = lowest.barrier.iz if z_bohr is None else nearest_index(potential.z_bohr, z_bohr)
    ir = lowest.barrier.ir if r_bohr is None else nearest_index(potential.r_bohr, r_bohr)
    return iz, ir


def orientation_at(potential: PotentialData, point: PointOfInterest, iz: int, ir: int) -> tuple[int, int]:
    """(itheta, iphi) of the lowest V over the orientations at a point's (Z, r) (v_min_node; the first of equal
    nodes)."""
    node = potential.v_min_node[iz, ir, point.ix, point.iy]
    return int(node[0]), int(node[1])


def node_edges_deg(theta_rad: np.ndarray, phi_rad: np.ndarray, homonuclear: bool) -> tuple[np.ndarray, np.ndarray]:
    """The edges (deg) of the cells around the DVR nodes: theta from 0 to 90 (homonuclear) or 180, midway between
    nodes; phi half a step around each node (the nodes are equally spaced)."""
    theta = np.degrees(theta_rad)
    theta_edges = np.concatenate([[0.0], (theta[1:] + theta[:-1]) / 2, [90.0 if homonuclear else 180.0]])
    phi = np.degrees(phi_rad)
    step = 360.0 / len(phi)
    return theta_edges, np.concatenate([phi - step / 2, [phi[-1] + step / 2]])


def count_not_finite(values: np.ndarray) -> int:
    """How many values are NaN or inf (drawn hatched, like the walls above the cap)."""
    return int(np.count_nonzero(~np.isfinite(np.asarray(values))))


def saddle_heights_hartree(barriers: list[list[Barrier]]) -> np.ndarray:
    """(nx, ny): the barrier heights of barriers_on_grid, NaN where the pass has no saddle (or r_c is not reached)."""
    return np.array([[b.height_hartree if b.saddle else np.nan for b in row] for row in barriers])
