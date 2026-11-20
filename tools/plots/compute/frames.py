"""Quantities of frames/ (readers/frames.py): the density over (Z, r) at one (x, y) grid point with V_min's contours
there, the probability over the (x, y) cell in a window of (Z, r), the populations of the rotational states and of
J, and how much of the wave packet sits at the edges of the basis (the top J and |mJ|, the first and last r points,
the first Z points).

Every probability here is per grid cell; compute.grids.density_per_volume makes a density of it. The functions take
the records before any NaN or inf (usable_frames).
"""

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from compute.potential import ENERGY_CAP_EV, FRAME_CONTOUR_STEP_EV, relative_ev
from compute.series import usable_records
from readers.frames import Frames, read_frames
from readers.potential import PotentialData

EDGE_POINTS = 2  # heuristic: the first and last r points and the first Z points counted as the edge of the grid
MAX_FRAME_SITES = 4  # the frame movies show the first this many sites of PlotSites


def usable_frames(runs: list[Path]) -> tuple[Frames, tuple[str, ...]]:
    """frames/ of a run's segments with its records before any NaN or inf, and a note when some were left out
    (compute.series.usable_records; ReadError when none is left)."""
    frames = read_frames(runs)
    series, notes = usable_records(frames.series)
    return replace(frames, series=series), notes


def zr_probability_at(frames: Frames, ix: int, iy: int) -> np.ndarray:
    """(records, nz, nr): the probability per (Z, r, x, y) cell at the one (x, y) grid point (ix, iy)."""
    return np.asarray(frames.series.arrays["prob_zrxy"][:, :, :, ix, iy])


def xy_probability(frames: Frames, record: int, z_slice: slice, r_slice: slice) -> np.ndarray:
    """(nx, ny): the probability per (x, y) cell of one record, summed over the (Z, r) cells of the slices (Z
    indices of the scattering grid)."""
    return np.asarray(frames.series.arrays["prob_zrxy"][record, z_slice, r_slice]).sum(axis=(0, 1))


def xy_probability_records(frames: Frames, z_slice: slice, r_slice: slice) -> np.ndarray:
    """(records, nx, ny): xy_probability of every record (a full pass over the window's cells)."""
    return np.asarray(frames.series.arrays["prob_zrxy"][:, z_slice, r_slice]).sum(axis=(1, 2))


def index_window(points_bohr: np.ndarray, low_bohr: float, high_bohr: float) -> slice:
    """The index slice of the grid points from the one nearest low_bohr to the one nearest high_bohr, both
    included, at least one point (a window whose ends meet holds that point)."""
    low = int(np.argmin(np.abs(points_bohr - low_bohr)))
    high = int(np.argmin(np.abs(points_bohr - high_bohr)))
    return slice(min(low, high), max(low, high) + 1)


def window_probability(frames: Frames, z_slice: slice, r_slice: slice) -> np.ndarray:
    """(records,): the probability in a (Z, r) window of the scattering grid, summed over (x, y), per record."""
    return np.asarray(frames.series.arrays["prob_zrxy"][:, z_slice, r_slice]).sum(axis=(1, 2, 3, 4))


def nearest_frame(times_au: np.ndarray, t_au: float) -> int:
    """The record whose time is nearest t_au (the first of two equally near)."""
    return int(np.argmin(np.abs(np.asarray(times_au) - t_au)))


def most_at_surface(frames: Frames, nz_surface: int) -> int:
    """The record with the most probability at the first nz_surface Z points (those below zasymp)."""
    return int(np.argmax(window_probability(frames, slice(0, nz_surface), slice(None))))


def state_populations_z(frames: Frames, record: int) -> np.ndarray:
    """(nz, nrot): the probability per (Z, rotational state) of one record on the scattering grid, in the basis's
    order."""
    return np.asarray(frames.series.arrays["pop_jmj_z"][record])


def j_populations(frames: Frames) -> tuple[np.ndarray, np.ndarray]:
    """The J values of the basis and, per record, the probability of each J on the grids: pop_jmj_z summed over mJ
    and Z, plus chi beyond the scattering grid (prob_tail_zr) in the initial state's J. Returns (J values,
    (records, n_J))."""
    pop = np.asarray(frames.series.arrays["pop_jmj_z"]).sum(axis=1)  # (records, nrot)
    tail = np.asarray(frames.series.arrays["prob_tail_zr"]).sum(axis=(1, 2))
    j_values = np.unique(frames.states[:, 0])
    per_j = np.stack([pop[:, frames.states[:, 0] == j].sum(axis=1) for j in j_values], axis=1)
    per_j[:, np.searchsorted(j_values, frames.states[frames.initial_index, 0])] += tail
    return j_values, per_j


@dataclass(frozen=True, kw_only=True)
class FrameContours:
    """V_min(Z, r) at one (x, y) point in eV relative to min V_asy, with its contour levels, for a frame movie;
    V_min ends at zasymp."""

    z_bohr: np.ndarray  # (nz_potential,)
    r_bohr: np.ndarray
    v_ev: np.ndarray  # (nz_potential, nr)
    levels_ev: np.ndarray  # every FRAME_CONTOUR_STEP_EV within the energy cap, 0 among them
    z_asymptotic_bohr: float


def contours_at(
    potential: PotentialData,
    reference_hartree: float,
    hartree_to_ev: float,
    ix: int,
    iy: int,
    z_asymptotic_bohr: float,
) -> FrameContours:
    """V_min at the grid point (ix, iy) relative to min V_asy (reference_hartree) and its contour levels."""
    v_ev = relative_ev(potential.v_min_hartree[:, :, ix, iy], reference_hartree, hartree_to_ev)
    n = int(round(ENERGY_CAP_EV / FRAME_CONTOUR_STEP_EV))
    levels = FRAME_CONTOUR_STEP_EV * np.arange(-n, n + 1)  # exact multiples: the 0 level is exactly 0
    return FrameContours(
        z_bohr=potential.z_bohr,
        r_bohr=potential.r_bohr,
        v_ev=v_ev,
        levels_ev=levels,
        z_asymptotic_bohr=z_asymptotic_bohr,
    )


@dataclass(frozen=True, kw_only=True)
class BasisEdges:
    """Per record, the probability at each edge of the basis relative to the initial norm (the first record's
    total: prob_zrxy and prob_tail_zr; the run's start when its segments are joined). The top |mJ| is an edge only
    when it is below the top J (else its states are among the top J's); the top J is none when it is the only J. The
    r edges hold both grids' cells there (prob_tail_zr included: a packet that moves from the specular grid onto the
    scattering grid must not look like one that reaches the edge)."""

    times_au: np.ndarray
    top_j: np.ndarray | None  # the top J of the basis (JMax, or JMax - 1 for a homonuclear basis of odd parity)
    top_mj: np.ndarray | None  # |mJ| = the basis's largest |mJ| (the effective MJMax)
    first_r: np.ndarray  # the first EDGE_POINTS r points, both grids
    last_r: np.ndarray  # the last EDGE_POINTS r points, both grids
    first_z: np.ndarray  # the first EDGE_POINTS Z points
    reference: float  # the initial norm


def share_of_states(pop: np.ndarray, in_edge: np.ndarray) -> np.ndarray:
    """The populations (records, nrot) summed over the states with in_edge true."""
    return pop[:, in_edge].sum(axis=1)


def r_edge(frames: Frames, r_slice: slice) -> np.ndarray:
    """(records,): the probability at the r points of r_slice on both grids (prob_zrxy and prob_tail_zr)."""
    scattering = np.asarray(frames.series.arrays["prob_zrxy"][:, :, r_slice]).sum(axis=(1, 2, 3, 4))
    return scattering + np.asarray(frames.series.arrays["prob_tail_zr"][:, :, r_slice]).sum(axis=(1, 2))


def edge_shares(frames: Frames) -> BasisEdges:
    """The probability at the basis edges in every record (see BasisEdges)."""
    prob = frames.series.arrays["prob_zrxy"]
    pop = np.asarray(frames.series.arrays["pop_jmj_z"]).sum(axis=1)  # (records, nrot)
    reference = float(np.sum(prob[0]) + np.sum(frames.series.arrays["prob_tail_zr"][0]))
    j, mj = frames.states[:, 0], np.abs(frames.states[:, 1])
    top_j = share_of_states(pop, j == j.max()) / reference if j.min() < j.max() else None
    top_mj = share_of_states(pop, mj == mj.max()) / reference if mj.max() < j.max() else None
    return BasisEdges(
        times_au=frames.series.times_au,
        top_j=top_j,
        top_mj=top_mj,
        first_r=r_edge(frames, slice(0, EDGE_POINTS)) / reference,
        last_r=r_edge(frames, slice(-EDGE_POINTS, None)) / reference,
        first_z=np.asarray(prob[:, :EDGE_POINTS]).sum(axis=(1, 2, 3, 4)) / reference,
        reference=reference,
    )


@dataclass(frozen=True, kw_only=True)
class EdgeSeries:
    """One edge of the basis against time, relative to the initial norm."""

    label: str
    values: np.ndarray  # (records,)


def rotational_edges(edges: BasisEdges) -> list[EdgeSeries]:
    """The edges of the rotational basis that exist (the top J, the top |mJ|): their shares, populated while the
    packet passes the surface without that being an error."""
    series = [EdgeSeries(label="top J", values=edges.top_j), EdgeSeries(label="top |mJ|", values=edges.top_mj)]
    return [s for s in series if s.values is not None]


def grid_edge_rises(edges: BasisEdges) -> list[EdgeSeries]:
    """The grid edges (the first and last r points, the first Z points) as the rise of their share above the first
    record's (the initial state's own tail there): what reaches the edge during the run, where the FFT wraps it
    around. The verdict compares these with FRAME_EDGE_LIMIT."""
    shares = [("first r points", edges.first_r), ("last r points", edges.last_r), ("first Z points", edges.first_z)]
    return [EdgeSeries(label=label, values=values - values[0]) for label, values in shares]


def worst_record(series: EdgeSeries) -> int:
    """The record where an edge series is largest (a NaN counts as the largest)."""
    return int(np.argmax(np.nan_to_num(series.values, nan=np.inf)))


def largest_edge(series: list[EdgeSeries]) -> EdgeSeries:
    """The series that reaches the largest value in any record (a NaN counts as the largest)."""
    peaks = [edge.values[worst_record(edge)] for edge in series]
    return series[int(np.argmax(np.nan_to_num(peaks, nan=np.inf)))]
