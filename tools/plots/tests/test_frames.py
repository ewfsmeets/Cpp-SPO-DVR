"""compute/frames.py on a synthetic frames/: the (x, y) probability in a (Z, r) window, the J populations and the
basis edges, each against a hand count, with a slip that a test must catch."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from compute.frames import (
    EDGE_POINTS,
    MAX_FRAME_SITES,
    EdgeSeries,
    edge_shares,
    grid_edge_rises,
    j_populations,
    largest_edge,
    most_at_surface,
    rotational_edges,
    worst_record,
    xy_probability,
    zr_probability_at,
)
from figures import FramePoint, frame_points
from readers.arrays import Series
from readers.frames import Frames
from readers.sites import Site


def synthetic_frames() -> Frames:
    """Two records on 4 Z x 3 r x 2 x 1 (x, y) points and the states (J, mJ) = (0, 0), (1, -1), (1, 0), (1, 1):
    prob_zrxy = (1 + iz + 10 ir) in record 0 and twice that in record 1 (a Z and an r window give different
    sums); the populations 1, 2, 3, 4 per Z in both records (deliberately not prob_zrxy's sums: each function
    reads one dataset); prob_tail_zr 0.5 per cell."""
    iz, ir = np.meshgrid(np.arange(4.0), np.arange(3.0), indexing="ij")
    cells = np.repeat((1.0 + iz + 10.0 * ir)[:, :, None, None], 2, axis=2)  # (4, 3, 2, 1)
    arrays = {
        "prob_zrxy": np.stack([cells, 2 * cells]),
        "prob_tail_zr": 0.5 * np.ones((2, 1, 3)),
        "pop_jmj_z": np.tile(np.array([1.0, 2.0, 3.0, 4.0]), (2, 4, 1)),
        "surface_index": np.zeros(2),
    }
    series = Series(Path("frames"), {}, np.array([0.0, 100.0]), arrays, None)
    return Frames(
        series=series,
        z_bohr=np.arange(5.0),
        r_bohr=np.arange(3.0),
        nz_scattering=4,
        nx=2,
        ny=1,
        states=np.array([[0, 0], [1, -1], [1, 0], [1, 1]]),
        initial_index=0,
        sites=[],
        starts_au=[0.0],
    )


def test_frame_points_are_the_sites_or_the_origin() -> None:
    """The frame movies show the first MAX_FRAME_SITES sites, their stems numbered (a name may hold any character),
    or, without sites, the cell's origin."""
    frames = synthetic_frames()
    assert frame_points(frames) == [FramePoint(stem="frame_movie_origin", name="the cell's origin", ix=0, iy=0)]
    sites = [Site(f"site/{k}", 0.0, 0.0, k % 2, 0) for k in range(MAX_FRAME_SITES + 1)]
    points = frame_points(replace(frames, sites=sites))
    assert [point.stem for point in points] == [f"frame_movie_site{k + 1}" for k in range(MAX_FRAME_SITES)]
    assert points[1] == FramePoint(stem="frame_movie_site2", name="site/1", ix=1, iy=0)


def test_windows_and_points() -> None:
    """Record 1 over Z indices 1..2 and r index 2: per (x, y) point 2 (1 + 1 + 20) + 2 (1 + 2 + 20) = 90; with
    Z and r swapped (Z index 2, r indices 1..2: a slip) 2 (1 + 2 + 10) + 2 (1 + 2 + 20) = 72. At (x, y) point
    (1, 0) one value per record, Z and r; the record with the most probability in the first 2 Z points is 1."""
    frames = synthetic_frames()
    assert xy_probability(frames, 1, slice(1, 3), slice(2, 3)).tolist() == [[90.0], [90.0]]
    assert xy_probability(frames, 1, slice(2, 3), slice(1, 3)).tolist() == [[72.0], [72.0]]
    assert zr_probability_at(frames, 1, 0).shape == (2, 4, 3)
    assert most_at_surface(frames, 2) == 1


def test_j_populations() -> None:
    """J = 0 holds 1 per Z (4 Z: 4) plus chi beyond the scattering grid (prob_tail_zr: 3 x 0.5), since the initial
    state (index 0) has J = 0: 5.5; J = 1 holds 2 + 3 + 4 per Z (36); in both records. Leaving the tail out (a
    slip) gives 4, adding it to the wrong J 37.5; summing over mJ instead of J would give three values."""
    j_values, per_j = j_populations(synthetic_frames())
    assert j_values.tolist() == [0, 1]
    assert per_j.tolist() == [[5.5, 36.0], [5.5, 36.0]]


def test_edge_shares() -> None:
    """The reference is the initial norm (record 0's cells and tail: 300 + 1.5); the top J (J = 1) holds 2 + 3 + 4
    per Z (36); its |mJ| = 1 is the top J, so no |mJ| edge (its states are the top J's); the first and last
    EDGE_POINTS r points hold their cells' sums on both grids (the tail adds 2 x 0.5: leaving it out, a slip, makes a
    packet that moves onto the scattering grid look like one reaching the edge), the first EDGE_POINTS Z points
    their cells' sums."""
    frames = synthetic_frames()
    edges = edge_shares(frames)
    prob = frames.series.arrays["prob_zrxy"]
    initial = float(prob[0].sum()) + 1.5
    assert edges.reference == initial == 301.5
    assert edges.top_j.tolist() == [36 / initial] * 2 and edges.top_mj is None
    assert edges.first_r.tolist() == [(float(prob[k, :, :EDGE_POINTS].sum()) + 1.0) / initial for k in (0, 1)]
    assert edges.last_r.tolist() == [(float(prob[k, :, -EDGE_POINTS:].sum()) + 1.0) / initial for k in (0, 1)]
    assert edges.first_z.tolist() == [float(prob[k, :EDGE_POINTS].sum()) / initial for k in (0, 1)]
    assert [edge.label for edge in rotational_edges(edges)] == ["top J"]


def test_grid_edge_rises_and_the_largest_edge() -> None:
    """A grid edge's rise is its share minus the first record's (record 1 doubles the cells, so the rise is the
    record-0 cell share again; the share itself, a slip, would be twice that); the largest series and its worst
    record are found with NaN counted as the largest."""
    edges = edge_shares(synthetic_frames())
    rises = grid_edge_rises(edges)
    assert [edge.label for edge in rises] == ["first r points", "last r points", "first Z points"]
    assert rises[2].values.tolist() == [0.0, edges.first_z[0]]
    series = [
        EdgeSeries(label="a", values=np.array([0.0, 2.0, 1.0])),
        EdgeSeries(label="b", values=np.array([3.0, 0.0, 0.0])),
    ]
    assert largest_edge(series).label == "b" and worst_record(series[0]) == 1
    series[0].values[2] = np.nan
    assert largest_edge(series).label == "a" and worst_record(series[0]) == 2


def test_mj_edge_below_the_top_j() -> None:
    """States (0, 0), (1, 1), (2, 0), (2, 1) (|mJ| <= 1 < J = 2): the top J holds 3 + 4 per Z (28), the top |mJ|
    = 1 holds 2 + 4 (24); taking |mJ| = 0 or the J = 1 states (slips) gives other sums. A basis of one J (3, with
    |mJ| up to 2) has an |mJ| edge (4 per Z: 16) and no J edge (None, left out of rotational_edges)."""
    frames = replace(synthetic_frames(), states=np.array([[0, 0], [1, 1], [2, 0], [2, 1]]))
    edges = edge_shares(frames)
    assert edges.top_j.tolist() == [28 / edges.reference] * 2
    assert edges.top_mj.tolist() == [24 / edges.reference] * 2
    one_j = edge_shares(replace(frames, states=np.array([[3, -1], [3, 0], [3, 1], [3, 2]])))
    assert one_j.top_j is None and one_j.top_mj.tolist() == [16 / one_j.reference] * 2
    assert [edge.label for edge in rotational_edges(one_j)] == ["top |mJ|"]
