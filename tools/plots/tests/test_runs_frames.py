"""frames/ (dtFrames) on the fixture runs, each check with its negative test: the frames sum to the norm of the
densities record at the same time; summed over (x, y) they are the densities' (Z, r) cells (psi_scat has no
component along s2s in the state (JIn, MJIn): the program projects it onto the specular grid every step, so the
cross terms with chi vanish in that sum, not point by point); the rotational populations sum to the cells;
the frames equal psi_full built here in numpy from checkpoint.bin (at the end of a run and in the middle, the first
segment of a split run), with s2s from its Cartesian formula; a split run's frames join to the unsplit run's."""

import json
from pathlib import Path

import numpy as np
import pytest

from readers.binary import read_checkpoint
from readers.densities import read_densities
from readers.errors import ReadError
from readers.frames import Frames, read_frames
from readers.run_info import read_run_info
from tests.frame_helpers import psi_full_of, relative_misses, run_with_frames_manifest, s2s_of

RUNS = ["hardwall_tiny", "dailight_oblique"]
SUM_TOLERANCE = 1e-13  # relative: sums over many cells in another order than the program's
NUMPY_TOLERANCE = 1e-12  # relative to the record's largest cell: psi_full rebuilt from the checkpoint


@pytest.mark.parametrize("run", RUNS)
def test_frames_sum_to_the_norm(fixture_runs: Path, run: str) -> None:
    """Every record: the sum of prob_zrxy and prob_tail_zr is the densities' norm at that time (frames and
    densities share the times, dtFrames = dtAve), and the analysis-surface index is the densities'. Negative: one
    cell off by 10 times the tolerance, or NaN, fails exactly its record."""
    frames = read_frames([fixture_runs / run])
    densities, _, _ = read_densities([fixture_runs / run])
    assert np.array_equal(frames.series.times_au, densities.times_au)
    assert np.array_equal(frames.series.arrays["surface_index"], densities.arrays["surface_index"])
    norm = np.asarray(densities.arrays["norm"])
    prob = np.array(frames.series.arrays["prob_zrxy"])
    tail = np.asarray(frames.series.arrays["prob_tail_zr"])
    totals = prob.sum(axis=(1, 2, 3, 4)) + tail.sum(axis=(1, 2))
    assert relative_misses(totals, norm, norm, SUM_TOLERANCE) == 0
    peak = np.unravel_index(np.argmax(prob[5]), prob[5].shape)
    for planted in (prob[5][peak] + 10 * SUM_TOLERANCE * norm[5], np.nan):
        wrong = prob.copy()
        wrong[5][peak] = planted
        totals = wrong.sum(axis=(1, 2, 3, 4)) + tail.sum(axis=(1, 2))
        assert relative_misses(totals, norm, norm, SUM_TOLERANCE) == 1


@pytest.mark.parametrize("run", RUNS)
def test_frames_summed_over_xy_are_the_density_cells(fixture_runs: Path, run: str) -> None:
    """Every record: prob_zrxy summed over (x, y) equals the densities' scattering plus specular cell at each
    (Z, r) of the scattering grid to 1e-13 of the record's largest cell, and prob_tail_zr the specular cells beyond
    it bitwise. Negative: one cell off by 10 times the tolerance fails exactly there; a tail cell one ulp off fails
    the bitwise check."""
    frames = read_frames([fixture_runs / run])
    densities, _, _ = read_densities([fixture_runs / run])
    nz = frames.nz_scattering
    summed = np.asarray(frames.series.arrays["prob_zrxy"]).sum(axis=(3, 4))
    cells = np.asarray(densities.arrays["prob_scat_zr"]) + np.asarray(densities.arrays["prob_spec_zr"])[:, :nz]
    scale = cells.max(axis=(1, 2))[:, None, None]
    assert relative_misses(summed, cells, scale, SUM_TOLERANCE) == 0
    tail = np.asarray(frames.series.arrays["prob_tail_zr"])
    spec_beyond = np.asarray(densities.arrays["prob_spec_zr"])[:, nz:]
    assert np.array_equal(tail, spec_beyond)  # the same |chi|^2 dZ dr, computed as the densities do: bitwise
    off_by_an_ulp = tail.copy()
    off_by_an_ulp[2, 0, 0] = np.nextafter(off_by_an_ulp[2, 0, 0], np.inf)
    assert not np.array_equal(off_by_an_ulp, spec_beyond)
    wrong = summed.copy()
    wrong[3, 0, 0] += 10 * SUM_TOLERANCE * scale[3, 0, 0]
    assert relative_misses(wrong, cells, scale, SUM_TOLERANCE) == 1


@pytest.mark.parametrize("run", RUNS)
def test_state_populations_sum_to_the_cells(fixture_runs: Path, run: str) -> None:
    """Every record and Z: pop_jmj_z summed over the rotational states is prob_zrxy summed over r, x, y (1e-13 of
    the record's largest Z sum). Negative: one population off by 10 times that fails exactly its Z."""
    frames = read_frames([fixture_runs / run])
    by_state = np.asarray(frames.series.arrays["pop_jmj_z"]).sum(axis=2)
    by_cell = np.asarray(frames.series.arrays["prob_zrxy"]).sum(axis=(2, 3, 4))
    scale = by_cell.max(axis=1)[:, None]
    assert relative_misses(by_state, by_cell, scale, SUM_TOLERANCE) == 0
    wrong = by_state.copy()
    wrong[-1, 10] += 10 * SUM_TOLERANCE * scale[-1, 0]
    assert relative_misses(wrong, by_cell, scale, SUM_TOLERANCE) == 1


def frame_from_checkpoint(run: Path, frames: Frames, s2s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(prob_zrxy, pop_jmj_z) of psi_full built from the run's checkpoint.bin with the given s2s (nx, ny)."""
    psi, _ = psi_full_of(run, frames, s2s)
    g = read_run_info(run).data["grids"]
    density = np.abs(psi) ** 2 * (g["dz"] * g["dr"] * g["xy_volume_element"])  # (nz, nr, nx, ny, nrot)
    return density.sum(axis=4), density.sum(axis=(1, 2, 3))


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique", "dailight_oblique_split/first"])
def test_frames_equal_numpy_on_the_checkpoint(fixture_runs: Path, run: str) -> None:
    """The last record (at the checkpoint's time: tEnd, or mid-run for the split run's first segment) equals
    psi_full rebuilt in numpy, to 1e-12 of the record's largest value."""
    directory = fixture_runs / run
    frames = read_frames([directory])
    info = read_run_info(directory)
    assert frames.series.times_au[-1] == read_checkpoint(directory / "checkpoint.bin").time_au
    prob, pop = frame_from_checkpoint(directory, frames, s2s_of(info, frames.nx, frames.ny))
    got_prob = np.asarray(frames.series.arrays["prob_zrxy"][-1])
    got_pop = np.asarray(frames.series.arrays["pop_jmj_z"][-1])
    assert relative_misses(got_prob, prob, got_prob.max(), NUMPY_TOLERANCE) == 0
    assert relative_misses(got_pop, pop, got_pop.max(), NUMPY_TOLERANCE) == 0


def test_conjugated_s2s_fails_where_the_cross_terms_matter(fixture_runs: Path) -> None:
    """Negative test of the numpy check on the oblique run (k_par != 0, mid-run: the packet is at the surface): with
    conj(s2s) the frame differs exactly at the cells where the cross term 2 Re(psi* s2s chi) changes by more than the
    tolerance, and the largest change is a sizeable part of the largest cell (more than 1e-2 of it)."""
    directory = fixture_runs / "dailight_oblique_split/first"
    frames = read_frames([directory])
    s2s = s2s_of(read_run_info(directory), frames.nx, frames.ny)
    right, _ = frame_from_checkpoint(directory, frames, s2s)
    wrong, _ = frame_from_checkpoint(directory, frames, np.conj(s2s))
    got = np.asarray(frames.series.arrays["prob_zrxy"][-1])
    scale = got.max()
    changed = ~(np.abs(wrong - right) <= NUMPY_TOLERANCE * scale)
    assert np.max(np.abs(wrong - right)) > 1e-2 * scale
    assert relative_misses(got, wrong, scale, NUMPY_TOLERANCE) == np.count_nonzero(changed)


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique"])
def test_split_run_joins_to_the_unsplit_frames(fixture_runs: Path, run: str) -> None:
    """The two segments of a split run, joined, hold the unsplit run's records bitwise (the record at the split
    time appears once)."""
    split = fixture_runs / f"{run}_split"
    joined = read_frames([split / "second", split / "first"])  # out of order: the reader sorts by tStart
    whole = read_frames([fixture_runs / run])
    assert np.array_equal(joined.series.times_au, whole.series.times_au)
    for name in ("prob_zrxy", "prob_tail_zr", "pop_jmj_z", "surface_index"):
        assert np.array_equal(np.asarray(joined.series.arrays[name]), np.asarray(whole.series.arrays[name])), name


def test_reader_refuses_frames_it_cannot_index(fixture_runs: Path, tmp_path: Path) -> None:
    """frames/ whose grids do not give the datasets' shapes (nx and ny swapped on the oblique 4 x 4 grid cannot show
    it, so nx doubled), without its rotational states or with a dataset removed is a ReadError naming what is wrong;
    the unchanged copy reads."""
    fixture = fixture_runs / "dailight_oblique"
    manifest = json.loads((fixture / "frames" / "manifest.json").read_text(encoding="utf-8"))
    read_frames([run_with_frames_manifest(fixture, tmp_path / "same", {})])
    datasets = {k: v for k, v in manifest["datasets"].items() if k != "pop_jmj_z"}
    cases = {
        "prob_zrxy has records of shape": {"grids": {**manifest["grids"], "nx": 2 * manifest["grids"]["nx"]}},
        "rotational_states": {"rotational_states": None},
        "no pop_jmj_z": {"datasets": datasets},
    }
    for k, (message, change) in enumerate(cases.items()):
        run = run_with_frames_manifest(fixture, tmp_path / f"case{k}", change)
        with pytest.raises(ReadError, match=message):
            read_frames([run])
