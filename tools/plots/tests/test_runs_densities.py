"""densities/ (writeDensities) on the fixture runs, each check with its negative test: the cells sum to the
averages' norm, give its mean positions, and partition it into the regions; the records equal the
checkpoints' wave functions; a split run's series joins to the unsplit one's, bitwise."""

import shutil
from pathlib import Path

import numpy as np
import pytest

from compute.densities import (
    DensityGrids,
    band_edges,
    band_norms,
    combined_probability,
    density_grids,
    mean_positions_bohr,
    r_absorber_first_index,
    region_norms,
)
from readers.arrays import read_series
from readers.binary import read_checkpoint
from readers.errors import ReadError
from readers.run_info import RunInfo, read_run_info
from readers.segments import read_joined_series
from readers.tables import read_averages

RUNS = ["hardwall_tiny", "dailight_oblique"]
SUM_TOLERANCE = 1e-13  # relative: the program sums the cells in another order and scales the sum once
MEAN_TOLERANCE = 1e-12  # relative: <Z>, <r> are ratios of such sums
EDGE_BOHR = 1e-9  # a grid point this close to a region's edge lies on it (grid points are rounded)


def relative_norm_error(total: np.ndarray, norm: np.ndarray) -> np.ndarray:
    """|sum of cells - norm| / norm per record (NaN stays NaN)."""
    return np.abs(np.sum(total, axis=(1, 2)) - norm) / norm


@pytest.mark.parametrize("run", RUNS)
def test_cells_sum_to_the_norm(fixture_runs: Path, run: str) -> None:
    """Both grids' cells sum to averages.dat's norm at every record (same times, bitwise); negative tests: one
    cell scaled by 1 + 1e-10, or NaN, fails exactly its record."""
    series = read_series(fixture_runs / run / "densities")
    averages = read_averages(fixture_runs / run / "averages.dat")
    assert series.times_au.tobytes() == averages.times_au.tobytes()
    total = combined_probability(series)
    assert np.max(relative_norm_error(total, averages.norm)) <= SUM_TOLERANCE
    peak = np.unravel_index(np.argmax(total[3]), total[3].shape)
    for planted in (total[3][peak] * (1 + 1e-10), np.nan):
        wrong = total.copy()
        wrong[3][peak] = planted
        failing = ~(relative_norm_error(wrong, averages.norm) <= SUM_TOLERANCE)
        assert list(np.flatnonzero(failing)) == [3]


@pytest.mark.parametrize("run", RUNS)
def test_cells_give_the_mean_positions(fixture_runs: Path, run: str) -> None:
    """<Z> and <r> from the cells equal averages.dat's; negative test: the specular cells beyond the scattering
    grid given the scattering grid's last Z (a realistic slip) fail where the specular wave is out there."""
    series = read_series(fixture_runs / run / "densities")
    averages = read_averages(fixture_runs / run / "averages.dat")
    grids = density_grids(series)
    total = combined_probability(series)
    z_mean, r_mean = mean_positions_bohr(total, grids)
    assert np.max(np.abs(z_mean - averages.z_mean_bohr) / np.abs(averages.z_mean_bohr)) <= MEAN_TOLERANCE
    assert np.max(np.abs(r_mean - averages.r_mean_bohr) / averages.r_mean_bohr) <= MEAN_TOLERANCE
    nz = grids.nz_scattering
    grids.z_bohr[nz:] = grids.z_bohr[nz - 1]
    wrong_z, _ = mean_positions_bohr(total, grids)
    beyond = np.sum(series.arrays["prob_spec_zr"][:, nz:, :], axis=(1, 2)) > 1e-6
    assert beyond.any() and np.all(
        np.abs(wrong_z - averages.z_mean_bohr)[beyond] / averages.z_mean_bohr[beyond] > MEAN_TOLERANCE
    )


def regions_by_value(total: np.ndarray, grids: DensityGrids, info: RunInfo) -> list[np.ndarray]:
    """The regions of region_norms selected here by coordinate values (run_info.json's ROptStart, zasymp and
    Zana, the scattering grid's last point), in its order: r >= r_c, then for r < r_c the four Z bands."""
    r_c = next(a["start"] for a in info.data["absorbers"]["list"] if a["coordinate"] == "r")
    outer = grids.r_bohr >= r_c - EDGE_BOHR
    z, z_asymptotic, z_analysis = grids.z_bohr, info.data["z_asymptotic"], info.data["analysis"]["z_ana"]
    z_end = z[grids.nz_scattering - 1]
    bands = [
        z < z_asymptotic - EDGE_BOHR,
        (z >= z_asymptotic - EDGE_BOHR) & (z < z_analysis - EDGE_BOHR),
        (z >= z_analysis - EDGE_BOHR) & (z <= z_end + EDGE_BOHR),
        z > z_end + EDGE_BOHR,
    ]
    return [np.sum(total[:, :, outer], axis=(1, 2))] + [
        np.sum(total[:, band][:, :, ~outer], axis=(1, 2)) for band in bands
    ]


@pytest.mark.parametrize("run", RUNS)
def test_regions_partition_the_norm(fixture_runs: Path, run: str) -> None:
    """The regions sum to the norm, and each equals the probability selected by coordinate values; negative
    test: the edge between the bands zasymp..Zana and Zana..end moved one point down, through the region code,
    fails exactly at the records where that point holds probability (r < r_c)."""
    series = read_series(fixture_runs / run / "densities")
    info = read_run_info(fixture_runs / run)
    grids = density_grids(series)
    norm = np.asarray(series.arrays["norm"])
    regions = region_norms(series, info, grids)
    assert np.max(np.abs(sum(regions.regions.values()) - norm) / norm) <= SUM_TOLERANCE
    total = combined_probability(series)
    for computed, expected in zip(regions.regions.values(), regions_by_value(total, grids, info), strict=True):
        assert np.max(np.abs(computed - expected) / norm) <= SUM_TOLERANCE
    ir_c = r_absorber_first_index(info, len(grids.r_bohr))
    bands = band_edges(info, grids)
    (low, first, _), (high, edge, end) = bands[1], bands[2]
    moved = [*bands[:1], (low, first, edge - 1), (high, edge - 1, end), *bands[3:]]
    wrong = list(band_norms(total, moved, ir_c).values())
    expected = regions_by_value(total, grids, info)[1:]
    failing = np.zeros(len(norm), dtype=bool)
    for computed, want in zip(wrong, expected, strict=True):
        failing |= ~(np.abs(computed - want) / norm <= SUM_TOLERANCE)
    holds = np.sum(total[:, edge - 1, :ir_c], axis=1) / norm > SUM_TOLERANCE
    assert holds.any() and np.array_equal(failing, holds)


def checkpoint_cells(run: Path) -> np.ndarray:
    """The combined probability per (Z, r) cell of the run's checkpoint.bin, computed here from psi."""
    info = read_run_info(run)
    checkpoint = read_checkpoint(run / "checkpoint.bin")
    g = info.data["grids"]
    scattering = np.sum(np.abs(np.asarray(checkpoint.scattering)) ** 2, axis=(1, 2)).reshape(g["nz"], g["nr"])
    cells = np.abs(np.asarray(checkpoint.specular)) ** 2 * g["dz"] * g["dr"]
    cells[: g["nz"]] += scattering * g["dz"] * g["dr"] * g["xy_volume_element"]
    return cells


@pytest.mark.parametrize("run", [*RUNS, "dailight_oblique_split/first"])
def test_last_record_is_the_checkpoints_wave(fixture_runs: Path, run: str) -> None:
    """The last record is at the checkpoint's time (tEnd, a multiple of dtAve in the fixtures) and equals the
    probability per cell of checkpoint.bin's wave function, to 1e-12 of its largest cell (the first segment of
    the split run ends while the packet is at the surface). Negative: the largest cell off by 10 times the
    tolerance, or NaN, fails exactly that cell."""
    series = read_series(fixture_runs / run / "densities")
    assert series.times_au[-1] == read_checkpoint(fixture_runs / run / "checkpoint.bin").time_au
    expected = checkpoint_cells(fixture_runs / run)
    record = np.asarray(combined_probability(series)[-1])
    tolerance = 1e-12 * np.max(expected)
    assert np.max(np.abs(record - expected)) <= tolerance
    peak = np.unravel_index(np.argmax(expected), expected.shape)
    for planted in (record[peak] + 10 * tolerance, np.nan):
        wrong = record.copy()
        wrong[peak] = planted
        assert np.count_nonzero(~(np.abs(wrong - expected) <= tolerance)) == 1


@pytest.mark.parametrize("name", RUNS)
def test_split_series_joins_to_the_unsplit_one(fixture_runs: Path, tmp_path: Path, name: str) -> None:
    """The two segments' densities/, given in the wrong order (the join orders them), joined equal the unsplit
    run's bitwise, every dataset; the second segment's first record repeats the first's last (the split lies on
    the record grid), so the repeat is dropped. Negative test: that repeated record changed in one cell is
    refused."""
    split = fixture_runs / f"{name}_split"
    first, second = read_series(split / "first" / "densities"), read_series(split / "second" / "densities")
    assert first.times_au[-1] == second.times_au[0]
    joined, starts = read_joined_series([split / "second", split / "first"], "densities")
    whole = read_series(fixture_runs / name / "densities")
    assert joined.times_au.tobytes() == whole.times_au.tobytes() and len(starts) == 2
    for dataset, values in whole.arrays.items():
        assert np.asarray(joined.arrays[dataset]).tobytes() == np.asarray(values).tobytes(), dataset
    for segment in ("first", "second"):
        shutil.copytree(split / segment, tmp_path / segment, ignore=shutil.ignore_patterns("checkpoint.bin", "fort.*"))
    changed = tmp_path / "second" / "densities" / "prob_spec_zr.bin"
    cells = np.fromfile(changed, "<f8")
    cells[np.argmax(cells[: second.arrays["prob_spec_zr"][0].size])] *= 1 + 1e-15
    cells.tofile(changed)
    with pytest.raises(ReadError, match="differs between two segments"):
        read_joined_series([tmp_path / "first", tmp_path / "second"], "densities")
