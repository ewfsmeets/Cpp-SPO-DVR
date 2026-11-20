"""norm_budget/ (NormBudget) on the fixture runs, each check with its negative test: the losses and their error
indicators equal a numpy recomputation from the densities/ cells and the absorbers of run_info.json (hardwall_budget:
a densities record at every window's end); the budget closes (the residual within its error estimate) on the
hard-wall fixtures, and the indicator falls with dtAna squared; on the coupled Dai-Light fixture it opens while the
packet is at the surface (probability lost in the angular transform, docs/output_files.md "norm_budget/"); the window
norms chain and equal the densities' norms; a split run joins to the unsplit one."""

import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pytest

from compute.norm_budget import budget_check, excess_over_error, summed_budget
from readers.densities import read_densities
from readers.errors import ReadError
from readers.norm_budget import read_norm_budget
from readers.run_info import RunInfo, read_run_info

NUMPY_TOLERANCE = 1e-12  # relative to the largest loss of the run: sums over cells in another order
HARD_WALL_RUNS = ["hardwall_tiny", "hardwall_budget", "hardwall_budget_dt4"]  # no angular coupling
SURFACE_REACHED_AU = 2000.0  # dailight_oblique: the packet reaches Z < zasymp after this (densities/)
SURFACE_LEFT_AU = 6000.0  # ... and has left it by this time


def absorbers_on(info: RunInfo, coordinate: str, grid: str) -> list[dict]:
    """The absorbers of run_info.json on `coordinate` ("Z" or "r") that act on `grid` (ROpt's grid is "both")."""
    return [a for a in info.data["absorbers"]["list"] if a["coordinate"] == coordinate and a["grid"] in (grid, "both")]


def absorber_factors(info: RunInfo, absorbers: list[dict], points: np.ndarray) -> np.ndarray:
    """exp(-V_opt dt) per time step at the grid points (1 outside), from each absorber's start, end, strength and the
    power p: V_opt = strength (p + 1) / 2 x^p, x rising from 0 at the start of a high-end absorber to 1 at its end
    (falling for a low-end one, *Low; a later absorber overwrites an earlier one where they meet). The formula of
    src/propagation/propagation_tables.cpp set_optical_factors (Fortran: Absorption.f90 SetupAbsorpPot)."""
    power, dt = info.data["absorbers"]["power"], info.data["times"]["dt_step"]
    factors = np.ones(len(points))
    for a in absorbers:
        inside = (points >= a["start"]) & (points <= a["end"])
        length = a["end"] - a["start"]
        x = (a["end"] - points) / length if a["name"].endswith("Low") else (points - a["start"]) / length
        factors = np.where(inside, np.exp(-(a["strength"] * 0.5 * (power + 1) * x**power) * dt), factors)
    return factors


@dataclass(frozen=True, kw_only=True)
class GridLosses:
    """One grid's losses over a window: the Z absorbers' part per Z point and the r absorber's part, for the loss, its
    error indicator and the spread of the two one-sided estimates (a first-order quantity, for a negative test)."""

    loss_by_z: np.ndarray  # (nz,)
    loss_r: float
    error_by_z: np.ndarray  # (nz,)
    error_r: float
    spread: float  # summed over the grid's cells: |end (1/F - 1) - start (1 - F)|


def grid_losses(start: np.ndarray, end: np.ndarray, f_z: np.ndarray, f_r: np.ndarray, factor_power: int) -> GridLosses:
    """The losses of one grid's cells (nz, nr) from their probabilities at the window's start and end: gamma = -ln F
    with F = (f_z f_r)^factor_power the density's factor over the window; per cell the constant-inflow loss
    start (1 - F) + (end - F start) (gamma / (1 - F) - 1) and its indicator gamma^2 / 6 |end - start + loss|
    (docs/output_files.md "norm_budget/"), shared between the Z and the r absorber as ln f_z : ln f_r."""
    f = np.outer(f_z, f_r) ** factor_power
    gamma = -np.log(f)
    absorbing = gamma > 0
    weight = np.divide(gamma, 1.0 - f, out=np.ones_like(f), where=absorbing) - 1.0
    loss = start * (1.0 - f) + (end - f * start) * weight
    error = gamma**2 / 6.0 * np.abs(end - start + loss)
    v_z, v_r = -np.log(f_z)[:, None], -np.log(f_r)[None, :]
    z_share = np.divide(v_z, v_z + v_r, out=np.zeros_like(f), where=absorbing)
    return GridLosses(
        loss_by_z=(loss * z_share).sum(axis=1),
        loss_r=float((loss * (1.0 - z_share)).sum()),
        error_by_z=(error * z_share).sum(axis=1),
        error_r=float((error * (1.0 - z_share)).sum()),
        spread=float(np.abs(end * (1.0 / f - 1.0) - start * (1.0 - f)).sum()),
    )


def recomputed_losses(run: Path, factor_per_step: int = 2) -> tuple[np.ndarray, np.ndarray, float]:
    """(windows, absorbers) losses and error indicators recomputed from densities/ (record k is window k's start,
    record k + 1 its end) and run_info.json, in the order of norm_budget/'s absorbers, and the one-sided estimates'
    spread summed over the run. factor_per_step 2: the density shrinks by f^2 per time step (4 applies the factor
    twice: the negative test)."""
    info, records = read_run_info(run), read_norm_budget([run])
    series, grids, _ = read_densities([run])
    assert np.array_equal(series.times_au[1:], records.series.times_au)
    z_by_grid = {"scattering": grids.z_bohr[: grids.nz_scattering], "specular": grids.z_bohr}
    cells_by_grid = {"scattering": series.arrays["prob_scat_zr"], "specular": series.arrays["prob_spec_zr"]}
    loss = np.zeros((len(records.series.times_au), len(records.absorbers)))
    error = np.zeros_like(loss)
    spread = 0.0
    for grid, z in z_by_grid.items():
        z_absorbers = absorbers_on(info, "Z", grid)
        f_z = absorber_factors(info, z_absorbers, z)
        f_r = absorber_factors(info, absorbers_on(info, "r", grid), grids.r_bohr)
        for k in range(len(records.series.times_au)):
            cells = cells_by_grid[grid]
            exponent = factor_per_step * records.steps_per_window
            g = grid_losses(np.asarray(cells[k]), np.asarray(cells[k + 1]), f_z, f_r, exponent)
            spread += g.spread
            for j, absorber in enumerate(records.absorbers):
                if absorber.grid != grid:
                    continue
                if absorber.name == "ROpt":
                    loss[k, j], error[k, j] = g.loss_r, g.error_r
                    continue
                (a,) = [a for a in z_absorbers if a["name"] == absorber.name]
                points = np.flatnonzero((z >= a["start"]) & (z <= a["end"]))
                assert [points[0], points[-1]] == [a["points"]["first"], a["points"]["last"]]  # run_info's points
                loss[k, j], error[k, j] = g.loss_by_z[points].sum(), g.error_by_z[points].sum()
    return loss, error, spread


def test_losses_equal_numpy_on_the_densities(fixture_runs: Path) -> None:
    """hardwall_budget (dtAve = dtAna: densities/ holds every window's start and end cells): the loss and its error
    indicator of every window and absorber equal the numpy recomputation to 1e-12 of the run's largest loss, and the
    fixture loses something (the check is not empty). Negative: the factor applied twice per step fails at every
    window with a loss above 1e-6 of the largest."""
    run = fixture_runs / "hardwall_budget"
    records = read_norm_budget([run])
    loss, error, _ = recomputed_losses(run)
    got_loss = np.asarray(records.series.arrays["loss"])
    got_error = np.asarray(records.series.arrays["loss_error"])
    scale = float(np.max(got_loss))
    assert scale > 1e-6
    assert np.max(np.abs(got_loss - loss)) <= NUMPY_TOLERANCE * scale
    assert np.max(np.abs(got_error - error)) <= NUMPY_TOLERANCE * scale
    twice, _, _ = recomputed_losses(run, factor_per_step=4)
    missed = np.abs(got_loss - twice) > NUMPY_TOLERANCE * scale
    losing = got_loss > 1e-6 * scale
    assert np.array_equal(missed & losing, losing) and np.any(losing)


@pytest.mark.parametrize("run", HARD_WALL_RUNS)
def test_budget_closes(fixture_runs: Path, run: str) -> None:
    """The residual (the norm that left the grids minus the absorbers' summed losses) closes as the verdict checks it
    (budget_check: within the error estimate plus BUDGET_LIMIT), and, stricter, within the error estimate plus 1e-8
    at every window end (the time stepping's own error stays below that here; the verdict's limit alone would hide
    a wrong indicator). Negative: leaving out the absorber that took the most opens it."""
    budget = summed_budget(read_norm_budget([fixture_runs / run]))
    assert budget_check(budget).closes
    assert np.max(np.abs(budget.residual) - budget.error_estimate) <= 1e-8
    largest = int(np.argmax(budget.lost_by_absorber[-1]))
    without = budget.left_grids - np.delete(budget.lost_by_absorber, largest, axis=1).sum(axis=1)
    assert not budget_check(replace(budget, residual=without)).closes


def test_coupled_budget_opens_at_the_surface(fixture_runs: Path) -> None:
    """dailight_oblique (JMax 2): the budget closes until the packet reaches the surface, then the residual grows by
    more than 1e-3 of the norm by the time it has left and exceeds its error estimate: probability that the absorbers
    did not take (the angular transform's, see the module's docstring)."""
    budget = summed_budget(read_norm_budget([fixture_runs / "dailight_oblique"]))
    excess = excess_over_error(budget)
    reached = int(np.argmin(np.abs(budget.times_au - SURFACE_REACHED_AU)))
    left = int(np.argmin(np.abs(budget.times_au - SURFACE_LEFT_AU)))
    assert np.max(excess[: reached + 1]) <= 0.0
    assert not budget_check(budget).closes
    assert budget.residual[left] - budget.residual[reached] > 1e-3


def test_error_estimate_falls_with_dt_analysis_squared(fixture_runs: Path) -> None:
    """The error indicator summed over the hard-wall run falls by about (10 / 4)^2 = 6.25 from dtAna 10 to dtAna 4
    (accepted 4 to 9): gamma^2 / 6 of the inflow. The residual of the constant-inflow loss is smaller than the
    indicator in both runs. Negative: the one-sided estimates' spread, a first-order quantity recomputed from the
    same densities, falls only about 10 / 4, outside the accepted range."""
    long_windows = summed_budget(read_norm_budget([fixture_runs / "hardwall_budget"]))
    short_windows = summed_budget(read_norm_budget([fixture_runs / "hardwall_budget_dt4"]))
    ratio = long_windows.error_estimate[-1] / short_windows.error_estimate[-1]
    assert 4.0 <= ratio <= 9.0, ratio
    for budget in (long_windows, short_windows):
        assert abs(budget.residual[-1]) < budget.error_estimate[-1]
    _, _, long_spread = recomputed_losses(fixture_runs / "hardwall_budget")
    _, _, short_spread = recomputed_losses(fixture_runs / "hardwall_budget_dt4")
    assert not 4.0 <= long_spread / short_spread <= 9.0, long_spread / short_spread


@pytest.mark.parametrize("run", ["hardwall_budget", "dailight_oblique"])
def test_window_norms_chain_and_equal_the_densities(fixture_runs: Path, run: str) -> None:
    """Each window starts where the one before ended (the same cells: bitwise); the norm at a window end equals the
    densities' norm at that time where densities/ has a record (every window for hardwall_budget, every dtAve for
    dailight_oblique) to 1e-13 (the same cells summed in another order). Negative: a norm one part in 1e-12 off
    fails."""
    records = read_norm_budget([fixture_runs / run])
    arrays = records.series.arrays
    assert np.array_equal(np.asarray(arrays["norm_start"])[1:], np.asarray(arrays["norm_end"])[:-1])
    densities, _, _ = read_densities([fixture_runs / run])
    common, at_budget, at_densities = np.intersect1d(
        records.series.times_au, densities.times_au, assume_unique=True, return_indices=True
    )
    assert len(common) >= 10
    norm_end = np.asarray(arrays["norm_end"])[at_budget]
    norm = np.asarray(densities.arrays["norm"])[at_densities]
    assert np.max(np.abs(norm_end - norm) / norm) <= 1e-13
    off = norm_end.copy()
    off[3] *= 1 + 1e-12
    assert np.count_nonzero(~(np.abs(off - norm) / norm <= 1e-13)) == 1


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique"])
def test_split_run_joins_to_the_unsplit_budget(fixture_runs: Path, run: str) -> None:
    """The two segments of a split run, joined, hold the unsplit run's windows bitwise (the second segment's first
    window starts at its tStart, so no window repeats)."""
    split = fixture_runs / f"{run}_split"
    joined = read_norm_budget([split / "second", split / "first"])
    whole = read_norm_budget([fixture_runs / run])
    assert np.array_equal(joined.series.times_au, whole.series.times_au)
    for name in ("loss", "loss_error", "norm_start", "norm_end"):
        assert np.array_equal(np.asarray(joined.series.arrays[name]), np.asarray(whole.series.arrays[name])), name


def budget_copy(fixture: Path, work: Path, manifest: dict) -> Path:
    """A copy of a fixture run (its files linked) whose norm_budget/manifest.json is `manifest`."""
    work.mkdir()
    for path in fixture.iterdir():
        if path.name != "norm_budget":
            (work / path.name).symlink_to(path.resolve())
    (work / "norm_budget").mkdir()
    for path in (fixture / "norm_budget").iterdir():
        if path.name != "manifest.json":
            (work / "norm_budget" / path.name).symlink_to(path.resolve())
    (work / "norm_budget" / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return work


def test_reader_refuses_a_budget_it_cannot_index(fixture_runs: Path, tmp_path: Path) -> None:
    """norm_budget/ without its absorbers, with an absorber too many (the loss datasets then hold one value too few)
    or without a dataset is a ReadError naming what is wrong; the unchanged copy reads."""
    fixture = fixture_runs / "hardwall_tiny"
    manifest = json.loads((fixture / "norm_budget" / "manifest.json").read_text(encoding="utf-8"))
    read_norm_budget([budget_copy(fixture, tmp_path / "same", manifest)])
    without_absorbers = {k: v for k, v in manifest.items() if k != "absorbers"}
    extra_absorber = {**manifest, "absorbers": [*manifest["absorbers"], {"name": "ZOptLow", "grid": "scattering"}]}
    without_norm_end = {**manifest, "datasets": {k: v for k, v in manifest["datasets"].items() if k != "norm_end"}}
    cases = {"absorbers": without_absorbers, "records of shape": extra_absorber, "no norm_end": without_norm_end}
    for k, (message, changed) in enumerate(cases.items()):
        with pytest.raises(ReadError, match=message):
            read_norm_budget([budget_copy(fixture, tmp_path / f"case{k}", changed)])
