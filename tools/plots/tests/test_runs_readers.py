"""The readers on the fixture runs: every file reads, with the record counts its run implies, and the
checkpoint's wave function gives the norm of the last averages.dat row."""

from pathlib import Path

import numpy as np
import pytest

from readers.binary import Checkpoint, read_bk_coefficients, read_checkpoint
from readers.fine_analysis import read_fine_grid
from readers.run_info import read_run_info
from readers.tables import read_averages, read_gas_phase_states, read_scattering_summary

RUNS = ["hardwall_tiny", "dailight_oblique", "hardwall_tiny_split/second", "dailight_oblique_split/first"]


@pytest.mark.parametrize("run", RUNS)
def test_record_counts(fixture_runs: Path, run: str) -> None:
    """averages.dat has one row per dtAve from tStart to tEnd, the BK file one record per analysis step up to
    tEnd (a restart's file holds the earlier segments' records too), and nothing is NaN."""
    directory = fixture_runs / run
    info = read_run_info(directory)
    times = info.data["times"]
    averages = read_averages(directory / "averages.dat")
    expected_rows = round((info.t_end_au - info.t_start_au) / times["dt_averages"]) + 1
    assert len(averages.times_au) == expected_rows
    assert averages.first_non_finite is None
    bk = read_bk_coefficients(directory / "bk_coefficients.bin")
    dt_analysis = times["steps_per_analysis"] * times["dt_step"]
    assert len(bk.times_au) == round(info.t_end_au / dt_analysis)
    assert bk.first_non_finite is None
    assert bk.header.system_id == info.data["system_id"]
    summary = read_scattering_summary(directory / "scattering_summary.tsv")
    assert len(summary["E_coll"]) == len(info.data["analysis_energies"]["collision_ev"])
    states = read_gas_phase_states(directory / "gas_phase_states.dat")
    assert states.wave_functions.shape[1] == info.data["grids"]["nr"]


def checkpoint_norm(checkpoint: Checkpoint, info_data: dict, with_dr: bool = True) -> float:
    """The norm of the checkpoint's wave function, as averages.dat computes it: sum |psi|^2 dZ dr dxy on the
    scattering grid plus sum |chi|^2 dZ dr on the specular grid. with_dr=False leaves dr out (the slip of the
    negative test)."""
    g = info_data["grids"]
    dr = g["dr"] if with_dr else 1.0
    scattering = float(np.sum(np.abs(checkpoint.scattering) ** 2)) * g["dz"] * dr * g["xy_volume_element"]
    specular = float(np.sum(np.abs(checkpoint.specular) ** 2)) * g["dz"] * dr
    return scattering + specular


@pytest.mark.parametrize("run", RUNS)
def test_checkpoint_norm_is_the_last_averages_norm(fixture_runs: Path, run: str) -> None:
    """The checkpoint at tEnd and the last averages row (also at tEnd) have the same norm to 1e-12 relative
    (the two sums run in different orders; the program writes averages.dat with 16 digits)."""
    directory = fixture_runs / run
    info = read_run_info(directory)
    checkpoint = read_checkpoint(directory / "checkpoint.bin")
    averages = read_averages(directory / "averages.dat")
    assert checkpoint.time_au == info.t_end_au == averages.times_au[-1]
    relative = abs(checkpoint_norm(checkpoint, info.data) - averages.norm[-1]) / averages.norm[-1]
    assert relative <= 1e-12


def test_norm_without_dr_fails(fixture_runs: Path) -> None:
    """Negative test: the volume element without dr (a realistic slip) fails the 1e-12 of the check above."""
    directory = fixture_runs / "hardwall_tiny"
    info = read_run_info(directory)
    norm = checkpoint_norm(read_checkpoint(directory / "checkpoint.bin"), info.data, with_dr=False)
    expected = read_averages(directory / "averages.dat").norm[-1]
    assert abs(norm - expected) / expected > 1e-12


def test_fine_grid_of_a_run(fixture_runs: Path) -> None:
    """grid.json of dailight_oblique's fine analysis: one spectrum row per grid point first..last, and every
    analysis energy either interpolated or listed as skipped."""
    fine = fixture_runs / "dailight_oblique" / "fine_analysis"
    grid = read_fine_grid(fine / "grid.json")
    spectrum = read_scattering_summary(fine / "spectrum" / "scattering_summary.tsv")
    assert len(spectrum["E_z"]) == grid.last_index - grid.first_index + 1
    interpolated = read_scattering_summary(fine / "scattering_summary.tsv")
    info = read_run_info(fixture_runs / "dailight_oblique")
    assert len(interpolated["E_coll"]) + len(grid.skipped) == len(info.data["analysis_energies"]["collision_ev"])
