"""pytest set-up of the plot tool's tests: figures without a window (matplotlib's Agg backend, chosen before
anything imports matplotlib), and the runs that the tests named test_runs_* read: --fixture-runs, written by
the ctest PlotRuns.Fixtures (tests/test_plot_fixtures.cpp of the repository), and a copy of one as a run
that stopped early."""

import json
import os
import shutil
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

from readers.arrays import record_bytes  # noqa: E402 (after the backend is chosen)
from readers.binary import read_bk_coefficients  # noqa: E402

REPOSITORY = Path(__file__).resolve().parents[3]


def pytest_addoption(parser: pytest.Parser) -> None:
    """--fixture-runs: the directory of the fixture runs (ctest passes it; default: the build tree's)."""
    parser.addoption(
        "--fixture-runs",
        type=Path,
        default=REPOSITORY / "build" / "tests" / "plot_fixtures",
        help="the runs written by ctest PlotRuns.Fixtures",
    )


@pytest.fixture(scope="session")
def fixture_runs(request: pytest.FixtureRequest) -> Path:
    """The fixture runs; fails (never skips) when they have not been written."""
    path = request.config.getoption("--fixture-runs")
    if not (path / "hardwall_tiny" / "run_info.json").is_file():
        pytest.fail(f"no fixture runs in {path}: write them with `ctest --test-dir build -R '^PlotRuns\\.'`")
    return path


STOPPED_RECORDS = 50  # BK records kept in the stopped run (one per analysis step: up to t = 50 dtAna)


@pytest.fixture(scope="session")
def stopped_run(fixture_runs: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """hardwall_tiny as a run that stopped (or is still running) at t = 500 a.u.: only the files the program
    writes before the propagation (run_info.json, the gas-phase files) and bk_coefficients.bin cut after
    STOPPED_RECORDS records plus part of the next (the program appends one record per analysis step), and
    densities/ cut after its records up to the same time (tStart and every dtAve)."""
    whole = fixture_runs / "hardwall_tiny"
    run = tmp_path_factory.mktemp("stopped") / "hardwall_tiny_stopped"
    run.mkdir()
    for name in ("run_info.json", "gas_phase_states.dat", "asymptotic_potential.dat"):
        shutil.copy(whole / name, run / name)
    bk = read_bk_coefficients(whole / "bk_coefficients.bin")
    n_times, n_v, n_states = bk.coefficients.shape
    bk_record_bytes = 8 + 16 * n_v * n_states
    data = (whole / "bk_coefficients.bin").read_bytes()
    records_start = len(data) - n_times * bk_record_bytes
    (run / "bk_coefficients.bin").write_bytes(
        data[: records_start + STOPPED_RECORDS * bk_record_bytes + bk_record_bytes // 2]
    )
    times = json.loads((whole / "run_info.json").read_text(encoding="utf-8"))["times"]
    density_records = 1 + STOPPED_RECORDS // times["analyses_per_average"]
    copy_densities_up_to(whole / "densities", run / "densities", density_records)
    return run


def copy_densities_up_to(source: Path, target: Path, records: int) -> None:
    """densities/ as the program leaves it when it stops after `records` records: every dataset file cut there,
    the next record of the first dataset half written (its time not yet)."""
    target.mkdir()
    manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    shutil.copy(source / "manifest.json", target / "manifest.json")
    (target / "times.bin").write_bytes((source / "times.bin").read_bytes()[: 8 * records])
    for k, spec in enumerate(manifest["datasets"].values()):
        size = record_bytes(spec)
        kept = records * size + (size // 2 if k == 0 else 0)
        (target / spec["file"]).write_bytes((source / spec["file"]).read_bytes()[:kept])


@pytest.fixture(scope="session")
def fresh_run(stopped_run: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The stopped run before its first analysis step: bk_coefficients.bin holds its header, its state list and
    half a record, no complete record."""
    run = tmp_path_factory.mktemp("fresh") / "hardwall_tiny_fresh"
    shutil.copytree(stopped_run, run)
    shutil.rmtree(run / "densities")
    copy_densities_up_to(stopped_run / "densities", run / "densities", 1)  # only the record at tStart
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    n_times, n_v, n_states = bk.coefficients.shape
    bk_record_bytes = 8 + 16 * n_v * n_states
    data = (run / "bk_coefficients.bin").read_bytes()
    (run / "bk_coefficients.bin").write_bytes(data[: len(data) - n_times * bk_record_bytes])  # the half record stays
    return run


@pytest.fixture(scope="session")
def no_density_record_run(fresh_run: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The fresh run stopped before its first averages time: densities/ holds its manifest and half a record."""
    run = tmp_path_factory.mktemp("no_density_record") / "hardwall_tiny_no_density_record"
    shutil.copytree(fresh_run, run)
    shutil.rmtree(run / "densities")
    copy_densities_up_to(fresh_run / "densities", run / "densities", 0)
    return run
