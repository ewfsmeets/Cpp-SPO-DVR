"""pytest set-up of the plot tool's tests: figures without a window (matplotlib's Agg backend, chosen before
anything imports matplotlib), and the runs that the tests named test_runs_* read: --fixture-runs, written by
the ctest PlotRuns.Fixtures (tests/test_plot_fixtures.cpp of the repository), and a copy of one as a run
that stopped early."""

import os
import shutil
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

from readers.binary import read_bk_coefficients  # noqa: E402 (after the backend is chosen)

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


STOPPED_RECORDS = 50  # BK records kept in the stopped run


@pytest.fixture(scope="session")
def stopped_run(fixture_runs: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """hardwall_tiny as a run that stopped (or is still running) at t = 500 a.u.: only the files the program
    writes before the propagation (run_info.json, the gas-phase files) and bk_coefficients.bin cut after
    STOPPED_RECORDS records plus part of the next (the program appends one record per analysis step)."""
    whole = fixture_runs / "hardwall_tiny"
    run = tmp_path_factory.mktemp("stopped") / "hardwall_tiny_stopped"
    run.mkdir()
    for name in ("run_info.json", "gas_phase_states.dat", "asymptotic_potential.dat"):
        shutil.copy(whole / name, run / name)
    bk = read_bk_coefficients(whole / "bk_coefficients.bin")
    n_times, n_v, n_states = bk.coefficients.shape
    record_bytes = 8 + 16 * n_v * n_states
    data = (whole / "bk_coefficients.bin").read_bytes()
    records_start = len(data) - n_times * record_bytes
    (run / "bk_coefficients.bin").write_bytes(
        data[: records_start + STOPPED_RECORDS * record_bytes + record_bytes // 2]
    )
    return run


@pytest.fixture(scope="session")
def fresh_run(stopped_run: Path, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The stopped run before its first analysis step: bk_coefficients.bin holds its header, its state list and
    half a record, no complete record."""
    run = tmp_path_factory.mktemp("fresh") / "hardwall_tiny_fresh"
    shutil.copytree(stopped_run, run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    n_times, n_v, n_states = bk.coefficients.shape
    record_bytes = 8 + 16 * n_v * n_states
    data = (run / "bk_coefficients.bin").read_bytes()
    (run / "bk_coefficients.bin").write_bytes(data[: len(data) - n_times * record_bytes])  # the half record stays
    return run
