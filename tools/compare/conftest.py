"""pytest set-up of the tests of tools/compare/: the scratch directory `work` of each test, and the
spodvr_compare_wave executable that the CompareWave tests run (option --compare-wave-tool, which
ctest passes; a test that needs it fails without it, it never skips)."""

from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--compare-wave-tool", type=Path, default=None, help="the spodvr_compare_wave executable")


@pytest.fixture
def work(tmp_path: Path) -> Path:
    """The test's scratch directory (pytest's tmp_path, under --basetemp)."""
    return tmp_path


@pytest.fixture
def tool(request: pytest.FixtureRequest) -> Path:
    """The spodvr_compare_wave executable of --compare-wave-tool."""
    path = request.config.getoption("--compare-wave-tool")
    if path is None:
        pytest.fail("--compare-wave-tool is missing: give build/tools/spodvr_compare_wave (ctest does)")
    if not path.is_file():
        pytest.fail(f"{path} is not a file; build the target spodvr_compare_wave")
    return path
