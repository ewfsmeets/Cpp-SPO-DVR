"""pytest set-up of the Python tests of tools/fortran_reference/: the scratch directory `work` of
each test."""

from pathlib import Path

import pytest


@pytest.fixture
def work(tmp_path: Path) -> Path:
    """The test's scratch directory (pytest's tmp_path, under --basetemp)."""
    return tmp_path
