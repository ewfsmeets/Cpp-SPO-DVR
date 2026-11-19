"""pytest set-up of the designer's tests: figures without a window (matplotlib's Agg backend, chosen
before anything imports matplotlib; ctest also sets it and a scratch matplotlib cache), and the
scratch directory `work` of each test."""

import os
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")


@pytest.fixture
def work(tmp_path: Path) -> Path:
    """The test's scratch directory (pytest's tmp_path, under --basetemp)."""
    return tmp_path
