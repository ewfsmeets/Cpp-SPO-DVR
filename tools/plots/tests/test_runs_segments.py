"""Runs split into restart segments, joined: the same averages as the unsplit run, bitwise."""

from pathlib import Path

import numpy as np
import pytest

from readers.errors import ReadError
from readers.segments import ordered_segments, read_joined_averages
from readers.tables import Averages, read_averages


def columns(a: Averages) -> np.ndarray:
    """The data of averages.dat as one array (times first), for a bitwise comparison."""
    return np.column_stack([a.times_au, a.norm, a.z_mean_bohr, a.r_mean_bohr, a.rest_norm])


@pytest.mark.parametrize("name", ["hardwall_tiny", "dailight_oblique"])
def test_split_run_joins_to_the_unsplit_run(fixture_runs: Path, name: str) -> None:
    """The two segments' averages.dat, joined (the row at the split time once), equal the unsplit run's
    bitwise, and the second segment starts at the split time the fixture writer recorded."""
    split = fixture_runs / f"{name}_split"
    joined, starts = read_joined_averages([split / "second", split / "first"])  # any order
    whole = read_averages(fixture_runs / name / "averages.dat")
    assert starts == [0.0, float((split / "split_time.txt").read_text(encoding="utf-8"))]
    assert columns(joined).tobytes() == columns(whole).tobytes()


def test_segments_that_do_not_meet_are_refused(fixture_runs: Path) -> None:
    """A segment that does not start where the previous one ends (the unsplit run, then the second segment)."""
    with pytest.raises(ReadError, match="starts at 4000"):
        ordered_segments([fixture_runs / "hardwall_tiny", fixture_runs / "hardwall_tiny_split" / "second"])
