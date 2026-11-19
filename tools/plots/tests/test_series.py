"""compute/series.py on synthetic series: the records before a NaN or inf, and records evenly spaced."""

from pathlib import Path

import numpy as np
import pytest

from compute.series import evenly_spaced_records, usable_records
from readers.arrays import Series
from readers.errors import ReadError


def series_of(values: list[float], first_non_finite: int | None) -> Series:
    """A one-dataset series with a record every 100 a.u."""
    times = 100.0 * np.arange(len(values))
    return Series(Path("s"), {}, times, {"x": np.array(values)}, first_non_finite)


def test_usable_records() -> None:
    """All records when every one is finite; the records before the first non-finite one, with a note; a
    ReadError when none is left."""
    finite = series_of([1.0, 2.0, 3.0], None)
    same, notes = usable_records(finite)
    assert same is finite and notes == ()
    cut, notes = usable_records(series_of([1.0, 2.0, np.nan], 2))
    assert cut.arrays["x"].tolist() == [1.0, 2.0] and cut.times_au.tolist() == [0.0, 100.0]
    assert cut.first_non_finite is None and "from t = 200 a.u." in notes[0]
    for empty in (series_of([], None), series_of([np.nan], 0)):
        with pytest.raises(ReadError, match="no finite record"):
            usable_records(empty)


def test_evenly_spaced_records() -> None:
    """At most the asked number, the first and the last always; all records when there are fewer; at least 2."""
    assert evenly_spaced_records(3, 30).tolist() == [0, 1, 2]
    assert evenly_spaced_records(101, 5).tolist() == [0, 25, 50, 75, 100]
    assert evenly_spaced_records(10, 2).tolist() == [0, 9]
    assert evenly_spaced_records(0, 2).tolist() == []
    with pytest.raises(ValueError, match="at least 2"):
        evenly_spaced_records(10, 1)
