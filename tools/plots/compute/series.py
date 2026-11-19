"""Helpers for any series directory (readers/arrays.py): the records that can be drawn, and records chosen
evenly for a movie."""

from dataclasses import replace

import numpy as np

from readers.arrays import Series
from readers.errors import ReadError


def usable_records(series: Series) -> tuple[Series, tuple[str, ...]]:
    """The records before the first one with a NaN or inf, and a note when some were left out. ReadError when no
    record is left (the run stopped before its first record, or diverged at once)."""
    n = len(series.times_au) if series.first_non_finite is None else series.first_non_finite
    if n == 0:
        raise ReadError(f"{series.directory}: no finite record; the run stopped before its first record or diverged")
    if n == len(series.times_au):
        return series, ()
    arrays = {name: values[:n] for name, values in series.arrays.items()}
    cut = replace(series, times_au=series.times_au[:n], arrays=arrays, first_non_finite=None)
    first_bad = f"t = {series.times_au[n]:g} a.u." if np.isfinite(series.times_au[n]) else f"record {n}"
    return cut, (f"the records from {first_bad} on hold a NaN or inf (the propagation diverged): not drawn",)


def evenly_spaced_records(n_records: int, at_most: int) -> np.ndarray:
    """At most `at_most` (>= 2) record indices, evenly spaced, the first and the last included."""
    if at_most < 2:
        raise ValueError(f"at_most {at_most}: give at least 2 (the first and the last record)")
    if n_records <= at_most:
        return np.arange(n_records)
    return np.unique(np.round(np.linspace(0, n_records - 1, at_most)).astype(int))
