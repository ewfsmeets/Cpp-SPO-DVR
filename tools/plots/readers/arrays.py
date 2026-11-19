"""Series directories (densities/ and the later plot data): manifest.json, times.bin and one record-major file
per dataset (docs/output_files.md, "Series directories"). The program appends a record's data before its
time, so the complete records are n = min(len(times), file size / record bytes) of every dataset.

Shadows: src/io/series_writer.cpp (SeriesWriter::SeriesWriter, SeriesWriter::append)
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.errors import ReadError

SERIES_FORMAT_VERSION = 1  # manifest.json's format_version (io/series_writer.cpp)


@dataclass
class Series:
    """A series directory: its manifest, the record times and every dataset as an array (n_records first)."""

    directory: Path
    manifest: dict
    times_au: np.ndarray  # (n_records,)
    arrays: dict[str, np.ndarray]  # dataset name -> (n_records, *shape of one record)
    first_non_finite: int | None  # the first record with a NaN or inf in its time or any dataset


def read_manifest(directory: Path) -> dict:
    """manifest.json of a series directory; ReadError when it is missing, not JSON or of another version."""
    path = directory / "manifest.json"
    if not path.is_file():
        raise ReadError(f"{path}: no such file; give a run written with the key that writes {directory.name}/")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ReadError(f"{path}: not JSON ({error}); the file is damaged") from None
    if not isinstance(manifest, dict):
        raise ReadError(f"{path}: not a JSON object; the file is damaged")
    if manifest.get("format_version") != SERIES_FORMAT_VERSION:
        raise ReadError(
            f"{path}: format_version {manifest.get('format_version')}, not {SERIES_FORMAT_VERSION}; update the tool"
        )
    return manifest


def record_bytes(spec: dict) -> int:
    """The bytes of one record of a dataset of the manifest (its dtype times the elements of its shape)."""
    return np.dtype(spec["dtype"]).itemsize * int(np.prod(spec["shape"], dtype=int))


def complete_records(directory: Path, manifest: dict) -> int:
    """The number of records whose time and data are all written."""
    counts = [(directory / "times.bin").stat().st_size // 8]
    for spec in manifest["datasets"].values():
        counts.append((directory / spec["file"]).stat().st_size // record_bytes(spec))
    return min(counts)


def first_non_finite_record(times_au: np.ndarray, arrays: dict[str, np.ndarray]) -> int | None:
    """The first record with a NaN or inf in its time or any dataset, or None; reads record by record and stops
    there (a finite series is read once)."""
    for k in range(len(times_au)):
        if not np.isfinite(times_au[k]) or not all(np.isfinite(values[k]).all() for values in arrays.values()):
            return k
    return None


def read_series(directory: Path) -> Series:
    """A series directory, the datasets memory-mapped (the check for non-finite values reads them once, record
    by record; later only the records in use are paged in).

    Source: <run>/<series>/ (docs/output_files.md, "Series directories").
    """
    manifest = read_manifest(directory)
    try:
        n = complete_records(directory, manifest)
    except (FileNotFoundError, KeyError, TypeError) as error:
        raise ReadError(
            f"{directory}: a dataset file or field is missing or wrong ({error}); the series is damaged"
        ) from None
    except ZeroDivisionError:
        raise ReadError(f"{directory}: a dataset of zero bytes per record; the series is damaged") from None
    times = np.fromfile(directory / "times.bin", dtype="<f8", count=n)
    arrays = {}
    for name, spec in manifest["datasets"].items():
        shape = (n, *spec["shape"])
        if n:
            arrays[name] = np.memmap(directory / spec["file"], dtype=spec["dtype"], mode="r", shape=shape)
        else:
            arrays[name] = np.zeros(shape, dtype=spec["dtype"])
    return Series(directory, manifest, times, arrays, first_non_finite_record(times, arrays))
