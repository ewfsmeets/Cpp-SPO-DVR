"""The readers on small synthetic files: damaged and empty files, wrong headers, NaN data, segment joins."""

import json
import struct
from pathlib import Path

import numpy as np
import pytest

from copies.real_compare import rne
from readers.arrays import read_series
from readers.binary import BK_MAGIC, CHECKPOINT_MAGIC, FORMAT_VERSION, read_bk_coefficients, read_checkpoint
from readers.errors import ReadError
from readers.fine_analysis import read_fine_grid
from readers.run_info import read_run_info
from readers.segments import join_records
from readers.tables import read_averages, read_gas_phase_states, read_table

HEADER = "# Wave packet averages (Fortran: averages.out)\n# t(a.u.)  norm  <Z>(bohr)  <r>(bohr)  rest_norm\n"


def test_damaged_row_is_refused(tmp_path: Path) -> None:
    """averages.dat is written once, completely: a row with too few numbers (or a non-number) means a damaged
    file, and the error names its line."""
    path = tmp_path / "averages.dat"
    path.write_text(HEADER + "0 1 11.4 1.47 1\n100 0.99 11.2 1.47\n200 0.98 11.0 1.47 0.98\n", encoding="utf-8")
    with pytest.raises(ReadError, match="line 4: 4 fields, not 5"):
        read_averages(path)
    path.write_text(HEADER + "0 1 11.4 1.47 x\n", encoding="utf-8")
    with pytest.raises(ReadError, match="line 3"):
        read_averages(path)


def test_file_without_rows_reads_as_empty(tmp_path: Path) -> None:
    """An averages.dat with only its '#' lines: no rows, no error."""
    path = tmp_path / "averages.dat"
    path.write_text(HEADER, encoding="utf-8")
    averages = read_averages(path)
    assert len(averages.times_au) == 0 and averages.first_non_finite is None


def test_nan_is_kept_and_recorded(tmp_path: Path) -> None:
    """A NaN stays in the data and the first row holding one is recorded (a check decides, not the reader)."""
    path = tmp_path / "averages.dat"
    path.write_text(HEADER + "0 1 11.4 1.47 1\n100 nan 11.2 1.47 0.99\n200 0.98 11.0 1.47 0.98\n", encoding="utf-8")
    averages = read_averages(path)
    assert np.isnan(averages.norm[1]) and averages.first_non_finite == 1


def test_table_without_a_required_column_or_with_a_short_row_is_refused(tmp_path: Path) -> None:
    """A result table without a needed column, or with a row missing a cell: ReadError naming the file."""
    path = tmp_path / "scattering_summary.tsv"
    path.write_text("# comment\nE_coll\tP\n0.2\t0.5\n", encoding="utf-8")
    with pytest.raises(ReadError, match="no column P_scattered"):
        read_table(path, ("E_coll", "P_scattered"))
    path.write_text("# comment\nE_coll\tP_scattered\n0.2\t0.5\n0.3\t\n", encoding="utf-8")
    with pytest.raises(ReadError, match="line 4"):
        read_table(path, ("E_coll", "P_scattered"))


def test_missing_file_is_refused(tmp_path: Path) -> None:
    """A file that does not exist: ReadError saying what to give instead."""
    with pytest.raises(ReadError, match="no such file"):
        read_averages(tmp_path / "averages.dat")


def test_gas_phase_file_needs_its_r_grid_line(tmp_path: Path) -> None:
    """gas_phase_states.dat without its '# r grid:' line is from another program version."""
    path = tmp_path / "gas_phase_states.dat"
    path.write_text("# Gas-phase rovibrational states\n0 0 -0.16 0.1 0.2\n", encoding="utf-8")
    with pytest.raises(ReadError, match="0 '# r grid:' lines, not 1"):
        read_gas_phase_states(path)


def header(magic: bytes, version: int = FORMAT_VERSION) -> bytes:
    """The common header with SystemID "Test" and no parameters."""
    return magic + struct.pack("<I", version) + struct.pack("<Q", 4) + b"Test" + struct.pack("<Q", 0)


CHANNEL_ENERGY_HARTREE = 0.01  # of the single analysed state


def bk_file(path: Path, records: int, extra_bytes: int = 0) -> None:
    """A BK file of 1 state, 1 level and `records` records (coefficient = record number), then extra bytes."""
    data = header(BK_MAGIC) + struct.pack("<2i", 1, 1) + struct.pack("<5i", 0, 0, 0, 0, 0)
    data += struct.pack("<d", CHANNEL_ENERGY_HARTREE)
    for k in range(records):
        data += struct.pack("<d", 10.0 * (k + 1)) + struct.pack("<2d", float(k), 0.0)
    path.write_bytes(data + b"\0" * extra_bytes)


def test_bk_record_cut_short_is_left_out(tmp_path: Path) -> None:
    """A BK file whose last record is incomplete (a run still writing): the complete records are read."""
    path = tmp_path / "bk_coefficients.bin"
    bk_file(path, 3, extra_bytes=11)
    bk = read_bk_coefficients(path)
    assert list(bk.times_au) == [10.0, 20.0, 30.0] and bk.coefficients[2, 0, 0] == 2.0
    assert bk.channel_energy_hartree[0, 0] == CHANNEL_ENERGY_HARTREE


def test_bk_file_without_records(tmp_path: Path) -> None:
    """A BK file with only its header (a run before its first analysis step): no records, no error."""
    path = tmp_path / "bk_coefficients.bin"
    bk_file(path, 0)
    bk = read_bk_coefficients(path)
    assert bk.coefficients.shape == (0, 1, 1) and bk.first_non_finite is None


def test_wrong_binary_files_are_refused(tmp_path: Path) -> None:
    """Another file kind, another format version, a header cut short, a checkpoint without its waves."""
    path = tmp_path / "bk_coefficients.bin"
    path.write_bytes(header(CHECKPOINT_MAGIC))
    with pytest.raises(ReadError, match="not a SPODVR-BKCOEFFS file"):
        read_bk_coefficients(path)
    path.write_bytes(header(BK_MAGIC, version=FORMAT_VERSION + 1))
    with pytest.raises(ReadError, match=f"format version {FORMAT_VERSION + 1}"):
        read_bk_coefficients(path)
    path.write_bytes(header(BK_MAGIC)[:30])
    with pytest.raises(ReadError, match="ends inside its header"):
        read_bk_coefficients(path)
    checkpoint = tmp_path / "checkpoint.bin"
    sizes = struct.pack("<d5i2iBd", 100.0, 1, 1, 1, 1, 1, 0, 0, 1, 1e-6)
    checkpoint.write_bytes(header(CHECKPOINT_MAGIC) + sizes)
    with pytest.raises(ReadError, match="the file is incomplete or damaged"):
        read_checkpoint(checkpoint)


def test_damaged_run_info_is_refused(tmp_path: Path) -> None:
    """run_info.json that is not JSON, of another format version, or without a field the tool needs."""
    path = tmp_path / "run_info.json"
    path.write_text("{", encoding="utf-8")
    with pytest.raises(ReadError, match="not JSON"):
        read_run_info(tmp_path)
    path.write_text(json.dumps({"format_version": 99}), encoding="utf-8")
    with pytest.raises(ReadError, match="format_version 99"):
        read_run_info(tmp_path)
    path.write_text(json.dumps({"format_version": 1, "times": {}}), encoding="utf-8")
    with pytest.raises(ReadError, match="no field"):
        read_run_info(tmp_path)


def test_rne_is_the_programs() -> None:
    """Equal below 1e-12 apart, different at 1e-12 and for NaN (the C++ !(|a - b| < 1e-12))."""
    assert not rne(1.0, 1.0 + 0.5e-12) and rne(0.0, 1e-12) and rne(float("nan"), float("nan"))


def test_join_refuses_a_shared_record_that_differs() -> None:
    """Two segments whose shared record differs are not one run."""
    with pytest.raises(ReadError, match="differs between two segments"):
        join_records([np.array([0.0, 1.0]), np.array([1.0, 2.0])], [np.array([1.0, 2.0]), np.array([2.5, 3.0])], "x")


def test_join_drops_the_repeated_record() -> None:
    """The record at the shared time appears once."""
    times, values = join_records(
        [np.array([0.0, 1.0]), np.array([1.0, 2.0])], [np.array([1.0, 2.0]), np.array([2.0, 3.0])], "x"
    )
    assert list(times) == [0.0, 1.0, 2.0] and list(values) == [1.0, 2.0, 3.0]


def test_fine_grid(tmp_path: Path) -> None:
    """grid.json with its fields and a skipped energy; a missing field is a ReadError."""
    fields = {
        "spectral_resolution_eV": 0.006,
        "first_index": 33,
        "last_index": 97,
        "effective_time": 13860.0,
        "records_used": 1386,
        "skipped_energies": [{"E_z": 0.2, "reason": "outside the fine-analysis grid"}],
    }
    path = tmp_path / "grid.json"
    path.write_text(json.dumps(fields), encoding="utf-8")
    grid = read_fine_grid(path)
    assert (grid.first_index, grid.last_index, grid.skipped) == (33, 97, [(0.2, "outside the fine-analysis grid")])
    del fields["first_index"]
    path.write_text(json.dumps(fields), encoding="utf-8")
    with pytest.raises(ReadError, match="first_index"):
        read_fine_grid(path)


def series_directory(path: Path, times: list[float], values: list[float], extra_bytes: int = 0) -> Path:
    """A series directory with a scalar <f8 dataset "x" and a scalar <i4 dataset "k" (k = the record's index),
    and `extra_bytes` of x from a record being written."""
    path.mkdir(parents=True)
    datasets = {
        "x": {"file": "x.bin", "dtype": "<f8", "shape": []},
        "k": {"file": "k.bin", "dtype": "<i4", "shape": []},
    }
    (path / "manifest.json").write_text(json.dumps({"format_version": 1, "datasets": datasets}), encoding="utf-8")
    (path / "times.bin").write_bytes(np.array(times, "<f8").tobytes())
    (path / "x.bin").write_bytes(np.array(values, "<f8").tobytes() + b"\0" * extra_bytes)
    (path / "k.bin").write_bytes(np.arange(len(values), dtype="<i4").tobytes())
    return path


def test_series_reads_only_complete_records(tmp_path: Path) -> None:
    """Records a program still writing leaves incomplete are left out: (a) the data of a third record and part of
    a fourth written, but not their times; (b) a third time without its data (the program writes data first,
    but the reader does not rely on it)."""
    series = read_series(series_directory(tmp_path / "a", [0.0, 100.0], [1.0, 2.0, 3.0], extra_bytes=4))
    assert list(series.times_au) == [0.0, 100.0] and list(series.arrays["x"]) == [1.0, 2.0]
    series = read_series(series_directory(tmp_path / "b", [0.0, 100.0, 200.0], [1.0, 2.0]))
    assert len(series.times_au) == 2 and series.arrays["k"].dtype == np.dtype("<i4")


def test_series_without_records(tmp_path: Path) -> None:
    """A run stopped before its first record: zero records, each dataset of its own dtype, nothing non-finite."""
    series = read_series(series_directory(tmp_path / "empty", [], []))
    assert (
        len(series.times_au) == 0 and series.arrays["x"].shape == (0,) and series.arrays["k"].dtype == np.dtype("<i4")
    )
    assert series.first_non_finite is None


def test_series_records_the_first_non_finite_record(tmp_path: Path) -> None:
    """A NaN (or inf) in a dataset marks its record; the data are returned as written."""
    series = read_series(series_directory(tmp_path / "nan", [0.0, 100.0, 200.0], [1.0, np.nan, np.inf]))
    assert series.first_non_finite == 1 and np.isnan(series.arrays["x"][1])
    assert read_series(series_directory(tmp_path / "finite", [0.0, 100.0], [1.0, 2.0])).first_non_finite is None


def test_records_of_segments_off_the_record_grid_concatenate() -> None:
    """Segments whose boundary falls between records share no record: join_records concatenates them."""
    times, values = join_records(
        [np.array([0.0, 100.0]), np.array([150.0, 250.0])], [np.array([1.0, 2.0]), np.array([3.0, 4.0])], "x"
    )
    assert list(times) == [0.0, 100.0, 150.0, 250.0] and list(values) == [1.0, 2.0, 3.0, 4.0]
