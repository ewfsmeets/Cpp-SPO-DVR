"""The readers on small synthetic files: damaged and empty files, wrong headers, NaN data, segment joins."""

import json
import struct
from pathlib import Path

import numpy as np
import pytest

from copies.real_compare import rne
from readers.binary import BK_MAGIC, CHECKPOINT_MAGIC, FORMAT_VERSION, read_bk_coefficients, read_checkpoint
from readers.errors import ReadError
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
