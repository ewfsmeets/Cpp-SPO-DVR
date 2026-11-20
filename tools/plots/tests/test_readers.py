"""The readers on small synthetic files: damaged and empty files, wrong headers, NaN data, segment joins."""

import json
import struct
from pathlib import Path

import numpy as np
import pytest

from copies.real_compare import rne
from readers.arrays import read_arrays, read_series
from readers.binary import BK_MAGIC, CHECKPOINT_MAGIC, FORMAT_VERSION, read_bk_coefficients, read_checkpoint
from readers.errors import ReadError
from readers.fine_analysis import read_fine_grid
from readers.pes_cut import read_pes_cut
from readers.potential import read_potential
from readers.run_info import read_run_info
from readers.segments import join_records, joined_first_non_finite
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


TOLERANCE_CASES = [  # the repeated record [0.25, 4.0] changed: (new record, accepted)
    ([0.25, 4.0 * (1 + 1e-14)], True),  # rounding, as on another device
    ([0.25 + 2e-12, 4.0], True),  # 5e-13 of the largest value, in the small element (8e-12 of that element itself)
    ([0.25 + 8e-12, 4.0], False),  # 2e-12 of the largest value
    ([np.nan, 4.0], False),
]


@pytest.mark.parametrize(("repeated", "accepted"), TOLERANCE_CASES)
def test_join_compares_a_repeated_record_to_its_largest_value(repeated: list[float], accepted: bool) -> None:
    """A repeated record must equal the kept one to 1e-12 of the record's largest |value| (not of each element):
    then it is dropped and the kept one stays; else, or with a NaN, the join is refused."""
    times = [np.array([0.0, 1.0]), np.array([1.0, 2.0])]
    earlier, later = np.array([[0.5, 2.0], [0.25, 4.0]]), np.array([repeated, [0.1, 5.0]])
    if accepted:
        _, joined = join_records(times, [earlier, later], "x")
        assert joined.tolist() == [[0.5, 2.0], [0.25, 4.0], [0.1, 5.0]]
    else:
        with pytest.raises(ReadError, match="differs between two segments|NaN or inf"):
            join_records(times, [earlier, later], "x")


def test_joined_first_non_finite_counts_the_dropped_repeats() -> None:
    """Segments of 3 and 4 records, the second's first one a dropped repeat: its record 2 is the joined record 4;
    a first segment's record comes first; none gives None."""
    assert joined_first_non_finite([None, 2], [3, 4], [0, 1]) == 4
    assert joined_first_non_finite([1, 2], [3, 4], [0, 1]) == 1
    assert joined_first_non_finite([None, None], [3, 4], [0, 1]) is None


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


def test_series_manifest_errors(tmp_path: Path) -> None:
    """A missing manifest, one that is not JSON or not an object, another format version, a missing dataset file
    and a dataset of zero bytes per record: each a ReadError with its own message."""
    path = series_directory(tmp_path / "s", [0.0], [1.0])
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    cases = {
        "no such file": None,
        "not JSON": "{",
        "not a JSON object": "[1, 2]",
        "format_version 2": json.dumps({**manifest, "format_version": 2}),
        "missing or wrong": json.dumps({**manifest, "datasets": {"y": {"file": "y.bin", "dtype": "<f8", "shape": []}}}),
        "zero bytes per record": json.dumps(
            {**manifest, "datasets": {"x": {"file": "x.bin", "dtype": "<f8", "shape": [0]}}}
        ),
    }
    for message, text in cases.items():
        (path / "manifest.json").unlink(missing_ok=True)
        if text is not None:
            (path / "manifest.json").write_text(text, encoding="utf-8")
        with pytest.raises(ReadError, match=message):
            read_series(path)


def test_arrays_written_once(tmp_path: Path) -> None:
    """An array directory (time_dependent false): each dataset one array of its shape; a file that is cut, a time
    series and a missing file are each a ReadError."""
    path = tmp_path / "arrays"
    path.mkdir()
    datasets = {"v": {"file": "v.bin", "dtype": "<f8", "shape": [2, 3]}}
    manifest = {"format_version": 1, "time_dependent": False, "datasets": datasets}
    (path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (path / "v.bin").write_bytes(np.arange(6, dtype="<f8").tobytes())
    assert read_arrays(path).arrays["v"].tolist() == [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]
    (path / "v.bin").write_bytes(np.arange(5, dtype="<f8").tobytes())
    with pytest.raises(ReadError, match="40 bytes, not 48"):
        read_arrays(path)
    (path / "manifest.json").write_text(json.dumps({**manifest, "time_dependent": True}), encoding="utf-8")
    with pytest.raises(ReadError, match="time_dependent is not false"):
        read_arrays(path)
    (path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (path / "v.bin").unlink()
    with pytest.raises(ReadError, match="no such file"):
        read_arrays(path)


def small_potential_manifest() -> dict:
    """The manifest of a small valid potential/: 2 Z x 2 r points, a 2 x 1 (x, y) grid, 1 x 2 nodes, one site."""
    datasets = {
        "v_min": {"file": "v_min.bin", "dtype": "<f8", "shape": [2, 2, 2, 1], "axes": ["Z", "r", "x", "y"]},
        "v_min_node": {
            "file": "v_min_node.bin",
            "dtype": "<i4",
            "shape": [2, 2, 2, 1, 2],
            "axes": ["Z", "r", "x", "y", "node"],
        },
        "v_sites": {
            "file": "v_sites.bin",
            "dtype": "<f8",
            "shape": [1, 2, 2, 1, 2],
            "axes": ["site", "Z", "r", "theta", "phi"],
        },
    }
    limits = {"use_switching": False, "z_min_switch": 6.0, "z_max_switch": 8.0, "impose_v_min": False}
    limits |= {"v_min_hartree": -0.2, "impose_v_max": False, "v_max_hartree": 0.4}
    return {
        "format_version": 1,
        "time_dependent": False,
        "datasets": datasets,
        "grids": {"z_start": 0.0, "dz": 1.0, "nz_potential": 2, "r_start": 1.0, "dr": 0.5, "nr": 2, "nx": 2, "ny": 1},
        "theta_nodes_rad": [1.5707963267948966],
        "phi_nodes_rad": [0.0, 3.141592653589793],
        "homonuclear": True,
        "sites": [{"name": "top", "x": 0.0, "y": 0.0, "ix": 0, "iy": 0, "ixy": 0}],
        "limits": limits,
    }


def changed(manifest: dict, path: str, value: object) -> None:
    """Replaces the field at a dotted path (list items by index) by value; None deletes it."""
    *parents, last = path.split(".")
    target = manifest
    for key in parents:
        target = target[int(key)] if isinstance(target, list) else target[key]
    if value is None:
        del target[last]
    else:
        target[int(last) if isinstance(target, list) else last] = value


def potential_directory(run: Path, change: dict | None = None) -> Path:
    """small_potential_manifest's potential/ in `run`, with fields changed (changed), and zero-filled datasets."""
    directory = run / "potential"
    directory.mkdir(parents=True)
    manifest = small_potential_manifest()
    for path, value in (change or {}).items():
        changed(manifest, path, value)
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    for spec in manifest.get("datasets", {}).values():
        count = int(np.prod(spec["shape"]))
        (directory / spec["file"]).write_bytes(np.zeros(count, dtype=spec["dtype"]).tobytes())
    return run


def test_potential_reader_refuses_what_it_cannot_index(tmp_path: Path) -> None:
    """A valid potential/ is read; a dropped field, swapped axes, a wrong shape, a wrong ixy and a missing v_sites
    with sites are each a ReadError naming the manifest (an x/y mix-up would otherwise go unnoticed)."""
    potential = read_potential(potential_directory(tmp_path / "valid"))
    assert potential.v_min_hartree.shape == (2, 2, 2, 1) and potential.sites[0].name == "top"
    cases = {
        "'impose_v_max'": {"limits.impose_v_max": None},
        "axes": {"datasets.v_min.axes": ["Z", "r", "y", "x"]},
        "shape": {"grids.nx": 1},
        "ixy": {"sites.0.ixy": 1},
        "datasets": {"datasets.v_sites": None},
        "'axes'": {"datasets.v_min.axes": None},
        "'datasets'": {"datasets": None},
    }
    for k, (reason, change) in enumerate(cases.items()):
        with pytest.raises(ReadError, match=reason):
            read_potential(potential_directory(tmp_path / f"case{k}", change))


PES_CUT_HEADER = (
    "# V(r, Z) of the PES at fixed x, y, theta, phi (write2DPES)\n"
    "# x = 1.5000000000000000e+00 bohr, y = 0.0000000000000000e+00 bohr (along the lattice vectors)\n"
    "# theta = 9.0000000000000000e+01 deg (1.5707963267948966e+00 rad), phi = 4.5000000000000000e+01 deg "
    "(7.8539816339744828e-01 rad)\n"
    "# r(bohr)  Z(bohr)  V(hartree)  V(eV); Z fastest, a blank line after each r\n"
)


def pes_cut_text(rows: list[tuple[float, float, float]], block: int) -> str:
    """A pes_cut_2d.dat of rows (r, Z, V hartree) with a blank line after every `block` rows."""
    lines = []
    for k, (r, z, v) in enumerate(rows):
        lines.append(f"{r:.16e} {z:.16e} {v:.16e} {v * 27.2114:.16e}")
        if (k + 1) % block == 0:
            lines.append("")
    return PES_CUT_HEADER + "\n".join(lines) + "\n"


def test_pes_cut(tmp_path: Path) -> None:
    """A cut of 2 r x 3 Z points: V(Z, r) with Z first (V = 10 r + Z here), the grids and the header's position and
    orientation. Each a ReadError: a file cut inside its first block or a block short (the count from run_info.json),
    a cut last line, a repeated block (r not increasing), a Z grid that differs between blocks, a missing header."""
    rows = [(r, z, 10.0 * r + z) for r in (1.0, 2.0) for z in (-1.0, 0.0, 1.0)]
    (tmp_path / "pes_cut_2d.dat").write_text(pes_cut_text(rows, 3), encoding="utf-8")
    cut = read_pes_cut(tmp_path, 3, 2)
    assert cut.v_hartree.shape == (3, 2) and cut.v_hartree[2, 1] == 21.0 and cut.v_hartree[0, 0] == 9.0
    assert cut.z_bohr.tolist() == [-1.0, 0.0, 1.0] and cut.r_bohr.tolist() == [1.0, 2.0]
    assert (cut.x_bohr, cut.y_bohr, cut.theta_deg, cut.phi_deg) == (1.5, 0.0, 90.0, 45.0)
    broken = [  # (the reason the error names, the file)
        ("2 points, not nR x nZ = 2 x 3", pes_cut_text(rows[:2], 3)),
        ("5 points, not nR x nZ", pes_cut_text(rows[:-1], 3)),
        ("fields, not 4", pes_cut_text(rows, 3).rstrip("\n").rsplit(" ", 1)[0] + "\n"),
        ("increasing Z points per increasing r", pes_cut_text([*rows[:3], *rows[:3]], 3)),
        ("increasing Z points per increasing r", pes_cut_text([*rows[:5], (2.0, 1.5, 0.0)], 3)),
        ("no '# x = ...'", pes_cut_text(rows, 3).replace("# x = ", "# u = ")),
    ]
    for k, (reason, text) in enumerate(broken):
        (tmp_path / f"case{k}").mkdir()
        (tmp_path / f"case{k}" / "pes_cut_2d.dat").write_text(text, encoding="utf-8")
        with pytest.raises(ReadError, match=reason):
            read_pes_cut(tmp_path / f"case{k}", 3, 2)
