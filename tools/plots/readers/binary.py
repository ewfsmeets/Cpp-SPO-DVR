"""The binary files of a run: bk_coefficients.bin and checkpoint.bin (little-endian). Both start with the same
header (magic, format version, SystemID, the restart parameters); the layouts are in docs/output_files.md,
"Checkpoint files", and each struct format below names the fields it reads.

Shadows: src/io/checkpoint.cpp (write_header, write_checkpoint, BkFileWriter::BkFileWriter, BkFileWriter::append)
"""

import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.errors import ReadError

BK_MAGIC = b"SPODVR-BKCOEFFS\0"
CHECKPOINT_MAGIC = b"SPODVR-CHECKPNT\0"
FORMAT_VERSION = 2
END_MARKER = 0x454E444F464350  # the last uint64 of a completely written checkpoint (checkpoint.cpp end_marker)
HEADER_BYTES_AT_MOST = 1 << 20  # the header holds at most 1000 named parameters (checkpoint.cpp max_stored_parameters)

VERSION = "<I"  # format version, after the 16-byte magic
LENGTH = "<Q"  # the length of a string, or the number of parameters
VALUE = "<d"  # a parameter's value
BK_SIZES = "<2i"  # n_states, n_v
LABEL_COLUMNS = ("j", "mJ", "n", "m", "order")  # int32 each, per analysed state
CHECKPOINT_FIELDS = "<d5i2iBd"  # time; n_rot, n_xy, n_zr, n_r, n_zsp; surface index, previous index, fixed, anlyTol


@dataclass
class Header:
    """The header both files start with."""

    system_id: str
    parameters: dict[str, float]  # the restart parameters, by name
    size: int  # bytes


def unpack(fields: str, data: bytes, offset: int, path: Path) -> tuple:
    """struct.unpack_from, with a ReadError when the file ends too early."""
    try:
        return struct.unpack_from(fields, data, offset)
    except struct.error:
        raise ReadError(f"{path}: ends inside its header; the file is incomplete or not spodvr's") from None


def read_header(data: bytes, magic: bytes, path: Path) -> Header:
    """Magic, format version, SystemID and the restart parameters (a uint64 length before each string)."""
    if data[:16] != magic:
        raise ReadError(f"{path}: not a {magic[:15].decode()} file; give the file spodvr wrote")
    (version,) = unpack(VERSION, data, 16, path)
    if version != FORMAT_VERSION:
        raise ReadError(f"{path}: format version {version}, not {FORMAT_VERSION}; read it with its program version")
    offset = 16 + struct.calcsize(VERSION)

    def string() -> str:
        """The next length-prefixed string (uint64 length, then the bytes) of the header."""
        nonlocal offset
        (length,) = unpack(LENGTH, data, offset, path)
        start = offset + struct.calcsize(LENGTH)
        offset = start + length
        try:
            return data[start:offset].decode()
        except UnicodeDecodeError:
            raise ReadError(f"{path}: a header string is not text; the file is damaged or not spodvr's") from None

    system_id = string()
    (count,) = unpack(LENGTH, data, offset, path)
    offset += struct.calcsize(LENGTH)
    parameters = {}
    for _ in range(count):
        name = string()
        (parameters[name],) = unpack(VALUE, data, offset, path)
        offset += struct.calcsize(VALUE)
    return Header(system_id, parameters, offset)


@dataclass
class BkCoefficients:
    """bk_coefficients.bin: the Balint-Kurti coefficients of every analysed state, per analysis step."""

    header: Header
    labels: np.ndarray  # (n_states, 5) int: the LABEL_COLUMNS j, mJ, n, m, diffraction order of each state
    channel_energy_hartree: np.ndarray  # (n_v, n_states)
    times_au: np.ndarray  # (n_times,)
    coefficients: np.ndarray  # (n_times, n_v, n_states) complex, in bohr^-1/2
    first_non_finite: int | None  # the first record with a NaN or inf


def read_bk_coefficients(path: Path) -> BkCoefficients:
    """bk_coefficients.bin. The program appends one record per analysis step, so a run that is still going
    (or stopped) leaves a last record cut short: it is left out, and a file without records has n_times = 0.

    Source: bk_coefficients.bin (docs/output_files.md, "Checkpoint files").
    Returns: BkCoefficients; the coefficient of analysed state i_state and level v at record t is
    coefficients[t, v, i_state].
    """
    if not path.is_file():
        raise ReadError(f"{path}: no such file; give the directory of a propagating run")
    data = path.read_bytes()
    header = read_header(data, BK_MAGIC, path)
    n_states, n_v = unpack(BK_SIZES, data, header.size, path)
    offset = header.size + struct.calcsize(BK_SIZES)
    n_labels = len(LABEL_COLUMNS) * n_states
    if len(data) < offset + 4 * n_labels + 8 * n_states * n_v:
        raise ReadError(f"{path}: ends inside its state list; the file is incomplete or not spodvr's")
    labels = np.frombuffer(data, dtype="<i4", count=n_labels, offset=offset).reshape(n_states, len(LABEL_COLUMNS))
    offset += 4 * n_labels
    energies = np.frombuffer(data, dtype="<f8", count=n_states * n_v, offset=offset).reshape(n_v, n_states)
    offset += 8 * n_states * n_v
    record_bytes = 8 + 16 * n_states * n_v  # time (double), then the coefficients (complex double)
    n_times = (len(data) - offset) // record_bytes
    records = np.frombuffer(data, dtype=np.uint8, count=n_times * record_bytes, offset=offset)
    records = records.reshape(n_times, record_bytes)
    times = records[:, :8].copy().view("<f8").reshape(n_times)
    values = records[:, 8:].copy().view("<c16").reshape(n_times, n_v, n_states)
    bad = ~(np.isfinite(times) & np.isfinite(values.reshape(n_times, n_v * n_states)).all(axis=1))
    first_bad = int(np.argmax(bad)) if bad.any() else None
    return BkCoefficients(header, labels.copy(), energies.copy(), times, values, first_bad)


@dataclass
class Checkpoint:
    """checkpoint.bin: the wave packet at the end of a run and the state of the analysis surface."""

    header: Header
    time_au: float
    surface_index: int  # 0-based Z index of the specular analysis surface
    surface_fixed: bool
    specular: np.ndarray  # (n_zsp, n_r) complex: chi(Z, r), in bohr^-1
    scattering: np.ndarray  # (n_zr, n_xy, n_rot) complex: psi(izr = ir + n_r * iz, ixy, irot)
    n_r: int


def read_checkpoint(path: Path) -> Checkpoint:
    """checkpoint.bin, the wave function memory-mapped (it can be as large as psi in memory).

    Source: checkpoint.bin (docs/output_files.md, "Checkpoint files"). The scattering wave keeps the
    program's memory order (irot fastest).
    """
    if not path.is_file():
        raise ReadError(f"{path}: no such file; give the directory of a propagating run")
    with path.open("rb") as f:
        head = f.read(HEADER_BYTES_AT_MOST)
    header = read_header(head, CHECKPOINT_MAGIC, path)
    fields = unpack(CHECKPOINT_FIELDS, head, header.size, path)
    time, (n_rot, n_xy, n_zr, n_r, n_zsp), (index, _previous, fixed) = fields[0], fields[1:6], fields[6:9]
    offset = header.size + struct.calcsize(CHECKPOINT_FIELDS)
    expected = offset + 16 * (n_r * n_zsp + n_rot * n_xy * n_zr) + 8  # the waves, then the end marker
    if path.stat().st_size != expected:
        raise ReadError(f"{path}: {path.stat().st_size} bytes, not {expected}; the file is incomplete or damaged")
    if last_word(path) != END_MARKER:
        raise ReadError(f"{path}: no end marker; the checkpoint was not written completely")
    specular = np.memmap(path, dtype="<c16", mode="r", offset=offset, shape=(n_zsp, n_r))
    scattering = np.memmap(path, dtype="<c16", mode="r", offset=offset + 16 * n_r * n_zsp, shape=(n_zr, n_xy, n_rot))
    return Checkpoint(header, time, index, bool(fixed), specular, scattering, n_r)


def last_word(path: Path) -> int:
    """The last uint64 of a file, without reading the rest."""
    with path.open("rb") as f:
        f.seek(-8, 2)
        return struct.unpack("<Q", f.read(8))[0]
