"""Shared helpers of the validation scripts: paths, the results store, Fortran number parsing,
printed-precision tolerances, statistics, and readers of the binary run files.

Every number that appears in the report is computed by a script from the raw run files and
stored in `scripts/results.json` (via `update_results`) or in a small extract under `data/`.
"""
from __future__ import annotations

import json
import math
import os
import re
import struct
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
DOC = SCRIPTS.parent
REPO = DOC.parents[1]
DATA = DOC / "data"
ITEMS = DATA / "items"          # per-item extracts (JSON), read by the default make_all mode
FIGURES = DOC / "figures"
TABLES = SCRIPTS / "tables"     # generated Markdown table fragments
RUNS = REPO / "runs"
WORK = RUNS / "thesis_fig5_3_dz_validation"   # git-ignored: helper binaries, big intermediates
INPUT_DIR = REPO / "inputs" / "thesis_fig5_3_dz"

RUN_DIRS = {
    "cpu": RUNS / "thesis_fig5_3_dz_cpu_fast",
    "gpu": RUNS / "thesis_fig5_3_dz_cuda_fast",
    "fortran": RUNS / "thesis_fig5_3_dz_fortran",
}
RAW_LOGS = {
    "cpu": RUNS / "thesis_fig5_3_dz_cpu_fast.log",
    "gpu": RUNS / "thesis_fig5_3_dz_cuda_fast.log",
    "fortran": RUNS / "thesis_fig5_3_dz_fortran.log",
}
CLEAN_LOGS = {k: WORK / "logs" / f"{k}.txt" for k in RAW_LOGS}
LABEL = {"cpu": "C++ CPU", "gpu": "C++ GPU", "fortran": "Fortran"}
PAIRS = [("cpu", "fortran"), ("gpu", "fortran"), ("gpu", "cpu")]
PAIR_LABEL = {("cpu", "fortran"): "C++−F", ("gpu", "fortran"): "GPU−F", ("gpu", "cpu"): "GPU−C++"}

RESULTS = SCRIPTS / "results.json"

# Units, as used by the codes (UnitConversion legacy / precise)
CHAEV = 27.2114            # legacy hartree -> eV
EV2HA_PRECISE = 0.036749308286


# ----------------------------------------------------------------------------- results store
def load_item(name: str) -> dict:
    path = ITEMS / f"{name}.json"
    return json.loads(path.read_text()) if path.exists() else {}


def save_item(name: str, data: dict) -> None:
    ITEMS.mkdir(parents=True, exist_ok=True)
    text = json.dumps(to_jsonable(data), indent=1, sort_keys=True, allow_nan=False)
    (ITEMS / f"{name}.json").write_text(text + "\n")


def to_jsonable(x):
    if isinstance(x, dict):
        return {str(k): to_jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [to_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        x = float(x)
    if isinstance(x, float):
        if math.isnan(x) or math.isinf(x):
            return str(x)
        return x
    if isinstance(x, np.ndarray):
        return to_jsonable(x.tolist())
    if isinstance(x, Path):
        return str(x)
    return x


# ----------------------------------------------------------------------------- Fortran numbers
FORTRAN_REAL = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([EeDd][+-]?\d+|[+-]\d{3})?$")


def parse_fortran_real(text: str) -> float:
    """Fortran E/D/F text -> float. 'D' -> 'E'; 3-digit exponents drop the letter (1.0-100)."""
    t = text.strip().replace("D", "E").replace("d", "e")
    if "*" in t:
        raise ValueError(f"Fortran overflow field: {text!r}")
    m = re.match(r"^([+-]?\d*\.\d*)([+-]\d{3})$", t)
    if m:
        t = m.group(1) + "E" + m.group(2)
    return float(t)


def printed_half_unit(text: str, floor: float = 1e-20) -> float:
    """1.01 x half a unit in the last printed digit of a Fortran E/D/F number (tests/analysis_reference.hpp)."""
    t = text.strip().replace("D", "E").replace("d", "E").replace("e", "E")
    m = re.match(r"^([+-]?\d*\.\d*)([+-]\d{3})$", t)
    if m:
        t = m.group(1) + "E" + m.group(2)
    if "E" in t:
        mantissa, exponent = t.split("E")
        exponent = int(exponent)
    else:
        mantissa, exponent = t, 0
    digits = len(mantissa) - mantissa.index(".") - 1 if "." in mantissa else 0
    return max(1.01 * 0.5 * 10.0 ** (exponent - digits), floor)


# ----------------------------------------------------------------------------- statistics
def stats(values, where=None) -> dict:
    """n, median, p95, max and the location of the max of |values| (where: list of labels)."""
    a = np.abs(np.asarray(values, dtype=float))
    if a.size == 0:
        return {"n": 0}
    i = int(np.argmax(a))
    out = {
        "n": int(a.size),
        "median": float(np.median(a)),
        "p95": float(np.percentile(a, 95)),
        "max": float(a[i]),
    }
    if where is not None:
        out["argmax"] = where[i] if not callable(where) else where(i)
    return out


def fmt(x: float, digits: int = 2) -> str:
    """Short scientific format used in tables: 1.2e-13; 0 -> '0'."""
    if x == 0:
        return "0"
    return f"{x:.{digits - 1}e}".replace("e-0", "e-").replace("e+0", "e+")


# ----------------------------------------------------------------------------- logs
def clean_log_lines(name: str) -> list[str]:
    return CLEAN_LOGS[name].read_text(encoding="utf-8").splitlines()


def dotted_value(lines: list[str], label: str) -> str:
    """Value of a C++ log line '  Label ..... value' (first match)."""
    for line in lines:
        s = line.strip()
        if s.startswith(label):
            rest = s[len(label):]
            rest = rest.lstrip(" .")
            return rest.strip()
    raise KeyError(label)


# ----------------------------------------------------------------------------- gfortran unformatted
def fortran_record_offsets(path: Path, max_records: int | None = None):
    """(offset, length) of the payload of every record of a gfortran sequential unformatted file."""
    out = []
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        pos = 0
        while pos < size:
            f.seek(pos)
            (length,) = struct.unpack("<i", f.read(4))
            if length < 0:
                raise ValueError("subrecords are not supported here")
            f.seek(pos + 4 + length)
            (trailer,) = struct.unpack("<i", f.read(4))
            if trailer != length:
                raise ValueError(f"corrupt record at {pos} in {path}")
            out.append((pos + 4, length))
            pos += length + 8
            if max_records and len(out) >= max_records:
                break
    return out


def read_record(path: Path, offset: int, length: int) -> bytes:
    with open(path, "rb") as f:
        f.seek(offset)
        return f.read(length)


class FortranBk:
    """fort.21: ten header records, then one record of nstata*nviba complex per analysis step."""

    def __init__(self, path: Path):
        self.path = Path(path)
        recs = fortran_record_offsets(self.path)
        self.records = recs
        h = [read_record(self.path, o, n) for o, n in recs[:10]]
        self.n_records, self.nstatw, self.ndftow, self.nrotw = struct.unpack("<4i", h[0])
        self.deltat, self.kz, self.z0, self.alpha, self.z_analysis, self.e_initial = struct.unpack("<6d", h[1])
        e, th, ph, jin, mjin, ichn = struct.unpack("<3d3i", h[2])
        self.e_translational, self.theta, self.phi, self.j_in, self.mj_in, self.ichn = e, th, ph, jin, mjin, ichn
        n = len(h[3]) // 4
        names = ["PntToRot", "PntToDif", "n", "m", "J", "mJ"]
        self.labels = {nm: np.frombuffer(h[3 + k], dtype="<i4") for k, nm in enumerate(names)}
        self.nstata = n
        self.enymo = np.frombuffer(h[9], dtype="<f8")
        self.nviba = self.enymo.size // n
        data = recs[10:]
        assert len(data) == self.n_records, (len(data), self.n_records)
        sizes = {length for _, length in data}
        assert sizes == {16 * self.nstata * self.nviba}, sizes
        self.first_offset = data[0][0]
        self.stride = data[1][0] - data[0][0] if len(data) > 1 else 0
        # memory map: records are payload + 8 bytes of markers apart
        raw = np.memmap(self.path, dtype=np.uint8, mode="r")
        self._raw = raw
        self.count = 2 * self.nstata * self.nviba

    def times(self):
        return self.deltat * np.arange(1, self.n_records + 1)

    def record(self, k: int) -> np.ndarray:
        """Coefficients of record k (0-based), shape (nviba, nstata) [v, is]."""
        o = self.first_offset + k * self.stride
        a = np.frombuffer(self._raw[o:o + 8 * self.count], dtype="<f8")
        c = a[0::2] + 1j * a[1::2]
        return c.reshape(self.nviba, self.nstata)


class NativeBk:
    """bk_coefficients.bin (docs/output_files.md)."""

    def __init__(self, path: Path):
        self.path = Path(path)
        with open(self.path, "rb") as f:
            magic = f.read(16)
            assert magic.startswith(b"SPODVR-BKCOEFFS"), magic
            (self.version,) = struct.unpack("<I", f.read(4))
            (n,) = struct.unpack("<Q", f.read(8))
            self.system_id = f.read(n).decode()
            (count,) = struct.unpack("<Q", f.read(8))
            self.parameters = {}
            for _ in range(count):
                (n,) = struct.unpack("<Q", f.read(8))
                name = f.read(n).decode()
                (value,) = struct.unpack("<d", f.read(8))
                self.parameters[name] = value
            self.nstata, self.nviba = struct.unpack("<2i", f.read(8))
            lab = np.frombuffer(f.read(4 * 5 * self.nstata), dtype="<i4").reshape(self.nstata, 5)
            self.labels = {"J": lab[:, 0], "mJ": lab[:, 1], "n": lab[:, 2], "m": lab[:, 3], "order": lab[:, 4]}
            self.energies = np.frombuffer(f.read(8 * self.nstata * self.nviba), dtype="<f8").reshape(self.nviba, self.nstata)
            self.first_offset = f.tell()
        size = os.path.getsize(self.path)
        self.record_bytes = 8 + 16 * self.nstata * self.nviba
        n_rec, rest = divmod(size - self.first_offset, self.record_bytes)
        assert rest == 0, rest
        self.n_records = n_rec
        self._raw = np.memmap(self.path, dtype=np.uint8, mode="r")

    def time(self, k: int) -> float:
        o = self.first_offset + k * self.record_bytes
        return float(np.frombuffer(self._raw[o:o + 8], dtype="<f8")[0])

    def times(self):
        return np.array([self.time(k) for k in range(self.n_records)])

    def record(self, k: int) -> np.ndarray:
        o = self.first_offset + k * self.record_bytes + 8
        a = np.frombuffer(self._raw[o:o + 16 * self.nstata * self.nviba], dtype="<f8")
        c = a[0::2] + 1j * a[1::2]
        return c.reshape(self.nviba, self.nstata)
