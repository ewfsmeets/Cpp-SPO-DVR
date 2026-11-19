"""The text files of a run: averages.dat, the tab-separated result tables, the gas-phase files. The program
writes each of them once, completely (averages.dat and the tables after the propagation, the gas-phase files
at set-up), so a line that does not parse means a damaged file: a ReadError, never a guess.

Shadows: src/report/propagation_report.cpp (write_averages_file), src/report/analysis_files.cpp (write_summary,
write_s_matrix), src/report/sections.cpp (write_gas_phase_file, write_asymptotic_potential_file)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.errors import ReadError


def first_non_finite(rows: np.ndarray) -> int | None:
    """The index of the first row (along axis 0) that holds a NaN or inf, or None (also for no rows)."""
    if len(rows) == 0:
        return None
    bad = ~np.isfinite(rows.reshape(len(rows), -1)).all(axis=1)
    return int(np.argmax(bad)) if bad.any() else None


def text_of(path: Path) -> str:
    """The file's text; ReadError when it does not exist."""
    if not path.is_file():
        raise ReadError(f"{path}: no such file; give the run directory that holds it")
    return path.read_text(encoding="utf-8")


def parse_row(fields: list[str], columns: int, path: Path, line_number: int) -> list[float]:
    """The first `columns` fields of a line as floats; ReadError (naming the line) when there are fewer or one
    is not a number."""
    try:
        if len(fields) < columns:
            raise ValueError(f"{len(fields)} fields, not {columns}")
        return [float(f) for f in fields[:columns]]
    except ValueError as error:
        raise ReadError(f"{path}: line {line_number}: {error}; the file is damaged or incomplete") from None


def numbers(path: Path, text: str, columns: int) -> np.ndarray:
    """The rows of a space-separated file after its '#' lines, as a (rows, columns) float array."""
    rows = [
        parse_row(line.split(), columns, path, number)
        for number, line in enumerate(text.splitlines(), start=1)
        if line.strip() and not line.startswith("#")
    ]
    return np.array(rows, dtype=float).reshape(len(rows), columns)


@dataclass
class Averages:
    """averages.dat: one row per averages time (every dtAve)."""

    times_au: np.ndarray
    norm: np.ndarray  # total norm (dimensionless)
    z_mean_bohr: np.ndarray
    r_mean_bohr: np.ndarray
    rest_norm: np.ndarray  # norm of the specular wave between the two analysis lines
    first_non_finite: int | None  # the first row with a NaN or inf


def read_averages(path: Path) -> Averages:
    """averages.dat of a run: columns t (a.u.), norm, <Z> (bohr), <r> (bohr), rest norm.

    Source: averages.dat (docs/output_files.md, "Data files").
    """
    rows = numbers(path, text_of(path), 5)
    return Averages(rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4], first_non_finite(rows))


@dataclass
class Table:
    """A tab-separated result table: its columns by header name."""

    columns: dict[str, np.ndarray]
    first_non_finite: int | None

    def __getitem__(self, name: str) -> np.ndarray:
        """The column of that header name."""
        return self.columns[name]


def read_table(path: Path, required: tuple[str, ...]) -> Table:
    """A .tsv result table (the '#' lines, one header line, then rows): every column as a float array.

    Source: scattering_summary.tsv, s_matrix.tsv, probabilities_*.tsv (docs/output_files.md, "Results").
    Raises ReadError when a required column is missing or a row does not parse.
    """
    lines = [(n, line) for n, line in enumerate(text_of(path).splitlines(), start=1) if line and line[0] != "#"]
    if not lines:
        raise ReadError(f"{path}: no header line; the file is empty or not a result table")
    header = lines[0][1].split("\t")
    missing = [name for name in required if name not in header]
    if missing:
        raise ReadError(f"{path}: no column {', '.join(missing)}; the file is from another program version")
    rows = [parse_row(line.split("\t"), len(header), path, n) for n, line in lines[1:]]
    values = np.array(rows, dtype=float).reshape(len(rows), len(header))
    return Table({name: values[:, i] for i, name in enumerate(header)}, first_non_finite(values))


SUMMARY_COLUMNS = ("E_coll", "E_total", "E_z", "P_scattered", "P_reaction", "E_vib_mean", "E_rot_mean", "E_dif_mean")
S_MATRIX_COLUMNS = ("E_coll", "E_total", "v", "j", "mJ", "n", "m", "order", "Re_S", "Im_S", "P", "phase")


def read_scattering_summary(path: Path) -> Table:
    """scattering_summary.tsv: one row per analysis energy (energies in eV; P_v<v> columns per level)."""
    return read_table(path, SUMMARY_COLUMNS)


def read_s_matrix(path: Path) -> Table:
    """s_matrix.tsv: one row per energy and open channel (energies in eV, S dimensionless)."""
    return read_table(path, S_MATRIX_COLUMNS)


@dataclass
class GasPhaseStates:
    """gas_phase_states.dat: every rovibrational state of the gas phase on the r grid."""

    r_bohr: np.ndarray  # (n_r,)
    v: np.ndarray  # (n_states,) int
    j: np.ndarray  # (n_states,) int
    energy_hartree: np.ndarray  # (n_states,)
    wave_functions: np.ndarray  # (n_states, n_r), in bohr^-1/2


R_GRID_LINE = "# r grid:"


def read_gas_phase_states(path: Path) -> GasPhaseStates:
    """gas_phase_states.dat: the r grid (its '# r grid:' line), then per state v, j, E and psi(r).

    Source: gas_phase_states.dat (docs/output_files.md, "Data files"). Units: bohr, hartree.
    """
    text = text_of(path)
    grid = [line for line in text.splitlines() if line.startswith(R_GRID_LINE)]
    if len(grid) != 1:
        raise ReadError(f"{path}: {len(grid)} '{R_GRID_LINE}' lines, not 1; the file is from another program version")
    r = np.array([float(v) for v in grid[0].removeprefix(R_GRID_LINE).split()])
    rows = numbers(path, text, 3 + len(r))
    return GasPhaseStates(r, rows[:, 0].astype(int), rows[:, 1].astype(int), rows[:, 2], rows[:, 3:])


@dataclass
class AsymptoticPotential:
    """asymptotic_potential.dat: the gas-phase potential V(r) at Z = zasymp, theta = pi/2."""

    r_bohr: np.ndarray
    v_hartree: np.ndarray
    v_ev: np.ndarray


def read_asymptotic_potential(path: Path) -> AsymptoticPotential:
    """asymptotic_potential.dat: columns ir, r (bohr), V (hartree), V (eV)."""
    rows = numbers(path, text_of(path), 4)
    return AsymptoticPotential(rows[:, 1], rows[:, 2], rows[:, 3])
