"""run_info.json, the metadata of a run.

Shadows: src/report/run_info_file.cpp (write_run_info)
"""

import json
from dataclasses import dataclass
from pathlib import Path

from readers.errors import ReadError


@dataclass
class RunInfo:
    """run_info.json of one run directory (or segment). Fields by the file's names, e.g.
    data["wave_packet"]["kz"]; atomic units unless the name ends in a unit (docs/output_files.md,
    "Run metadata"). The fields below are the ones nearly every figure needs."""

    directory: Path
    data: dict
    hartree_to_ev: float  # the program's factor (the legacy 27.2114)
    t_start_au: float
    t_end_au: float


FORMAT_VERSION = 1


def read_run_info(directory: Path) -> RunInfo:
    """run_info.json of a run directory (an analysis-only run: its fine_analysis/run_info.json).

    Source: run_info.json (docs/output_files.md, "Run metadata").
    """
    path = directory / "run_info.json"
    if not path.is_file():
        path = directory / "fine_analysis" / "run_info.json"
    if not path.is_file():
        raise ReadError(f"{directory}: no run_info.json; give a directory of a run made since run_info.json exists")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("format_version") != FORMAT_VERSION:
            raise ReadError(
                f"{path}: format_version {data.get('format_version')}, not {FORMAT_VERSION}; update the tool"
            )
        times = data["times"]
        return RunInfo(
            directory,
            data,
            float(data["unit_factors"]["hartree_to_ev"]),
            float(times["t_start"]),
            float(times["t_end"]),
        )
    except json.JSONDecodeError as error:
        raise ReadError(f"{path}: not JSON ({error}); the file is damaged") from None
    except KeyError as error:
        raise ReadError(f"{path}: no field {error}; the file is damaged or from another program version") from None
    except (TypeError, ValueError) as error:
        raise ReadError(f"{path}: a field has the wrong type ({error}); the file is damaged") from None


@dataclass
class Lattice:
    """The surface lattice of a run (run_info.json lattice): the lengths of the two lattice vectors and the skew
    angle between them."""

    a_bohr: float
    b_bohr: float
    skew_angle_rad: float


def lattice_of(info: RunInfo) -> Lattice:
    """The run's surface lattice."""
    lattice = info.data["lattice"]
    return Lattice(float(lattice["a"]), float(lattice["b"]), float(lattice["skew_angle_rad"]))
