"""fine_analysis/grid.json: the spectral grid of the fine analysis (its tables are read by readers/tables.py).

Shadows: src/report/fine_analysis_report.cpp (write_fine_analysis_files)
"""

import json
from dataclasses import dataclass
from pathlib import Path

from readers.errors import ReadError


@dataclass
class FineGrid:
    """The fine analysis's grid E_z = k * resolution, k = first_index .. last_index (the points inside its
    95 % confidence interval), and the analysis energies it could not interpolate."""

    spectral_resolution_ev: float
    first_index: int
    last_index: int
    effective_time_au: float
    records_used: int
    skipped: list[tuple[float, str]]  # (E_z in eV, reason)


def read_fine_grid(path: Path) -> FineGrid:
    """grid.json of a fine analysis.

    Source: fine_analysis/grid.json (docs/output_files.md, "Fine analysis").
    """
    if not path.is_file():
        raise ReadError(f"{path}: no such file; give a run with alsoAnalysis or onlyAnalysis")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        skipped = [(float(s["E_z"]), str(s["reason"])) for s in data["skipped_energies"]]
        return FineGrid(
            float(data["spectral_resolution_eV"]),
            int(data["first_index"]),
            int(data["last_index"]),
            float(data["effective_time"]),
            int(data["records_used"]),
            skipped,
        )
    except json.JSONDecodeError as error:
        raise ReadError(f"{path}: not JSON ({error}); the file is damaged") from None
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{path}: no or a wrong field {error}; the file is damaged or from another version") from None
