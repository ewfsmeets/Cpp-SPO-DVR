"""norm_budget/ (NormBudget): per analysis window (t - dtAna, t], what each absorber removed (each cell's loss for a
constant inflow during the window) with an error indicator, and the norm on both grids at the window's two ends; a
series directory read by readers/arrays.py (docs/output_files.md, "norm_budget/"). This module knows the directory's name,
its datasets (compute uses them by their C++ names) and the absorbers they list.

Shadows: src/report/norm_budget_files.cpp (budget_datasets, budget_manifest_fields)
"""

from dataclasses import dataclass
from pathlib import Path

from readers.arrays import Series
from readers.errors import ReadError
from readers.segments import read_joined_series

DIRECTORY = "norm_budget"
DATASETS = ("loss", "loss_error", "norm_start", "norm_end")


@dataclass(frozen=True, kw_only=True)
class BudgetAbsorber:
    """An entry of the loss datasets: the absorber's input name and the grid it acts on."""

    name: str  # ZOpt, ZOptLow, ROpt, ZSpOpt, ZSpOptLow
    grid: str  # scattering or specular (ROpt acts on both: one entry per grid)


@dataclass(frozen=True, kw_only=True)
class NormBudgetRecords:
    """norm_budget/ of a run (its segments joined): one record per window, at the window's end time."""

    series: Series  # times_au and the datasets by name, record first (see DATASETS)
    absorbers: list[BudgetAbsorber]  # in the order of the loss datasets' last axis
    steps_per_window: int
    starts_au: list[float]  # the segments' start times


def read_norm_budget(directories: list[Path]) -> NormBudgetRecords:
    """norm_budget/ of a run's segments, joined (readers.segments.read_joined_series; the windows of a restart
    continue those of the run before, none repeated); ReadError when a dataset or the absorbers are missing, or the
    loss datasets do not hold one value per absorber."""
    series, starts = read_joined_series(directories, DIRECTORY)
    source = series.directory / "manifest.json"
    missing = [d for d in DATASETS if d not in series.arrays]
    if missing:
        raise ReadError(f"{source}: no {', '.join(missing)}; give norm_budget/ written by spodvr with NormBudget")
    try:
        absorbers = [BudgetAbsorber(name=str(a["name"]), grid=str(a["grid"])) for a in series.manifest["absorbers"]]
    except (KeyError, TypeError) as error:
        raise ReadError(f"{source}: no valid absorbers ({error}); give norm_budget/ written by spodvr") from None
    try:
        steps = int(series.manifest["window"]["steps"])
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{source}: no valid window steps ({error}); give norm_budget/ written by spodvr") from None
    for name in ("loss", "loss_error"):
        if series.arrays[name].shape[1:] != (len(absorbers),):
            shape = series.arrays[name].shape[1:]
            raise ReadError(f"{source}: {name} has records of shape {shape}, not ({len(absorbers)},) (the absorbers)")
    return NormBudgetRecords(series=series, absorbers=absorbers, steps_per_window=steps, starts_au=starts)
