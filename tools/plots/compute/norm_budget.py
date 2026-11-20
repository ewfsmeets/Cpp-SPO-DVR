"""The norm budget of a run (norm_budget/, NormBudget): what each absorber took, summed over the analysis windows,
against the norm that left the grids, and how well the two agree.

A window's loss per absorber is the program's: for each cell, exact for a constant inflow during the window (docs/
output_files.md, "norm_budget/"). Its error indicator is gamma^2 / 6 of each cell's inflow, summed: the leading
distance between the loss and the mean of the two one-sided estimates (from the cells at the window's start and at
its end), so it falls with dtAna squared. The residual is the norm that left the grids minus the summed
losses. Beyond the indicator it holds the time stepping's own error (the absorbers act once per time step, between
the kinetic half steps: second order in dtStep, 4.1e-7 of the norm on the hard-wall fixture by t = 6000 a.u.) and
probability that left another way: the angular transform at the surface (docs/output_files.md). The verdict allows
BUDGET_LIMIT.
"""

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from compute.diagnostics import BUDGET_LIMIT
from compute.series import usable_records
from readers.norm_budget import BudgetAbsorber, NormBudgetRecords, read_norm_budget


@dataclass(frozen=True, kw_only=True)
class Budget:
    """The budget at each window's end (windows,): every absorber's loss so far, the norm on the grids, the residual
    and its error estimate, as probabilities."""

    times_au: np.ndarray  # (windows,): the window ends
    absorbers: list[BudgetAbsorber]  # in the order of lost_by_absorber's last axis
    lost_by_absorber: np.ndarray  # (windows, absorbers): each absorber's loss summed up to the window's end
    norm: np.ndarray  # (windows,): the probability on both grids at the window's end
    initial_norm: float  # at the first window's start
    left_grids: np.ndarray  # (windows,): initial_norm - norm, the norm that left the grids
    residual: np.ndarray  # (windows,): left_grids - lost_by_absorber summed over the absorbers
    error_estimate: np.ndarray  # (windows,): the windows' loss_error over the absorbers, summed up to there
    starts_au: list[float]  # the segments' start times


@dataclass(frozen=True, kw_only=True)
class BudgetCheck:
    """Whether the budget closes: |residual| within its error estimate plus BUDGET_LIMIT of the initial norm at every
    window end; the window where it is furthest from that and by how much (a probability, positive where it does not
    close; NaN counts as the furthest)."""

    closes: bool
    window: int
    excess: float


def usable_budget(runs: list[Path]) -> tuple[NormBudgetRecords, tuple[str, ...]]:
    """norm_budget/ of a run's segments with its records before any NaN or inf, and a note when some were left out
    (compute.series.usable_records; ReadError when none is left)."""
    records = read_norm_budget(runs)
    series, notes = usable_records(records.series)
    return replace(records, series=series), notes


def absorber_label(absorber: BudgetAbsorber) -> str:
    """An absorber as the panels name it: "ZOpt (scattering grid)"."""
    return f"{absorber.name} ({absorber.grid} grid)"


def summed_budget(records: NormBudgetRecords) -> Budget:
    """The budget summed over the windows (see Budget)."""
    arrays = records.series.arrays
    lost = np.cumsum(np.asarray(arrays["loss"]), axis=0)
    initial_norm = float(arrays["norm_start"][0])
    norm = np.asarray(arrays["norm_end"])
    left_grids = initial_norm - norm
    return Budget(
        times_au=records.series.times_au,
        absorbers=records.absorbers,
        lost_by_absorber=lost,
        norm=norm,
        initial_norm=initial_norm,
        left_grids=left_grids,
        residual=left_grids - lost.sum(axis=1),
        error_estimate=np.cumsum(np.asarray(arrays["loss_error"]).sum(axis=1)),
        starts_au=records.starts_au,
    )


def excess_over_error(budget: Budget) -> np.ndarray:
    """(windows,): how far |residual| exceeds its error estimate plus BUDGET_LIMIT of the initial norm, as a
    probability; positive where the budget does not close (NaN stays NaN)."""
    return np.abs(budget.residual) - budget.error_estimate - BUDGET_LIMIT * budget.initial_norm


def budget_check(budget: Budget) -> BudgetCheck:
    """Whether the budget closes, and its worst window (the one test the verdict and the figure's title share)."""
    excess = excess_over_error(budget)
    window = int(np.argmax(np.nan_to_num(excess, nan=np.inf)))
    return BudgetCheck(closes=bool(excess[window] <= 0.0), window=window, excess=float(excess[window]))
