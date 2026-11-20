"""compute/norm_budget.py and the verdict's budget finding on a synthetic norm_budget/: the summed losses, the
residual and its error by hand, each with a slip that a test must catch."""

from dataclasses import replace
from pathlib import Path

import numpy as np

from compute.diagnostics import BUDGET_LIMIT, budget_open
from compute.norm_budget import budget_check, summed_budget
from readers.arrays import Series
from readers.norm_budget import BudgetAbsorber, NormBudgetRecords


def synthetic_records() -> NormBudgetRecords:
    """Three windows, two absorbers. The norm falls 1.0 -> 0.9 -> 0.75 -> 0.7; the losses (0.07, 0.03), (0.11, 0.04),
    (0.03, 0.02) with the error indicators (0.01, 0), (0.01, 0), (0, 0.01)."""
    arrays = {
        "loss": np.array([[0.07, 0.03], [0.11, 0.04], [0.03, 0.02]]),
        "loss_error": np.array([[0.01, 0.0], [0.01, 0.0], [0.0, 0.01]]),
        "norm_start": np.array([1.0, 0.9, 0.75]),
        "norm_end": np.array([0.9, 0.75, 0.7]),
    }
    series = Series(Path("norm_budget"), {}, np.array([10.0, 20.0, 30.0]), arrays, None)
    absorbers = [BudgetAbsorber(name="ZOpt", grid="scattering"), BudgetAbsorber(name="ROpt", grid="specular")]
    return NormBudgetRecords(series=series, absorbers=absorbers, steps_per_window=5, starts_au=[0.0])


def test_budget_by_hand() -> None:
    """The losses summed: ZOpt 0.07, 0.18, 0.21 and ROpt 0.03, 0.07, 0.09; the norm that left, 0.1, 0.25, 0.3, minus
    their sums 0.1, 0.25, 0.30 is 0 at every window; the error estimates (the indicators summed over the absorbers and
    windows) 0.01, 0.02, 0.03. Summing one absorber only (a slip) leaves a residual of 0.03 at the first window."""
    budget = summed_budget(synthetic_records())
    assert np.allclose(budget.lost_by_absorber, [[0.07, 0.03], [0.18, 0.07], [0.21, 0.09]], rtol=0, atol=1e-15)
    assert np.allclose(budget.left_grids, [0.1, 0.25, 0.3], rtol=0, atol=1e-15)
    assert np.allclose(budget.residual, 0.0, rtol=0, atol=1e-15)
    assert np.allclose(budget.error_estimate, [0.01, 0.02, 0.03], rtol=0, atol=1e-15)
    assert [a.name for a in budget.absorbers] == ["ZOpt", "ROpt"] and budget.initial_norm == 1.0
    one_absorber = budget.left_grids - budget.lost_by_absorber[:, 0]
    assert abs(one_absorber[0] - 0.03) < 1e-15


def test_the_check_and_the_finding() -> None:
    """The budget closes when |residual| <= error estimate + BUDGET_LIMIT x the initial norm at every window. With
    the initial norm 0.5, a residual of -(0.02 + 0.4e-4) at the second window (error estimate 0.02) still closes, one
    of -(0.02 + 0.6e-4) does not (a limit not scaled by the initial norm, or |residual| without abs, slips here);
    NaN is the worst window. A closed budget gives no finding, an open one a warning naming its time."""
    budget = replace(summed_budget(synthetic_records()), initial_norm=0.5)
    assert budget_check(budget).closes
    inside = replace(budget, residual=np.array([0.0, -(0.02 + 0.4e-4), 0.0]))
    outside = replace(budget, residual=np.array([0.0, -(0.02 + 0.6e-4), 0.0]))
    assert budget_check(inside).closes
    check = budget_check(outside)
    assert not check.closes and check.window == 1 and abs(check.excess - 0.1e-4) < 1e-15
    assert budget_check(replace(budget, residual=np.array([0.0, np.nan, 0.0]))).window == 1
    assert budget_open(excess=-1e-3, residual=0.0, error_estimate=0.02, t_au=10.0) is None
    finding = budget_open(excess=1e-3, residual=0.002, error_estimate=0.001, t_au=20.0)
    assert finding.severity == "warning" and finding.name == "norm budget does not close"
    assert finding.value == "|residual| 2.00e-03 at t = 20 a.u., error estimate 1.00e-03"
    assert finding.limit == f"|residual| <= error estimate + {BUDGET_LIMIT:g} of the initial norm"
    assert budget_open(excess=float("nan"), residual=float("nan"), error_estimate=0.0, t_au=10.0) is not None
