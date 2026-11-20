"""Panels of the norm budget (norm_budget/): where the norm went against time, and whether the absorbers' losses
account for the norm that left the grids."""

from matplotlib.axes import Axes

from compute.norm_budget import Budget, absorber_label
from panels.time_series import mark_segments
from readers.norm_budget import BudgetAbsorber
from style import GRID_COLOR, INK, LINE_STYLES, MARKERS, MUTED, apply_scale, axis_label, masked, scale_and_limits


def absorber_color(absorber: BudgetAbsorber) -> str:
    """An absorber's colour by what it acts on: ROpt the r colour, a Z absorber its grid's."""
    return GRID_COLOR["r"] if absorber.name == "ROpt" else GRID_COLOR[absorber.grid]


def draw_where_the_norm_went(ax: Axes, legend_ax: Axes, budget: Budget) -> None:
    """Against time (window ends): the norm on the grids (black), the norm that left them, N(0) - N(t) (muted,
    dashed), and each absorber's loss summed so far (its colour by absorber_color, a line style each); y log when
    they span more than 1e3, values below its lower limit not drawn. The legend goes into legend_ax (its axis turned
    off)."""
    scale = scale_and_limits([budget.norm, budget.left_grids, budget.lost_by_absorber])
    every = max(1, len(budget.times_au) // 10)
    ax.plot(budget.times_au, masked(budget.norm, scale.lower), color=INK, label="on the grids")
    ax.plot(
        budget.times_au, masked(budget.left_grids, scale.lower), color=MUTED, linestyle="--", label="left the grids"
    )
    for k, absorber in enumerate(budget.absorbers):
        ax.plot(
            budget.times_au,
            masked(budget.lost_by_absorber[:, k], scale.lower),
            color=absorber_color(absorber),
            linestyle=LINE_STYLES[k % len(LINE_STYLES)],
            marker=MARKERS[k % len(MARKERS)],
            markevery=every,
            markersize=3,
            linewidth=1.2,
            label=f"taken by {absorber_label(absorber)}",
        )
    apply_scale(ax, scale)
    mark_segments(ax, budget.starts_au)
    ax.set_xlabel(axis_label("t", "a.u."))
    ax.set_ylabel(axis_label("Probability"))
    ax.set_title("Where did the norm go?")
    legend_ax.axis("off")
    legend_ax.legend(*ax.get_legend_handles_labels(), loc="center left", fontsize=7)


def draw_budget_residual(ax: Axes, budget: Budget, title: str) -> None:
    """Against time (window ends): the residual, the norm that left the grids minus the absorbers' summed losses
    (black), and its error estimate, the windows' loss_error summed (a muted band of plus and minus it); linear
    (signed), as probabilities."""
    error = budget.error_estimate
    ax.fill_between(budget.times_au, -error, error, color=MUTED, alpha=0.3, linewidth=0, label="error estimate")
    ax.plot(budget.times_au, budget.residual, color=INK, label="residual")
    ax.axhline(0.0, color=MUTED, linewidth=0.8)
    apply_scale(ax, scale_and_limits([budget.residual, error, -error], signed=True))
    mark_segments(ax, budget.starts_au)
    ax.legend(loc="best")
    ax.set_xlabel(axis_label("t", "a.u."))
    ax.set_ylabel(axis_label("Probability"))
    ax.set_title(title)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
