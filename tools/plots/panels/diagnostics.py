"""The verdict as a text panel."""

import textwrap

from matplotlib.axes import Axes

from compute.diagnostics import Finding


def verdict_line(f: Finding) -> str:
    """ "WARNING  name: value (limit ...). meaning [inputs/README.md section]", leaving out empty parts."""
    line = f"{f.severity.upper():12} {f.name}: {f.value}"
    line += f" (limit {f.limit})." if f.limit else "."
    line += f" {f.meaning}" if f.meaning else ""
    line += f" [inputs/README.md {f.section}]" if f.section else ""
    return line


def verdict_lines(findings: list[Finding]) -> list[str]:
    """One line per finding, worst first."""
    return [verdict_line(f) for f in findings] if findings else ["no findings"]


VERDICT_CHARACTERS = 118  # characters per line of the verdict panel (monospace 7.5 pt on a report-wide figure)
VERDICT_LINE_IN = 0.13  # height of one of its lines (inches)


def wrapped_verdict(findings: list[Finding]) -> list[str]:
    """The verdict's lines wrapped to VERDICT_CHARACTERS, continuation lines indented under the name."""
    lines = []
    for line in verdict_lines(findings):
        lines += textwrap.wrap(line, VERDICT_CHARACTERS, subsequent_indent=" " * 13)
    return lines


def draw_verdict(ax: Axes, lines: list[str]) -> None:
    """The verdict (wrapped_verdict) as text on an axes without axis lines."""
    ax.axis("off")
    ax.set_title("Verdict: is the run complete and converged?", loc="left")
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=7.5, family="monospace", transform=ax.transAxes)
