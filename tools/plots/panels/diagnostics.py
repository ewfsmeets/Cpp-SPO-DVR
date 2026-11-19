"""The verdict as a text panel."""

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


def draw_verdict(ax: Axes, findings: list[Finding]) -> None:
    """The verdict as text on an axes without axis lines."""
    ax.axis("off")
    ax.set_title("Verdict: is the run complete and converged?", loc="left")
    text = "\n".join(verdict_lines(findings))
    ax.text(0.0, 1.0, text, va="top", ha="left", fontsize=7.5, family="monospace", wrap=True, transform=ax.transAxes)
