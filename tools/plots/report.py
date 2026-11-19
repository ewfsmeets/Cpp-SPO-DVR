"""The report of a run (layer 6): report.md with the verdict first, then every figure that the run's files
allow, each with its title and explanation, and a line for each figure not drawn. The tool writes only into
the output directory, and there replaces only files it wrote before (listed in written_by_spodvr_plot.json):
every name is checked before anything is written (also by `plot`, through write_figures)."""

import json
from pathlib import Path

from matplotlib.figure import Figure

from explanations import FIGURES, TITLES
from figures import PLOTS, default_options, has_files
from panels.diagnostics import verdict_lines
from readers.errors import ReadError
from readers.run_info import RunInfo
from style import save_figure
from verdict import run_verdict

MANIFEST = "written_by_spodvr_plot.json"  # the files the tool wrote into an output directory


def own_files(out: Path) -> set[str]:
    """The files the tool wrote into `out` before (its manifest), which it may replace."""
    manifest = out / MANIFEST
    if not manifest.is_file():
        return set()
    try:
        return set(json.loads(manifest.read_text(encoding="utf-8"))["files"])
    except (ValueError, KeyError, TypeError):
        raise ReadError(f"{manifest}: damaged; delete it, or give another --out") from None


def check_writable(out: Path, names: list[str]) -> None:
    """ReadError when a file of these names exists in `out` and the tool did not write it."""
    allowed = own_files(out) | {MANIFEST}
    clash = sorted(n for n in names if (out / n).exists() and n not in allowed)
    if clash:
        raise ReadError(f"{out}: {', '.join(clash)} exist and were not written by this tool; give another --out")


def record_written(out: Path, names: list[str]) -> None:
    """Adds `names` to the manifest of `out` (keeping the earlier entries whose files still exist)."""
    kept = {n for n in own_files(out) if (out / n).exists()}
    (out / MANIFEST).write_text(json.dumps({"files": sorted(kept | set(names))}, indent=1) + "\n", encoding="utf-8")


def claim(out: Path, names: list[str]) -> None:
    """Creates `out`, checks that the tool may write every one of `names` there, and records them in the
    manifest before anything is written (so that a file written before a crash is still the tool's)."""
    out.mkdir(parents=True, exist_ok=True)
    check_writable(out, names)
    record_written(out, names)


def figure_files(figures: list[tuple[str, Figure]]) -> list[str]:
    """The file names save_figure writes for (stem, figure) pairs."""
    return [f"{stem}{suffix}" for stem, _ in figures for suffix in (".png", ".pdf")]


def write_figures(figures: list[tuple[str, Figure]], out: Path) -> list[Path]:
    """Saves (stem, figure) pairs as PNG and PDF into `out` (claimed first); returns the paths."""
    claim(out, figure_files(figures))
    return [path for stem, fig in figures for path in save_figure(fig, out / stem)]


def verdict_section(runs: list[Path], info: RunInfo) -> tuple[list[str], list[str]]:
    """The report's verdict section (Markdown lines), and the problem when the verdict could not be made."""
    try:
        lines, problems = verdict_lines(run_verdict(runs, info)), []
    except ReadError as error:
        lines, problems = [f"not made: {error}"], [f"verdict: {error}"]
    return ["## Verdict", "", "```text", *lines, "```", ""], problems


def make_figures(runs: list[Path], info: RunInfo) -> tuple[list[tuple[str, str, Figure]], list[str], list[str]]:
    """Every figure the run's files allow except the verdict's (the report prints it as text): (name, stem,
    figure) triples, the Markdown lines for those not drawn, and the problems (figures that failed)."""
    made, not_drawn, problems = [], [], []
    for name, figure, needed in PLOTS:
        if name == "verdict":
            continue
        if not has_files(runs[-1], needed):
            not_drawn.append(f"- {TITLES[name]}: not drawn, the run has no {' or '.join(needed)}")
            continue
        try:
            made += [(name, stem, fig) for stem, fig in figure(runs, info, default_options())]
        except ReadError as error:
            not_drawn.append(f"- {TITLES[name]}: not drawn, {error}")
            problems.append(f"{name}: {error}")
    return made, not_drawn, problems


def write_report(runs: list[Path], info: RunInfo, out: Path) -> list[str]:
    """report.md, the figures and the manifest in `out`; returns the problems (exit status 1). Raises ReadError
    before writing anything when `out` holds a file of the same name that the tool did not write."""
    verdict, problems = verdict_section(runs, info)
    made, not_drawn, figure_problems = make_figures(runs, info)
    figures = [(stem, fig) for _, stem, fig in made]
    claim(out, ["report.md", *figure_files(figures)])
    for stem, fig in figures:
        save_figure(fig, out / stem)
    lines = [f"# {info.data['title']}", "", *verdict]
    for name, stem, _ in made:
        lines += [f"## {TITLES[name]}", "", f"![{TITLES[name]}]({stem}.png)", "", FIGURES[name], ""]
    lines += ["## Not drawn", "", *not_drawn, ""] if not_drawn else []
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    return problems + figure_problems
