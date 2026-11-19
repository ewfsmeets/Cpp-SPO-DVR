"""Plots of one spodvr run: the command line.

Usage (from inside the repository, where uv finds the environment):
  uv run python -I -B tools/plots/spodvr_plot.py list
  uv run python -I -B tools/plots/spodvr_plot.py plot <name> <run directory>... --out DIR
  uv run python -I -B tools/plots/spodvr_plot.py report <run directory>... --out DIR
  uv run python -I -B tools/plots/spodvr_plot.py verdict <run directory>...
A run split into restart segments: give every segment's directory (any order).
Exit status: 0 done (the verdict is advice and never changes it); 1 a figure or the report's verdict could
not be made (the report says which); 2 bad arguments, an unreadable run, or an --out directory that holds a
file of the same name that the tool did not write (nothing is written then).
Environment: none read. Run time: seconds per figure. Needs: numpy, matplotlib (the repository's environment).
Docs: tools/plots/README.md
"""

import argparse
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's packages import from its directory

from figures import PLOTS, default_options, has_files  # noqa: E402
from panels.diagnostics import verdict_lines  # noqa: E402
from readers.errors import ReadError  # noqa: E402
from readers.run_info import RunInfo  # noqa: E402
from readers.segments import joined_info, ordered_segments  # noqa: E402
from report import write_figures, write_report  # noqa: E402
from verdict import run_verdict  # noqa: E402


def arguments(argv: list[str]) -> argparse.Namespace:
    """The command line."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="the figures and the files each needs")
    plot = commands.add_parser("plot", help="one figure")
    plot.add_argument("name", choices=[name for name, _, _ in PLOTS])
    report = commands.add_parser("report", help="report.md with the verdict and every figure")
    verdict = commands.add_parser("verdict", help="the verdict as text")
    for command in (plot, report, verdict):
        command.add_argument("runs", nargs="+", type=Path, help="the run directory, or every segment's")
    for command in (plot, report):
        command.add_argument("--out", type=Path, required=True, help="the directory to write into")
    return parser.parse_args(argv)


def list_figures() -> int:
    """Prints one line per figure: name, what it shows, the files it needs."""
    for name, figure, needed in PLOTS:
        files = ", ".join(needed) if needed else "what the run has"
        print(f"{name:12} {figure.__doc__.splitlines()[0]} (needs {files})")
    return 0


def plot(name: str, runs: list[Path], info: RunInfo, out: Path) -> int:
    """Draws one figure into `out` (replacing only files the tool wrote there) and prints the files written."""
    figure, needed = next((f, n) for p, f, n in PLOTS if p == name)
    if not has_files(runs[-1], needed):
        raise ReadError(f"{runs[-1]}: no {' or '.join(needed)}; the figure {name} needs it")
    for path in write_figures(figure(runs, info, default_options()), out):
        print(path)
    return 0


def report(runs: list[Path], info: RunInfo, out: Path) -> int:
    """Writes the report; prints its path and the problems."""
    problems = write_report(runs, info, out)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    print(out / "report.md")
    return 1 if problems else 0


def main(argv: list[str]) -> int:
    """Runs one command; returns the exit status."""
    args = arguments(argv)
    if args.command == "list":
        return list_figures()
    try:
        infos, warnings = ordered_segments(args.runs)
        for warning in warnings:
            print(f"warning: {warning}", file=sys.stderr)
        runs, info = [i.directory for i in infos], joined_info(infos)
        if args.command == "verdict":
            print("\n".join(verdict_lines(run_verdict(runs, info))))
            return 0
        if args.command == "plot":
            return plot(args.name, runs, info, args.out)
        return report(runs, info, args.out)
    except (ReadError, OSError) as error:  # OSError: e.g. --out names a file, or a directory not writable
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
