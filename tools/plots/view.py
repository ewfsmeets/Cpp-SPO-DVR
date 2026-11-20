#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib>=3.9", "PyQt6"]
# ///
"""The interactive viewer of one spodvr run's frames/ (dtFrames): the probability over the unit cell in a window of
(Z, r), and rho(Z, r) above a grid point picked by a click, against time (viewer.py).

Usage (uv installs numpy, matplotlib and PyQt6 once, from the header above):
  uv run tools/plots/view.py <run directory>... [--out DIR]
  uv run tools/plots/view.py <run directory>... --self-test DIR
A run split into restart segments: give every segment's directory (any order). Without a window (no display, or
MPLBACKEND set to a backend without one, e.g. Agg) it says so and writes the frame movies of the sites
(spodvr_plot.py's frame_movie) into --out instead (default: the working directory). --self-test opens the window,
drives it (as opened, the middle frame, a click on the second site, one Play step, V_min's contours at the
surface), saves a PNG after each of the five into DIR and closes it.
Exit status: 0 done; 1 the self-test did not write its five PNGs; 2 bad arguments, a run without frames/ or not
readable, an --out holding a file of the same name the tool did not write, or --self-test without a window.
Environment: DISPLAY or WAYLAND_DISPLAY (the window; WSLg sets both); MPLBACKEND (Agg: no window, the movies).
Needs: the PEP 723 header's packages. Docs: tools/plots/README.md ("The viewer")
"""

import argparse
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's packages import from its directory

from figures import default_options, figure_frame_movie  # noqa: E402
from readers.errors import ReadError  # noqa: E402
from readers.run_info import RunInfo  # noqa: E402
from readers.segments import joined_info, ordered_segments  # noqa: E402
from report import write_figures  # noqa: E402
from viewer import SELF_TEST_PNGS, has_window_backend, open_window, run_self_test, show  # noqa: E402


def arguments(argv: list[str]) -> argparse.Namespace:
    """The command line."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+", type=Path, help="the run's directory (every segment's)")
    parser.add_argument("--out", type=Path, default=Path("."), help="where the movies go without a window")
    parser.add_argument("--self-test", type=Path, metavar="DIR", help="drive the window, save PNGs into DIR, close")
    return parser.parse_args(argv)


def write_movies(runs: list[Path], info: RunInfo, out: Path) -> int:
    """Without a window: the frame movies of the sites into `out`, their paths printed."""
    print("no window here (no display or no GUI backend): writing the frame movies instead", file=sys.stderr)
    for path in write_figures(figure_frame_movie(runs, info, default_options()), out):
        print(path)
    return 0


def self_test(runs: list[Path], info: RunInfo, out: Path) -> int:
    """The window driven and closed; 1 unless all its PNGs were written."""
    saved = run_self_test(runs, info, out)
    for path in saved:
        print(path)
    if len(saved) != len(SELF_TEST_PNGS):
        print(f"error: the self-test wrote {len(saved)} of {len(SELF_TEST_PNGS)} PNGs", file=sys.stderr)
        return 1
    return 0


def main(argv: list[str]) -> int:
    """Opens the viewer (or writes the movies, or runs the self-test); returns the exit status."""
    args = arguments(argv)
    try:
        infos, warnings = ordered_segments(args.runs)
        for warning in warnings:
            print(f"warning: {warning}", file=sys.stderr)
        runs, info = [i.directory for i in infos], joined_info(infos)
        if not has_window_backend():
            if args.self_test is not None:
                print("error: no window here (no display or no GUI backend): --self-test needs one", file=sys.stderr)
                return 2
            return write_movies(runs, info, args.out)
        if args.self_test is not None:
            return self_test(runs, info, args.self_test)
        _figure_and_viewer = open_window(runs, info)  # kept: widgets without a reference stop responding
        show()
        return 0
    except (ReadError, OSError) as error:  # OSError: e.g. --out names a file, or a directory not writable
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
