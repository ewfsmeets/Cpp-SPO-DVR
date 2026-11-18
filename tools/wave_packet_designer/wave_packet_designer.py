#!/usr/bin/env python3
"""Wave-packet designer, command line: checks the Z and r grids, the initial wave packet and the absorbers.

Usage (from any directory):
  python3 -I -B tools/wave_packet_designer/wave_packet_designer.py check <input.json>
      [--zasymp Z --point-group G]     (both required for CRP data sets: SystemID CRPlibrary, CRPDistrtSurface)

Exit status: 0 spodvr accepts the input, 1 spodvr would stop on it, 2 the file or the options are unusable.
Needs: numpy. Docs: tools/wave_packet_designer/README.md
"""

import argparse
import sys
from pathlib import Path

sys.dont_write_bytecode = True  # no __pycache__ in the repository, also when not run with -B
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's modules, found next to this file

from evaluation import SETTINGS_STAGE, evaluate  # noqa: E402
from input_file import InputFileError, load_input_file  # noqa: E402
from program.json_input import parse_input  # noqa: E402
from report import format_check_report, format_reading_failure  # noqa: E402
from settings import Settings  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    args = parse_arguments(argv)
    return run_check(Path(args.input), Settings(z_asymptotic=args.zasymp, point_group=args.point_group))


def parse_arguments(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="what spodvr does with the input, and the packet, grids and absorbers")
    check.add_argument("input", help="input.json")
    check.add_argument("--zasymp", type=float, help="Z_asymptotic (bohr) of a CRP data set")
    check.add_argument("--point-group", help="point group of a CRP data set (C2v, C4v, C6v, C3v, C3, none)")
    return parser.parse_args(argv)


def run_check(path: Path, settings: Settings) -> int:
    """Prints the report of one input file; returns the exit status."""
    try:
        input_file = load_input_file(path)
    except InputFileError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    parsed = parse_input(input_file.values)
    if parsed.parameters is None:
        print(format_reading_failure(str(path), parsed.findings), end="")
        return 1
    evaluation = evaluate(parsed.parameters, settings)
    if evaluation.stopped_at == SETTINGS_STAGE:  # a CRP data set without its Z_asymptotic and point group
        print(
            f"error: SystemID {parsed.parameters.system_id} takes Z_asymptotic and the point group from its CRP "
            "data set: give both with --zasymp and --point-group (the program's set-up output prints them)",
            file=sys.stderr,
        )
        return 2
    print(format_check_report(evaluation, str(path), parsed.findings), end="")
    return 0 if evaluation.accepted else 1


if __name__ == "__main__":
    sys.exit(main())
