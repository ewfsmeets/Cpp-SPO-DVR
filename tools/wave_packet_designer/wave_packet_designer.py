#!/usr/bin/env python3
"""Wave-packet designer, command line: checks the Z and r grids, the initial wave packet and the absorbers,
and estimates the memory of a run.

Usage (from any directory):
  python3 -I -B tools/wave_packet_designer/wave_packet_designer.py check <input.json>
  python3 -I -B tools/wave_packet_designer/wave_packet_designer.py memory <input.json>
      [--threads N] [--representatives N [--operations N]]     (counts from the SURFACE SYMMETRY block)
Both take [--zasymp Z --point-group G], required for CRP data sets (SystemID CRPlibrary, CRPDistrtSurface).

Exit status: 0 spodvr accepts the input, 1 spodvr would stop on it, 2 the file or the options are unusable.
Needs: numpy. Docs: tools/wave_packet_designer/README.md
"""

import argparse
import sys
from collections.abc import Callable
from pathlib import Path

sys.dont_write_bytecode = True  # no __pycache__ in the repository, also when not run with -B
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's modules, found next to this file

from evaluation import SETTINGS_STAGE, Evaluation, evaluate  # noqa: E402
from input_file import InputFileError, load_input_file  # noqa: E402
from program.findings import Finding  # noqa: E402
from program.json_input import parse_input  # noqa: E402
from program.point_group import GROUP_ORDER  # noqa: E402
from report import format_check_report, format_memory_report, format_reading_failure  # noqa: E402
from settings import Settings  # noqa: E402

# The report of each command: (evaluation, input path, findings of the reading) -> text.
REPORTS: dict[str, Callable[[Evaluation, str, list[Finding]], str]] = {
    "check": format_check_report,
    "memory": format_memory_report,
}


def main(argv: list[str] | None = None) -> int:
    args = parse_arguments(argv)
    settings = Settings(
        z_asymptotic=args.zasymp,
        point_group=args.point_group,
        threads=getattr(args, "threads", Settings.threads),
        representatives=getattr(args, "representatives", None),
        operations=getattr(args, "operations", None),
    )
    return run_command(args.command, Path(args.input), settings)


def parse_arguments(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # The input and the values of a CRP data set, common to both commands.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("input", help="input.json")
    common.add_argument("--zasymp", type=float, help="Z_asymptotic (bohr) of a CRP data set")
    common.add_argument("--point-group", choices=sorted(GROUP_ORDER), help="point group of a CRP data set")

    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "check", parents=[common], help="what spodvr does with the input, and the packet, grids and absorbers"
    )
    memory = commands.add_parser("memory", parents=[common], help="the CPU and GPU memory of the run")
    memory.add_argument("--threads", type=positive, default=Settings.threads, help="OpenMP threads (default 8)")
    memory.add_argument(
        "--representatives",
        type=positive,
        help='with UseSymmetry: "(x, y) points computed" of the SURFACE SYMMETRY block',
    )
    memory.add_argument(
        "--operations",
        type=positive,
        help='the number of names under "Operations used", identity included (default: the point group\'s order)',
    )
    args = parser.parse_args(argv)
    if getattr(args, "operations", None) is not None and args.representatives is None:
        parser.error("--operations needs --representatives")
    return args


def positive(text: str) -> int:
    """An integer >= 1 (argparse type)."""
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"{text} is not a positive integer")
    return value


def run_command(command: str, path: Path, settings: Settings) -> int:
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
    if evaluation.stopped_at == SETTINGS_STAGE:  # a setting is missing or cannot hold for this input
        for finding in evaluation.findings:
            print(f"error: {finding.message}", file=sys.stderr)
        print(
            "(settings: --zasymp and --point-group; with memory also --representatives, --operations)", file=sys.stderr
        )
        return 2
    print(REPORTS[command](evaluation, str(path), parsed.findings), end="")
    return 0 if evaluation.accepted else 1


if __name__ == "__main__":
    sys.exit(main())
