#!/usr/bin/env python3
"""Tests of spodvr_compare_wave (tools/compare/compare_wave.cpp; ctest CompareWave.*).

Each case writes small wave-packet files (raw doubles, re/im pairs, as spodvr_dump_wave writes
them) into WORK, runs the tool on them and checks its exit status and output line:
  finite      identical files pass ("bitwise identical"); a difference of 1e-14 relative fails
              with the default tolerance 0 and passes with 1e-13; files of different sizes fail
  non-finite  a NaN in b, an inf in a, a NaN at the same element of both (bitwise equal files)
              and a NaN with tolerance 1e300: each fails, naming the file and the element
Usage: python3 -I -B tools/compare/test_compare_wave.py TOOL CASE WORK
  TOOL  the spodvr_compare_wave executable;  CASE  finite or non-finite;  WORK  scratch directory
  (recreated)
Exit status: 0 the case behaved as expected, 1 it did not, 2 bad arguments.
Environment: none read. Run time: well under a second.
Docs: docs/performance.md, section "Benchmarking checklist" (the final wave packet of two builds).
"""

import argparse
import math
import shutil
import struct
import subprocess
import sys
from pathlib import Path

CASES = ("finite", "non-finite")
WAVE = [1.0, 0.0, 0.5, -0.5, 0.25, 0.75]  # three complex values


def write(path, values):
    """Writes the values as native doubles."""
    path.write_bytes(struct.pack(f"={len(values)}d", *values))
    return path


def run(tool, a, b, tolerance=None):
    """Runs the tool; returns (exit status, stdout)."""
    command = [str(tool), str(a), str(b)] + ([] if tolerance is None else [tolerance])
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    return result.returncode, result.stdout.strip()


def expect(what, status, output, want_status, want_output):
    """Raises if the exit status differs or want_output is not in the output."""
    if status != want_status or want_output not in output:
        raise AssertionError(f"{what}: exit {status}, output {output!r}; expected exit {want_status} "
                             f"and {want_output!r}")


def finite_case(tool, work):
    a = write(work / "a.bin", WAVE)
    b = write(work / "b.bin", WAVE)
    expect("identical", *run(tool, a, b), 0, "bitwise identical: max |a - b| / max |a| = 0.000e+00")
    changed = list(WAVE)
    changed[3] = changed[3] * (1.0 + 1e-14)
    c = write(work / "c.bin", changed)
    expect("1e-14 apart, tolerance 0", *run(tool, a, c), 1, "different: max |a - b| / max |a|")
    expect("1e-14 apart, tolerance 1e-13", *run(tool, a, c, "1e-13"), 0, "different: max |a - b| / max |a|")
    d = write(work / "d.bin", WAVE[:4])
    expect("sizes", *run(tool, a, d), 1, "DIFFERENT: sizes 6 and 4")


def non_finite_case(tool, work):
    a = write(work / "a.bin", WAVE)
    with_nan = list(WAVE)
    with_nan[3] = math.nan
    b = write(work / "b.bin", with_nan)
    expect("NaN in b", *run(tool, a, b), 1,
           f"DIFFERENT: non-finite value nan in {b} at element 3 (complex value 1, imaginary part)")
    expect("NaN in b, tolerance 1e300", *run(tool, a, b, "1e300"), 1, "DIFFERENT: non-finite value nan")
    with_inf = list(WAVE)
    with_inf[0] = math.inf
    c = write(work / "c.bin", with_inf)
    expect("inf in a", *run(tool, c, a), 1,
           f"DIFFERENT: non-finite value inf in {c} at element 0 (complex value 0, real part)")
    expect("NaN at the same element of both", *run(tool, b, b), 1,
           f"DIFFERENT: non-finite value nan in {b} at element 3")


def main(argv):
    """Runs the named case in WORK."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tool", metavar="TOOL", type=Path, help="the spodvr_compare_wave executable")
    parser.add_argument("case", metavar="CASE", choices=CASES, help="the case to run")
    parser.add_argument("work", metavar="WORK", type=Path, help="scratch directory (recreated)")
    args = parser.parse_args(argv)
    if not args.tool.is_file():
        print(f"error: {args.tool} is not a file; build the target spodvr_compare_wave", file=sys.stderr)
        return 1
    work = args.work.resolve()
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    try:
        (finite_case if args.case == "finite" else non_finite_case)(args.tool, work)
    except AssertionError as error:
        print(f"FAIL: {args.case}: {error}")
        return 1
    print(f"PASS: {args.case}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
