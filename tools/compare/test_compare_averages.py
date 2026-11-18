#!/usr/bin/env python3
"""Tests of compare_averages.py (tools/compare/; ctest ComparePropagation.*), the comparison of
tools/compare_propagation.sh.

Each case writes small averages.dat files into WORK, runs the tool on them and checks its exit
status and verdict line:
  equal       equal files pass with tolerance 0; a relative change of 1e-14 passes with 1e-13 and
              fails with 1e-15, naming the column
  non-finite  a NaN in NEW, an inf in BASE, an inf at the same place in both (equal files) and a
              NaN with tolerance 1e300: each fails, naming the file, row and column
  rows        two files without rows and one without rows are refused ("no rows"), a row of four
              numbers is refused, different row counts fail
  arguments   a tolerance that is not a finite number >= 0 (abc, nan, inf, -1) and a missing
              argument: exit 2 with a usage line
Usage: python3 -I -B tools/compare/test_compare_averages.py CASE WORK
  CASE  one of the cases above;  WORK  scratch directory (recreated)
Exit status: 0 the case behaved as expected, 1 it did not, 2 bad arguments.
Environment: none read. Run time: well under a second.
Docs: README.md, section "Testing".
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

TOOL = Path(__file__).resolve().parent / "compare_averages.py"
CASES = ("equal", "non-finite", "rows", "arguments")
HEADER = "# Wave packet averages (Fortran: averages.out)\n# t(a.u.)  norm  <Z>(bohr)  <r>(bohr)  rest_norm\n"
ROWS = [[0.0, 1.0, 11.44, 1.47, 1.0], [100.0, 0.9999996, 11.17, 1.471, 0.9999995], [200.0, 0.99, 10.9, 1.475, 0.98]]


def write(work, name, rows):
    """An averages.dat with the given rows (numbers or strings such as "nan") in %.16e."""
    lines = [" ".join(f"{float(v):.16e}" for v in row) for row in rows]
    path = work / f"{name}.dat"
    path.write_text(HEADER + "".join(line + "\n" for line in lines), encoding="utf-8")
    return path


def changed(row_index, column, value):
    """ROWS with one value replaced."""
    rows = [list(row) for row in ROWS]
    rows[row_index][column] = value
    return rows


def run(*arguments):
    """(exit status, stdout + stderr) of the tool."""
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(TOOL), *map(str, arguments)], capture_output=True, text=True, check=False
    )
    return result.returncode, result.stdout + result.stderr


def expect(outcome, status, pattern):
    """The tool exited with status and printed a line matching pattern."""
    actual_status, output = outcome
    if actual_status != status or not re.search(pattern, output, re.MULTILINE):
        raise SystemExit(
            f"FAIL: expected exit {status} and a line matching {pattern!r}, got exit {actual_status}:\n{output}"
        )
    print(f"as expected (exit {status}): {pattern}")


def equal_case(work):
    base = write(work, "base", ROWS)
    expect(run(base, write(work, "same", ROWS), 0), 0, r"^OK: within 0 of the base$")
    near = write(work, "near", changed(1, 2, 11.17 * (1 + 1e-14)))
    expect(run(base, near, 1e-13), 0, r"^OK: within 1e-13 of the base$")
    expect(run(base, near, 1e-15), 1, r"^DIFFERENT: [0-9.]+e-1[45] > tolerance 1e-15$")
    expect(run(base, near, 1e-15), 1, r"^  <Z> +max relative difference [0-9.]+e-1[45]$")


def non_finite_case(work):
    base = write(work, "base", ROWS)
    expect(
        run(base, write(work, "nan", changed(2, 3, "nan")), 1e-13),
        1,
        r"^DIFFERENT: .*nan\.dat: row 3, <r>: nan is not finite$",
    )
    expect(
        run(write(work, "inf", changed(0, 1, "inf")), base, 1e-13),
        1,
        r"^DIFFERENT: .*inf\.dat: row 1, norm: inf is not finite$",
    )
    both = write(work, "both_inf", changed(1, 4, "-inf"))
    expect(run(both, both, 1e-13), 1, r"^DIFFERENT: .*both_inf\.dat: row 2, rest norm: -inf is not finite$")
    expect(run(base, work / "nan.dat", 1e300), 1, r"^DIFFERENT: .*nan\.dat: row 3, <r>: nan is not finite$")


def rows_case(work):
    empty = write(work, "empty", [])
    expect(run(empty, write(work, "empty2", []), 1e-13), 1, r"^error: .*empty\.dat has no rows; check")
    base = write(work, "base", ROWS)
    expect(run(base, empty, 1e-13), 1, r"^error: .*empty\.dat has no rows; check")
    short = work / "short.dat"
    short.write_text(HEADER + "0.0 1.0 11.44 1.47\n", encoding="utf-8")
    expect(run(base, short, 1e-13), 1, r"^error: .*short\.dat: row 1 has 4 numbers, not 5; give")
    expect(run(base, write(work, "two", ROWS[:2]), 1e-13), 1, r"^DIFFERENT: 3 rows in .*base\.dat, 2 in .*two\.dat$")


def arguments_case(work):
    base = write(work, "base", ROWS)
    for tolerance in ("abc", "nan", "inf", "-1"):
        expect(run(base, base, tolerance), 2, r"^usage: compare_averages\.py ")
    expect(run(base, base), 2, r"^usage: compare_averages\.py ")


def main(argv):
    """Runs the named case in WORK."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", metavar="CASE", choices=CASES, help="the case to run")
    parser.add_argument("work", metavar="WORK", type=Path, help="scratch directory (recreated)")
    args = parser.parse_args(argv)
    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    {"equal": equal_case, "non-finite": non_finite_case, "rows": rows_case, "arguments": arguments_case}[args.case](
        args.work
    )
    print(f"PASS: {args.case}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
