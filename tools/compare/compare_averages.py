#!/usr/bin/env python3
"""Compares two averages.dat files of spodvr (docs/output_files.md) column by column: the largest
relative difference |a - b| / max(|a|, |b|) of each column must not exceed the tolerance. A value
that is not finite (NaN, inf) fails, in either file and also where both files hold it; so do
different row counts. Run by tools/compare_propagation.sh.

Usage: python3 -I -B tools/compare/compare_averages.py BASE NEW TOLERANCE
  BASE, NEW  averages.dat files;  TOLERANCE  a finite number >= 0 (0: equal values only)
Output: the largest relative difference per column, then "OK: ..." or "DIFFERENT: ...".
Exit status: 0 within the tolerance, 1 different, or refused ("error: <what>; <fix>": a file
without rows, unreadable, or not five numbers per row), 2 bad arguments.
Environment: none read. Run time: well under a second.
Docs: README.md, section "Testing".
"""

import argparse
import math
import sys
from pathlib import Path

COLUMNS = ("time", "norm", "<Z>", "<r>", "rest norm")  # t, norm, <Z>, <r>, rest norm
TINY = 1e-300  # the denominator floor: two zeros differ by 0


class Refusal(Exception):
    """A file that cannot be compared, printed as "error: <what>; <fix>" with exit status 1."""


def parse_tolerance(text):
    """The tolerance argument: a finite number >= 0 (argparse type; tools/compare_propagation.sh
    also checks it with this function before its builds)."""
    value = float(text)
    if not (math.isfinite(value) and value >= 0):
        raise argparse.ArgumentTypeError(f"{text!r} is not a finite number >= 0")
    return value


def read_rows(path):
    """The rows of an averages.dat ('#' lines skipped), five numbers each."""
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        rows = [[float(x) for x in line.split()] for line in lines if line.strip() and not line.startswith("#")]
    except (OSError, ValueError) as error:
        raise Refusal(f"{path}: {error}; give two averages.dat files written by spodvr") from error
    for number, row in enumerate(rows, start=1):
        if len(row) != len(COLUMNS):
            raise Refusal(f"{path}: row {number} has {len(row)} numbers, not {len(COLUMNS)}; give an averages.dat")
    if not rows:
        raise Refusal(f"{path} has no rows; check that the run wrote its averages")
    return rows


def first_non_finite(rows):
    """(row number, column name, value) of the first value that is not finite, or None."""
    for number, row in enumerate(rows, start=1):
        for name, value in zip(COLUMNS, row):
            if not math.isfinite(value):
                return number, name, value
    return None


def largest_differences(base, new):
    """The largest relative difference of each column (finite values only, so never NaN)."""
    return [max(abs(a[k] - b[k]) / max(abs(a[k]), abs(b[k]), TINY) for a, b in zip(base, new))
            for k in range(len(COLUMNS))]


def compare(base_path, new_path, tolerance):
    """Prints the comparison; returns the exit status (0 within the tolerance, 1 different)."""
    base, new = read_rows(base_path), read_rows(new_path)
    if len(base) != len(new):
        print(f"DIFFERENT: {len(base)} rows in {base_path}, {len(new)} in {new_path}")
        return 1
    for path, rows in ((base_path, base), (new_path, new)):
        found = first_non_finite(rows)
        if found:
            print(f"DIFFERENT: {path}: row {found[0]}, {found[1]}: {found[2]} is not finite")
            return 1
    differences = largest_differences(base, new)
    for name, difference in zip(COLUMNS, differences):
        print(f"  {name:10s} max relative difference {difference:.3e}")
    worst = max(differences)
    if not worst <= tolerance:  # written so that a NaN would fail too
        print(f"DIFFERENT: {worst:.3e} > tolerance {tolerance:g}")
        return 1
    print(f"OK: within {tolerance:g} of the base")
    return 0


def main(argv):
    """Compares BASE and NEW; a refusal exits with status 1."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("base", metavar="BASE", help="averages.dat of the base build")
    parser.add_argument("new", metavar="NEW", help="averages.dat of the new build")
    parser.add_argument("tolerance", metavar="TOLERANCE", type=parse_tolerance, help="largest relative difference")
    args = parser.parse_args(argv)
    try:
        return compare(args.base, args.new, args.tolerance)
    except Refusal as refusal:
        sys.stdout.flush()
        print(f"error: {refusal}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
