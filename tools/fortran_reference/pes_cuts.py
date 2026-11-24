"""The write2DPES cuts of a reference case (tests/reference/<case>/pes_cuts.json) for the Fortran runs of
run_pes_cut_reference.sh and regenerate_all.sh: their number, the DVR.inp lines of one cut, and the check that the
case's DVR.inp sets none of the keys the script appends.

pes_cuts.json holds {"cuts": [{...}, ...]}: each object's keys are write2DPES coordinates (KEYS; x, y in bohr,
theta, phi in degrees), each optional (the Fortran default then applies). The C++ gate 1c reads the same file
(tests/test_regression_pes.cpp). A float is written as Python's repr, the shortest decimal that reads back to the
same double, so gfortran and the C++ get the same values.

Usage: python3 -I -B tools/fortran_reference/pes_cuts.py count FILE
       python3 -I -B tools/fortran_reference/pes_cuts.py lines FILE K     (K from 1)
       python3 -I -B tools/fortran_reference/pes_cuts.py check-dvr DVR_INP
Exit status: 0; 1 with "error: <what>; <fix>" for a file that is not JSON with a list "cuts" of objects, a key
that is not a coordinate, a value that is not a finite number, a K out of range, or a DVR.inp that sets
write2DPES or a coordinate itself; 2 for a bad command line.
Environment: none. Run time: instant. Needs: Python 3 (standard library). Docs: tools/fortran_reference/README.md.
"""

import json
import math
import sys

# The coordinates a cut may set (InputNew.f90:402-405; r1727: 446-449).
KEYS = ("x2DPES", "y2DPES", "theta2DPES", "phi2DPES")
APPENDED = ("write2DPES", *KEYS)  # what the script appends to DVR.inp: a second line of a key would be ambiguous


class Refusal(Exception):
    """A bad pes_cuts.json or argument: the message is "<what>; <fix>"."""


def read_cuts(path: str) -> list[dict]:
    """The cuts of a pes_cuts.json, each checked."""
    try:
        with open(path, encoding="utf-8") as file:
            cuts = json.load(file)["cuts"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise Refusal(f'{path} is not JSON with a list "cuts" ({error}); fix the file') from None
    if not isinstance(cuts, list) or not all(isinstance(cut, dict) for cut in cuts):
        raise Refusal(f'{path}: "cuts" is not a list of objects; fix the file')
    for k, cut in enumerate(cuts, start=1):
        for key, value in cut.items():
            if key not in KEYS:
                raise Refusal(f"{path}: cut {k}: {key} is not a write2DPES coordinate; use one of {', '.join(KEYS)}")
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise Refusal(f"{path}: cut {k}: {key} = {value!r} is not a finite number; give a number")
    return cuts


def dvr_lines(cut: dict) -> list[str]:
    """The DVR.inp lines of one cut: `Key : REAL : value;` (the Fortran reads these keys only with their type)."""
    return [f"{key} : REAL : {float(value)!r};" for key, value in cut.items()]


def check_dvr(path: str) -> None:
    """Refuses a DVR.inp with a line `Key :` of a key the script appends (Fortran key names ignore case)."""
    appended = {key.lower() for key in APPENDED}
    with open(path, encoding="utf-8") as file:
        for number, line in enumerate(file, start=1):
            name = line.split(":", 1)[0].strip()
            if ":" in line and name.lower() in appended:
                raise Refusal(f"{path}:{number} sets {name} itself; remove it (pes_cuts.json gives the cuts)")


def main(arguments: list[str]) -> int:
    """`count FILE` prints the number of cuts; `lines FILE K` prints the DVR.inp lines of cut K; `check-dvr
    DVR_INP` refuses a DVR.inp that sets an appended key."""
    if len(arguments) == 2 and arguments[0] == "check-dvr":
        check_dvr(arguments[1])
        return 0
    if len(arguments) == 2 and arguments[0] == "count":
        print(len(read_cuts(arguments[1])))
        return 0
    if len(arguments) == 3 and arguments[0] == "lines":
        cuts = read_cuts(arguments[1])
        if not arguments[2].isdigit() or not 1 <= int(arguments[2]) <= len(cuts):
            raise Refusal(f"cut {arguments[2]} is not one of 1..{len(cuts)} of {arguments[1]}; name a listed cut")
        print("\n".join(dvr_lines(cuts[int(arguments[2]) - 1])))
        return 0
    print("usage: pes_cuts.py count FILE | lines FILE K | check-dvr DVR_INP", file=sys.stderr)
    return 2


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Refusal as refusal:
        print(f"error: {refusal}", file=sys.stderr)
        sys.exit(1)
