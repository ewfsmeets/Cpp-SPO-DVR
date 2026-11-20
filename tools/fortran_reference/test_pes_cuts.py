"""Tests of pes_cuts.py (ctest FortranReference.PesCuts*), run as a program on files in the scratch directory
`work`:
  test_count_and_lines  the number of cuts and each cut's DVR.inp lines (the shortest decimal of each value; an
                        empty cut gives none); the committed pes_cuts.json of every case is read
  test_refusals         not JSON, no list "cuts", a key that is not a coordinate, a string, a bool and a NaN
                        value, a cut number out of range, a DVR.inp that sets an appended key (in any case):
                        exit 1 with "error: <what>; <fix>"; a bad command line: exit 2

Usage: uv run python -I -B -m pytest -c tools/fortran_reference/pytest.ini \
           --rootdir=tools/fortran_reference tools/fortran_reference/test_pes_cuts.py
Environment: none read. Run time: under a second. Docs: tools/fortran_reference/README.md.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

TOOLS_DIR = Path(__file__).resolve().parent
PES_CUTS = TOOLS_DIR / "pes_cuts.py"
REFERENCE = TOOLS_DIR.parents[1] / "tests" / "reference"


def run(*arguments: object) -> tuple[int, str, str]:
    """(exit status, stdout, stderr) of pes_cuts.py with arguments."""
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(PES_CUTS), *map(str, arguments)], capture_output=True, text=True
    )
    return result.returncode, result.stdout, result.stderr


def cuts_file(work: Path, text: str) -> Path:
    """A pes_cuts.json with the given text."""
    path = work / "pes_cuts.json"
    path.write_text(text, encoding="utf-8")
    return path


def test_count_and_lines(work: Path) -> None:
    """Two cuts: the first empty (no line: the Fortran defaults), the second with the shortest decimals of its
    values (an integer becomes a real); every committed pes_cuts.json reads."""
    path = cuts_file(work, json.dumps({"cuts": [{}, {"x2DPES": 1.452, "theta2DPES": 37, "phi2DPES": 0.1}]}))
    assert run("count", path) == (0, "2\n", "")
    assert run("lines", path, 1) == (0, "\n", "")
    expected = "x2DPES : REAL : 1.452;\ntheta2DPES : REAL : 37.0;\nphi2DPES : REAL : 0.1;\n"
    assert run("lines", path, 2) == (0, expected, "")
    committed = sorted(REFERENCE.glob("*/pes_cuts.json"))
    assert committed, f"no pes_cuts.json under {REFERENCE}"
    for file in committed:
        status, out, err = run("count", file)
        assert status == 0 and int(out) >= 1, f"{file}: {err}"


def expect_refusal(pattern: str, *arguments: object) -> None:
    """pes_cuts.py exits 1 with "error: <what>; <fix>" matching pattern."""
    status, _, err = run(*arguments)
    if status != 1 or not re.search(rf"^error: .*{pattern}.*; [a-z]", err, re.MULTILINE):
        pytest.fail(f"should refuse with {pattern!r} (exit {status}): {err}")


def test_refusals(work: Path) -> None:
    """Each bad file, cut number and DVR.inp refuses; a bad command line exits 2."""
    cases = [
        ('{"cuts": ', "is not JSON with a list"),
        ('{"cut": []}', "is not JSON with a list"),
        ('{"cuts": [1]}', "is not a list of objects"),
        ('{"cuts": [{"z2DPES": 1.0}]}', "cut 1: z2DPES is not a write2DPES coordinate"),
        ('{"cuts": [{"x2DPES": "1.0"}]}', "cut 1: x2DPES = '1.0' is not a finite number"),
        ('{"cuts": [{"x2DPES": true}]}', "cut 1: x2DPES = True is not a finite number"),
        ('{"cuts": [{"phi2DPES": NaN}]}', "cut 1: phi2DPES = nan is not a finite number"),
    ]
    for text, pattern in cases:
        expect_refusal(re.escape(pattern), "count", cuts_file(work, text))
    path = cuts_file(work, '{"cuts": [{}]}')
    expect_refusal("cut 2 is not one of 1..1", "lines", path, 2)
    expect_refusal("cut x is not one of 1..1", "lines", path, "x")
    dvr = work / "DVR.inp"
    dvr.write_text('Title : "x";\n  THETA2DPES : REAL : 1.0;\n', encoding="utf-8")
    expect_refusal(r"DVR\.inp:2 sets THETA2DPES itself", "check-dvr", dvr)
    dvr.write_text('Title : "theta2DPES";\n', encoding="utf-8")
    assert run("check-dvr", dvr) == (0, "", "")
    assert run("count")[0] == 2
