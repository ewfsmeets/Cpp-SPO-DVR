"""Tests of spodvr_compare_wave (tools/compare/compare_wave.cpp; ctest CompareWave.*).

Each test writes small wave-packet files (raw doubles, re/im pairs, as spodvr_dump_wave writes
them) into its scratch directory `work`, runs the tool (`tool`: option --compare-wave-tool,
conftest.py) on them and checks its exit status and output line:
  test_finite      identical files pass ("bitwise identical"); a difference of 1e-14 relative fails
                   with the default tolerance 0 and passes with 1e-13; files of different sizes fail
  test_non_finite  a NaN in b, an inf in a, a NaN at the same element of both (bitwise equal files)
                   and a NaN with tolerance 1e300: each fails, naming the file and the element
  test_arguments   one file, four arguments, a tolerance that is not a number, negative, NaN or
                   inf: exit 2 with a usage line; a missing, an empty and a 12-byte file: exit 1
                   with "error: <what>; <fix>"
Usage: uv run python -I -B -m pytest -c tools/compare/pytest.ini --rootdir=tools/compare \
           tools/compare/test_compare_wave.py --compare-wave-tool=build/tools/spodvr_compare_wave
Environment: none read. Run time: well under a second.
Docs: docs/performance.md, section "Benchmarking checklist" (the final wave packet of two builds).
"""

import math
import struct
import subprocess

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


def run_with(tool, *arguments):
    """Runs the tool with the given arguments; returns (exit status, stderr)."""
    result = subprocess.run([str(tool), *map(str, arguments)], capture_output=True, text=True, check=False)
    return result.returncode, result.stderr.strip()


def expect(what, status, output, want_status, want_output):
    """Raises if the exit status differs or want_output is not in the output."""
    if status != want_status or want_output not in output:
        raise AssertionError(
            f"{what}: exit {status}, output {output!r}; expected exit {want_status} and {want_output!r}"
        )


def test_finite(tool, work):
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


def test_non_finite(tool, work):
    a = write(work / "a.bin", WAVE)
    with_nan = list(WAVE)
    with_nan[3] = math.nan
    b = write(work / "b.bin", with_nan)
    expect(
        "NaN in b",
        *run(tool, a, b),
        1,
        f"DIFFERENT: non-finite value nan in {b} at element 3 (complex value 1, imaginary part)",
    )
    expect("NaN in b, tolerance 1e300", *run(tool, a, b, "1e300"), 1, "DIFFERENT: non-finite value nan")
    with_inf = list(WAVE)
    with_inf[0] = math.inf
    c = write(work / "c.bin", with_inf)
    expect(
        "inf in a",
        *run(tool, c, a),
        1,
        f"DIFFERENT: non-finite value inf in {c} at element 0 (complex value 0, real part)",
    )
    expect(
        "NaN at the same element of both", *run(tool, b, b), 1, f"DIFFERENT: non-finite value nan in {b} at element 3"
    )


def test_arguments(tool, work):
    """Bad command lines exit 2 with the usage line; files that cannot be compared are refused."""
    a = write(work / "a.bin", WAVE)
    for what, arguments in (
        ("one file", [a]),
        ("four arguments", [a, a, "0", "x"]),
        ("tolerance abc", [a, a, "abc"]),
        ("tolerance 1e-13x", [a, a, "1e-13x"]),
        ("tolerance -1e-13", [a, a, "-1e-13"]),
        ("tolerance nan", [a, a, "nan"]),
        ("tolerance inf", [a, a, "inf"]),
    ):
        expect(what, *run_with(tool, *arguments), 2, "usage: spodvr_compare_wave <a> <b> [tolerance]")
    empty = write(work / "empty.bin", [])
    partial = work / "partial.bin"
    partial.write_bytes(struct.pack("=3i", 1, 2, 3))  # 12 bytes: not a whole complex value
    for what, other in (("missing file", work / "missing.bin"), ("empty file", empty), ("12-byte file", partial)):
        status, stderr = run_with(tool, a, other)
        expect(what, status, stderr, 1, f"'{other}'")
        if not stderr.startswith("error: ") or "; give two files written by spodvr_dump_wave" not in stderr:
            raise AssertionError(f'{what}: {stderr!r} is not "error: <what>; <fix>"')
