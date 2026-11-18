"""Case cli: the `check` command on the repository's inputs, and its exit codes (wave_packet_designer.py)."""

import json
import subprocess
import sys
from pathlib import Path

from tests.helpers import DZ_INPUT, REPOSITORY, TOOL_DIR, check, check_equal

CLI = TOOL_DIR / "wave_packet_designer.py"
# Z_asymptotic of the SRP48 data set of inputs/scm_example: the last Z of Molecule/extrapol.dat, 6.5 angstrom,
# converted like the CRP library does (x * (1 / 0.529177249)); the data set is C3v.
SRP48_Z_ASYMPTOTIC = 6.5 * (1.0 / 0.529177249)


def run_check(*arguments: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", "-B", str(CLI), "check", *arguments], capture_output=True, text=True)


def test_repository_inputs_are_accepted(work: Path) -> None:
    """The four inputs of inputs/ pass (exit 0); the SCM example with its data set's Z_asymptotic and group."""
    for name in ("thesis_fig5_3_dz", "thesis_fig5_3_er", "hardwall_small"):
        result = run_check(str(REPOSITORY / "inputs" / name / "input.json"))
        check(result.returncode == 0 and "spodvr accepts" in result.stdout, f"{name}: {result.stdout}{result.stderr}")
    result = run_check(
        str(REPOSITORY / "inputs/scm_example/input.json"), "--zasymp", repr(SRP48_Z_ASYMPTOTIC), "--point-group", "C3v"
    )
    check(result.returncode == 0, f"scm_example: {result.stdout}{result.stderr}")


def test_rejected_input_exits_with_1(work: Path) -> None:
    """Z0 too close to Zana: the report names the TolSp error and the exit status is 1."""
    values = json.loads(DZ_INPUT.read_text()) | {"Z0": 16.78}
    path = work / "close.json"
    path.write_text(json.dumps(values))
    result = run_check(str(path))
    check(result.returncode == 1 and "TolSp" in result.stdout, f"{result.returncode}: {result.stdout}")


def test_data_set_system_needs_its_values(work: Path) -> None:
    """A CRP data set without --zasymp and --point-group: a usage problem (exit 2, not "spodvr stops"), and
    stderr names both options."""
    result = run_check(str(REPOSITORY / "inputs/scm_example/input.json"))
    check_equal(result.returncode, 2, "exit status")
    check("--zasymp" in result.stderr and "--point-group" in result.stderr, result.stderr)
    check(result.stdout == "", f"no report expected: {result.stdout}")


def test_missing_file_exits_with_2(work: Path) -> None:
    """An unreadable file is a usage problem (exit 2), reported on stderr."""
    result = run_check(str(work / "missing.json"))
    check_equal(result.returncode, 2, "exit status")
    check("cannot read" in result.stderr, result.stderr)


TESTS = [
    test_repository_inputs_are_accepted,
    test_rejected_input_exits_with_1,
    test_data_set_system_needs_its_values,
    test_missing_file_exits_with_2,
]
