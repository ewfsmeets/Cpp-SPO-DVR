"""Tests (ctest WavePacketDesigner.Cli): the `check` command on the repository's inputs, and its exit codes
(wave_packet_designer.py)."""

import json
import re
import subprocess
import sys
from pathlib import Path

from tests.helpers import DZ_INPUT, REPOSITORY, TOOL_DIR, check, check_equal

CLI = TOOL_DIR / "wave_packet_designer.py"
# Z_asymptotic of the SRP48 data set of inputs/scm_example: the last Z of Molecule/extrapol.dat, 6.5 angstrom,
# converted like the CRP library does (x * (1 / 0.529177249)); the data set is C3v.
SRP48_Z_ASYMPTOTIC = 6.5 * (1.0 / 0.529177249)
NUMBER = r" [0-9.]+(?:e[+-][0-9]+)?"  # one value of a table row, with the space before it


def run_check(*arguments: str) -> subprocess.CompletedProcess:
    return run_cli("check", *arguments)


def run_cli(*arguments: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", "-B", str(CLI), *arguments], capture_output=True, text=True)


def test_repository_inputs_are_accepted() -> None:
    """The four inputs of inputs/ pass (exit 0); the SCM example with its data set's Z_asymptotic and group."""
    for name in ("thesis_fig5_3_dz", "thesis_fig5_3_er", "hardwall_small"):
        result = run_check(str(REPOSITORY / "inputs" / name / "input.json"))
        check(result.returncode == 0 and "spodvr accepts" in result.stdout, f"{name}: {result.stdout}{result.stderr}")
    result = run_check(
        str(REPOSITORY / "inputs/scm_example/input.json"), "--zasymp", repr(SRP48_Z_ASYMPTOTIC), "--point-group", "C3v"
    )
    check(result.returncode == 0, f"scm_example: {result.stdout}{result.stderr}")


def test_check_prints_r_t_time_and_tails() -> None:
    """`check` on the DZ input prints R and T of each absorber besides R + T, the TIME section (the slowest
    analysed part's a.u. per bohr and the times it needs) and the packet's two momentum tails: the quantities
    inputs/README.md and inputs/thesis_fig5_3_dz/README.md quote."""
    result = run_check(str(DZ_INPUT))
    check(result.returncode == 0, f"{result.returncode}: {result.stdout}{result.stderr}")
    text = " ".join(result.stdout.split())  # whitespace collapsed: no dependence on columns or line wrapping

    def table_row(absorber: str, label: str) -> list[str]:
        """The values of one row (R + T, R or T) of an absorber in ABSORBERS: REFLECTION + TRANSMISSION, each row
        found by its label after the absorber's energy row."""
        rows = re.search(
            rf"\b{absorber} E \(eV\)(?:{NUMBER})+ R \+ T(?P<total>(?:{NUMBER})+) R(?P<R>(?:{NUMBER})+) "
            rf"T(?P<T>(?:{NUMBER})+)",
            text,
        )
        check(rows is not None, f"{absorber}: no E, R + T, R, T rows: {result.stdout}")
        return rows.group({"R + T": "total", "R": "R", "T": "T"}[label]).split()

    check_equal(table_row("ZOpt", "R")[0], "4.2e-05", "ZOpt R at the lowest energy")
    check_equal(table_row("ZOpt", "T")[-1], "8.3e-06", "ZOpt T at the highest energy")
    check_equal(table_row("ZSpOpt", "R + T")[-1], "2.7e-09", "ZSpOpt R + T at the highest energy")
    for pattern in (
        re.escape("TIME (inputs/README.md section 8)"),
        r"M / k_min \.+ 354 a\.u\. per bohr",
        re.escape("after 10252 a.u."),
        re.escape("needs 22438"),
        r"Packet moving outward, k > 0 \.+ 2\.9e-21",
        r"Wrapped to outward motion, beyond k_grid \.+ 1\.4e-20",
    ):
        check(re.search(pattern, text) is not None, f"{pattern!r} not found: {result.stdout}")


def test_rejected_input_exits_with_1(work: Path) -> None:
    """Z0 too close to Zana: the report names the TolSp error and the exit status is 1."""
    values = json.loads(DZ_INPUT.read_text()) | {"Z0": 16.78}
    path = work / "close.json"
    path.write_text(json.dumps(values))
    result = run_check(str(path))
    check(result.returncode == 1 and "TolSp" in result.stdout, f"{result.returncode}: {result.stdout}")


def test_data_set_system_needs_its_values() -> None:
    """A CRP data set without --zasymp and --point-group: a usage problem (exit 2, not "spodvr stops"), and
    stderr names both options."""
    result = run_check(str(REPOSITORY / "inputs/scm_example/input.json"))
    check_equal(result.returncode, 2, "exit status")
    check("--zasymp" in result.stderr and "--point-group" in result.stderr, result.stderr)
    check(result.stdout == "", f"no report expected: {result.stdout}")


def test_memory_of_the_dz_input() -> None:
    """`memory` with the DZ run's symmetry counts prints the program's totals (7.9 and 7.3 GiB on the CPU,
    8.8 GiB on the GPU; no other row has these sizes)."""
    result = run_cli("memory", str(DZ_INPUT), "--representatives", "19")
    check(result.returncode == 0, f"{result.returncode}: {result.stdout}{result.stderr}")
    for total in ("7.9 GiB", "7.3 GiB", "8.8 GiB"):
        check(total in result.stdout, f"{total} missing: {result.stdout}")


def test_impossible_memory_options_exit_with_2() -> None:
    """Usage errors (exit 2): 0 threads, more representatives than (x, y) points (DZ has 144), more operations
    than the point group has (C6v: 12), and --operations without --representatives."""
    for options in (
        ("--threads", "0"),
        ("--representatives", "1000"),
        ("--representatives", "19", "--operations", "13"),
        ("--operations", "12"),
    ):
        check_equal(run_cli("memory", str(DZ_INPUT), *options).returncode, 2, f"exit status of {options}")


def test_suggest_writes_an_accepted_input(work: Path) -> None:
    """`suggest --output` writes a copy of the input with the suggested Z keys, which `check` accepts; it
    refuses to replace an existing file (exit 2)."""
    target = work / "suggested.json"
    target.unlink(missing_ok=True)
    result = run_cli("suggest", str(DZ_INPUT), "--dz", "0.12", "--output", str(target))
    check(result.returncode == 0 and target.exists(), f"{result.returncode}: {result.stdout}{result.stderr}")
    check_equal(run_check(str(target)).returncode, 0, "exit status of check on the copy")
    check_equal(run_cli("suggest", str(DZ_INPUT), "--output", str(target)).returncode, 2, "existing file")


def test_missing_file_exits_with_2(work: Path) -> None:
    """An unreadable file is a usage problem (exit 2), reported on stderr."""
    result = run_check(str(work / "missing.json"))
    check_equal(result.returncode, 2, "exit status")
    check("cannot read" in result.stderr, result.stderr)
