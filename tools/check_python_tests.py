"""Checks that every pytest test of the repository is run by exactly one ctest (ctest Python.EveryTestRegistered).

tests/CMakeLists.txt writes one line per pytest ctest into <build>/pytest_selections.txt: the tool directory,
the selection (a test module, or module::test) and whether that ctest exists in this configuration ("built")
or only in others ("not-built", e.g. CompareWave.* without SPODVR_BUILD_TOOLS). This script collects each
tool's tests with pytest and requires every collected test to be selected exactly once, and every selection
to select at least one test.

--plant-faults (ctest Python.UnregisteredTestFails) is the negative test: it drops the first selection and
repeats the last one, in memory, and passes only if exactly the tests of those two selections are reported.

Usage: uv run python -I -B tools/check_python_tests.py <pytest_selections.txt> --root <repository> [--plant-faults]
Exit status: 0 every test is selected once (with --plant-faults: exactly the planted faults are found),
1 otherwise or pytest could not collect a tool's tests, 2 bad arguments.
Environment: none read (pytest's collection imports the test modules; ctest sets MPLBACKEND=Agg).
Run time: a few seconds (one pytest collection per tool). Needs: pytest in the Python that runs it.
Docs: docs/CONTRIBUTING.md, section "Testing".
"""

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Selection:
    """One line of pytest_selections.txt."""

    tool: str  # tool directory, relative to the repository
    select: str  # test module or module::test, relative to the tool directory
    built: bool  # the ctest exists in this configuration


def read_selections(path: Path) -> list[Selection]:
    """The selections of pytest_selections.txt (tab-separated; '#' lines are comments)."""
    selections = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 3 or fields[2] not in ("built", "not-built"):
            raise SystemExit(f"error: {path}:{number}: expected tool<TAB>selection<TAB>built|not-built; reconfigure")
        selections.append(Selection(fields[0], fields[1], fields[2] == "built"))
    return selections


def collected_tests(root: Path, tool: str) -> list[str]:
    """The node ids pytest collects for a tool (relative to its directory), with the tool's pytest.ini."""
    directory = root / tool
    command = [sys.executable, "-I", "-B", "-m", "pytest", "-c", str(directory / "pytest.ini")]
    command += ["--rootdir", str(directory), "-p", "no:cacheprovider", "--collect-only", "-q", str(directory)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        output = result.stdout + result.stderr
        raise SystemExit(f"error: pytest could not collect the tests of {tool} (exit {result.returncode}):\n{output}")
    return [line.strip() for line in result.stdout.splitlines() if "::" in line]


def selects(selection: str, node: str) -> bool:
    """Whether a selection (module or module::test) selects a node id (also a parametrised test_x[case])."""
    return node == selection or node.startswith((selection + "::", selection + "["))


def problems(selections: list[Selection], tests: dict[str, list[str]]) -> list[tuple[str, str, str]]:
    """(tool, test or selection, what is wrong) for every test not selected exactly once and every selection
    that selects nothing."""
    found = []
    for tool, nodes in sorted(tests.items()):
        for node in nodes:
            count = sum(1 for s in selections if s.tool == tool and selects(s.select, node))
            if count != 1:
                found.append((tool, node, f"selected by {count} ctests, not 1"))
    for s in selections:
        if not any(selects(s.select, node) for node in tests.get(s.tool, [])):
            found.append((s.tool, s.select, "selects no test"))
    return found


def planted_faults(selections: list[Selection], tests: dict[str, list[str]]) -> tuple[list, set[tuple[str, str]]]:
    """The problems found after dropping the first selection and repeating the last, and the (tool, test)
    pairs that those two faults must make reported."""
    dropped, repeated = selections[0], selections[-1]
    found = problems(selections[1:] + [repeated], tests)
    expected = {(s.tool, node) for s in (dropped, repeated) for node in tests[s.tool] if selects(s.select, node)}
    return found, expected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("selections", type=Path, help="<build>/pytest_selections.txt")
    parser.add_argument("--root", type=Path, required=True, help="the repository")
    parser.add_argument("--plant-faults", action="store_true", help="the negative test")
    args = parser.parse_args()
    selections = read_selections(args.selections)
    if len(selections) < 2:
        print(f"error: {args.selections} lists {len(selections)} selections; reconfigure the build")
        return 1
    tests = {tool: collected_tests(args.root, tool) for tool in sorted({s.tool for s in selections})}
    total = sum(len(nodes) for nodes in tests.values())
    if not args.plant_faults:
        found = problems(selections, tests)
        for tool, item, what in found:
            print(f"NOT REGISTERED ONCE: {tool}/{item}: {what}")
        if found:
            return 1
        built = sum(1 for s in selections if s.built)
        print(
            f"OK: {total} pytest tests of {len(tests)} tools, each selected exactly once: {built} ctests here, "
            f"{len(selections) - built} only in other configurations"
        )
        return 0
    found, expected = planted_faults(selections, tests)
    for tool, item, what in found:
        print(f"planted: {tool}/{item}: {what}")
    if {(tool, item) for tool, item, _ in found} != expected or len(found) != len(expected):
        print(f"FAIL: expected exactly the {len(expected)} tests of the two planted faults to be reported")
        return 1
    print(f"OK: the {len(expected)} tests of the dropped and the repeated selection are reported, nothing else")
    return 0


if __name__ == "__main__":
    sys.exit(main())
