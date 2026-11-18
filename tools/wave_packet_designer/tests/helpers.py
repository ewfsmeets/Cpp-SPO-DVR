"""Shared helpers of the test cases: the repository paths, the inputs, and the checks."""

import json
from dataclasses import replace
from pathlib import Path

from program.config import InputParameters
from program.findings import Finding, Severity
from program.json_input import parse_input

TOOL_DIR = Path(__file__).resolve().parent.parent
REPOSITORY = TOOL_DIR.parent.parent
DZ_INPUT = REPOSITORY / "inputs/thesis_fig5_3_dz/input.json"
# The program's values for DZ_INPUT that the tests compare with come from the set-up output of
#   build-cuda/spodvr inputs/thesis_fig5_3_dz/input.json --device cuda --threads 8      (2026-11-08)
# (blocks INPUT, GRIDS and ROTATIONAL BASIS AND INITIAL WAVE PACKET). All are printed before
# the propagation starts, so the run can be stopped there; the packet is built on the host, so a CPU run
# prints the same values.
ER_INPUT = REPOSITORY / "inputs/thesis_fig5_3_er/input.json"
HARDWALL_INPUT = REPOSITORY / "tests/reference/hardwall_small/input.json"


class TestFailure(Exception):
    """A failed check; the runner prints its message."""


def check(condition: bool, message: str) -> None:
    """Fails the test with `message` unless `condition` holds."""
    if not condition:
        raise TestFailure(message)


def check_equal(actual: object, expected: object, what: str) -> None:
    """Fails unless actual == expected (exact, also for floats)."""
    check(actual == expected, f"{what}: expected {expected!r}, got {actual!r}")


def load_values(path: Path) -> dict:
    """The JSON object of an input file."""
    return json.loads(path.read_text())


def load_parameters(path: Path, **changes: object) -> InputParameters:
    """The parameters of an input file (which must be readable), with `changes` to InputParameters fields."""
    parsed = parse_input(load_values(path))
    check(parsed.parameters is not None, f"{path.name} is not readable: {parsed.findings}")
    return replace(parsed.parameters, **changes)


def errors(findings: list[Finding] | tuple[Finding, ...]) -> list[Finding]:
    return [f for f in findings if f.severity is Severity.ERROR]


def check_has_error(findings: list[Finding] | tuple[Finding, ...], key: str, message_part: str, what: str) -> None:
    """Fails unless there is at least one error for `key` whose message contains `message_part`."""
    matching = [f for f in errors(findings) if f.key == key and message_part in f.message]
    check(len(matching) >= 1, f"{what}: no error {key}: '{message_part}' in {[(f.key, f.message) for f in findings]}")
