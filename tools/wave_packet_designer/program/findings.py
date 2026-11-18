"""Errors and warnings the program would give for an input, as the designer reports them.

Shadows: src/core/error.hpp (InputError, NotSupportedError)

The program stops at its first error (exit code 2 for InputError, 3 for NotSupportedError). The designer
continues within the same stage and lists every error of that stage, so one look shows all of them.
"""

from dataclasses import dataclass
from enum import Enum


class Severity(Enum):
    """How the program reacts to a finding."""

    ERROR = "error"  # the program stops
    WARNING = "warning"  # the program prints a warning and continues


@dataclass(frozen=True)
class Finding:
    """One error or warning of the program."""

    severity: Severity
    key: str  # the input key the message names ("" if none)
    message: str  # the program's message (wording may be shortened)
    origin: str  # the C++ function that raises it, e.g. "validate_input"


class FindingList:
    """Collects the findings of one stage, in the order the C++ checks them."""

    def __init__(self, origin: str) -> None:
        self.origin = origin
        self.findings: list[Finding] = []

    def error(self, key: str, message: str) -> None:
        """Records an error (where the C++ throws InputError)."""
        self.findings.append(Finding(Severity.ERROR, key, message, self.origin))

    def error_if(self, violated: bool, key: str, message: str) -> None:
        """Records an error when `violated` (the C++ `check(violated, key, message)`)."""
        if violated:
            self.error(key, message)

    def warning(self, key: str, message: str) -> None:
        """Records a warning."""
        self.findings.append(Finding(Severity.WARNING, key, message, self.origin))


def has_errors(findings: list[Finding]) -> bool:
    """True if any finding stops the program."""
    return any(f.severity is Severity.ERROR for f in findings)
