"""Tests (ctest WavePacketDesigner.RealCompare): the tolerant comparisons and NINT (program/real_compare.py)."""

from program.real_compare import nint, req, rge, rgt, rle, rlt, seq
from tests.helpers import check, check_equal


def test_tolerance_is_absolute_and_exclusive() -> None:
    """A difference of 1e-12 is neither equal nor greater (abs(d) < tol, d > tol); 9e-13 counts as equal."""
    check(not req(1e-12, 0.0) and not rgt(1e-12, 0.0) and not rge(1e-12, 0.0), "difference of exactly 1e-12")
    check(req(9e-13, 0.0) and rge(9e-13, 0.0) and rle(9e-13, 0.0), "difference of 9e-13")
    check(rge(15.0, 15.0 + 5e-13) and not rlt(15.0, 15.0 + 5e-13), "within the tolerance below")
    check(req(1e6, 1e6 + 1e-10) is False, "the tolerance does not scale with the values")


def test_nint_rounds_half_away_from_zero() -> None:
    """std::lround, not Python's round(): 2.5 -> 3, -2.5 -> -3, and 0.49999999999999994 -> 0."""
    for x, expected in ((0.5, 1.0), (1.5, 2.0), (2.5, 3.0), (-2.5, -3.0), (0.49999999999999994, 0.0), (-0.4, -0.0)):
        check_equal(nint(x), expected, f"nint({x!r})")


def test_string_comparison() -> None:
    """.SEQ. ignores case and blanks at both ends."""
    check(seq("  h2cu111dailight ", "H2Cu111DaiLight"), "case and blanks")
    check(not seq("H2Cu111Dai Light", "H2Cu111DaiLight"), "inner blanks count")
