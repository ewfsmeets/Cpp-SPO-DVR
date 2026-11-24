"""The Fortran's tolerant real comparisons and NINT, used wherever the program uses them.

Shadows: src/core/real_compare.hpp (fortran_compare::req, rgt, rge, seq), src/input/validation_helpers.hpp (nint)

The tolerance is absolute (1e-12), not relative. Some program checks use plain < and > instead; the
calling module says which one it copies. At a difference of exactly 1e-12, `rge` is false.
"""

import math

TOLERANCE = 1.0e-12  # fortran_compare::tolerance


def req(a: float, b: float) -> bool:
    """a .REQ. b: equal within the absolute tolerance."""
    return abs(a - b) < TOLERANCE


def rne(a: float, b: float) -> bool:
    """a .RNE. b."""
    return not req(a, b)


def rgt(a: float, b: float) -> bool:
    """a .RGT. b: greater by more than the tolerance."""
    return (a - b) > TOLERANCE


def rge(a: float, b: float) -> bool:
    """a .RGE. b: equal within the tolerance, or greater."""
    difference = a - b
    return abs(difference) < TOLERANCE or difference > TOLERANCE


def rle(a: float, b: float) -> bool:
    """a .RLE. b."""
    return not rgt(a, b)


def rlt(a: float, b: float) -> bool:
    """a .RLT. b."""
    return not rge(a, b)


def nint(x: float) -> float:
    """Fortran NINT as a double: rounds half away from zero (std::lround), unlike Python's round().

    The fraction is taken exactly (|x| - floor(|x|) has no rounding error); floor(|x| + 0.5) would round
    0.49999999999999994 up to 1.
    """
    whole = math.floor(abs(x))
    fraction = abs(x) - whole
    rounded = whole + 1 if fraction >= 0.5 else whole
    return math.copysign(float(rounded), x)


def seq(a: str, b: str) -> bool:
    """a .SEQ. b: equal ignoring case and leading/trailing blanks."""
    return a.strip(" ").lower() == b.strip(" ").lower()
