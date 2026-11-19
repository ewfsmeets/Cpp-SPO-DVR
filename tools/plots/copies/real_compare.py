"""The program's tolerant comparison of reals.

Shadows: src/core/real_compare.hpp (rne)
"""

TOLERANCE = 1e-12  # absolute, as the Fortran's .REQ./.RNE. (MyOperators.f90)


def rne(a: float, b: float) -> bool:
    """True unless |a - b| < 1e-12, the program's absolute tolerance; so NaN is never equal (C++: !req(a, b))."""
    return not abs(a - b) < TOLERANCE
