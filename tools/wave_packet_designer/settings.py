"""What the designer needs besides the input file: command-line options and front-end fields.

Part of the data layer (README, "Architecture"). Values the program takes from a CRP data set have no
built-in value here: the program's set-up output prints them.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Settings:
    """Values that are not input keys."""

    z_asymptotic: float | None = None  # bohr; needed for CRP data sets (SystemID CRPlibrary, CRPDistrtSurface)
    point_group: str | None = None  # likewise (C2v, C4v, C6v, C3v, C3 or none)
    threads: int = 8  # OpenMP threads of the run: each has its own work arrays (the repository runs with 8)
    # With UseSymmetry: the (x, y) points at which the program stores the potential, and the operations it
    # uses (its SURFACE SYMMETRY block prints both). Without them the memory estimate assumes no reduction.
    representatives: int | None = None
    operations: int | None = None  # None: the order of the point group (all its operations usable)
