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
