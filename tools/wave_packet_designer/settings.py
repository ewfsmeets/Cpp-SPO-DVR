"""What the designer needs besides the input file: command-line options and front-end fields.

Part of the data layer (README, "Architecture"). Values the program takes from a CRP data set have no
built-in value here: the program's set-up output prints them. this_machine() finds the memory of the machine
the front end runs on, for the memory checks.
"""

import os
import shutil
import subprocess
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
    # The memory the run may use (bytes): the physical RAM, and the free GPU memory (what the program checks
    # with --device cuda). None: unknown, no memory check.
    ram_bytes: int | None = None
    gpu_bytes: int | None = None
    # The time budget (inputs/README.md section 8): the turning point Z_turn (bohr) and the lowest outgoing
    # normal energy E_n,min (eV) tEnd must converge. None: the defaults of design/rules.py (2 bohr, 0.03 eV).
    z_turn: float | None = None
    lowest_outgoing_energy: float | None = None


def this_machine() -> tuple[int | None, int | None]:
    """The physical RAM of this machine and the free memory of its first GPU (nvidia-smi), in bytes; None
    where unknown."""
    try:
        ram = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (OSError, ValueError):
        ram = None
    gpu = None
    if shutil.which("nvidia-smi"):
        query = ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"]
        result = subprocess.run(query, capture_output=True, text=True, timeout=10)
        lines = result.stdout.split()
        if result.returncode == 0 and lines and lines[0].isdigit():
            gpu = int(lines[0]) * 2**20  # MiB
    return ram, gpu
