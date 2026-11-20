"""pes_cut_2d.dat of a write2DPES run: V(r, Z) of the PES itself (no switching, no VMin/VMax) at one fixed x, y,
theta, phi, on the r grid and the whole scattering Z grid (docs/output_files.md, "Data files"): '#' lines, then one
line "r Z V(hartree) V(eV)" per point, Z fastest, a blank line after each r. Where the cut was taken and the grid
sizes come from the run's run_info.json (its "pes_cut" and "grids"), not from the file's '#' lines.

Shadows: src/report/pes_cut_file.cpp (write_pes_cut_file), src/report/run_info_file.cpp (pes_cut)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.errors import ReadError
from readers.run_info import RunInfo
from readers.tables import numbers, text_of

FILE = "pes_cut_2d.dat"


@dataclass
class PesCut:
    """The cut: V on the (Z, r) grid and where it was taken."""

    v_hartree: np.ndarray  # (nz, nr): V(Z, r), the PES itself
    z_bohr: np.ndarray  # (nz,): the scattering Z grid, all of it (also Z >= zasymp)
    r_bohr: np.ndarray  # (nr,)
    x_bohr: float  # along the first lattice vector
    y_bohr: float  # along the second lattice vector
    theta_deg: float
    phi_deg: float


def read_pes_cut(run: Path, info: RunInfo) -> PesCut:
    """The cut of a write2DPES run, its position and orientation from run_info.json; ReadError when run_info.json
    is not a write2DPES run's, or the file is missing, lacks a number, or its points are not nR blocks of the same
    nZ increasing Z points at increasing r (a damaged or a foreign file).

    Source: <run>/pes_cut_2d.dat (docs/output_files.md, "Data files") and run_info.json ("Run metadata").
    """
    try:
        where, grids = info.data["pes_cut"], info.data["grids"]
        position = (float(where["x"]), float(where["y"]), float(where["theta_deg"]), float(where["phi_deg"]))
        n_z, n_r = int(grids["nz"]), int(grids["nr"])
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{info.directory}: run_info.json without a valid {error} of a write2DPES run") from None
    path = run / FILE
    r_bohr, z_bohr, v_hartree = blocks_of(numbers(path, text_of(path), 4), n_z, n_r, path)
    return PesCut(v_hartree, z_bohr, r_bohr, *position)


def blocks_of(columns: np.ndarray, n_z: int, n_r: int, path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The r grid, the Z grid and V(Z, r) (hartree) of the file's rows (r, Z, V hartree, V eV), which run over Z
    fastest in one block per r: n_r blocks of the same n_z Z points, both grids increasing; ReadError otherwise."""
    if len(columns) != n_z * n_r:
        raise ReadError(f"{path}: {len(columns)} points, not nR x nZ = {n_r} x {n_z} (run_info.json); run it again")
    r_blocks, z_blocks = columns[:, 0].reshape(n_r, n_z), columns[:, 1].reshape(n_r, n_z)
    r_bohr, z_bohr = r_blocks[:, 0].copy(), z_blocks[0].copy()
    same = np.all(r_blocks == r_bohr[:, None]) and np.all(z_blocks == z_bohr)
    if not (same and np.all(np.diff(r_bohr) > 0) and np.all(np.diff(z_bohr) > 0)):
        raise ReadError(
            f"{path}: the points are not one block of the same increasing Z points per increasing r; "
            "give a pes_cut_2d.dat written by spodvr"
        )
    return r_bohr, z_bohr, columns[:, 2].reshape(n_r, n_z).T.copy()
