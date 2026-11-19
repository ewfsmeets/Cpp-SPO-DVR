"""densities/ (writeDensities): the probability per (Z, r) cell of the scattering and the specular wave, the
specular analysis surface and the norm at tStart and every dtAve, a series directory read by readers/arrays.py
(docs/output_files.md, "densities/"). This module knows the directory's name, the datasets it must hold (compute
uses them by their C++ names) and its grids.

Shadows: src/report/densities_files.cpp (density_datasets, density_grids)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.arrays import Series
from readers.errors import ReadError
from readers.segments import read_joined_series

DIRECTORY = "densities"
DATASETS = ("prob_scat_zr", "prob_spec_zr", "surface_index", "norm")
GRID_FIELDS = ("z_start", "dz", "nz", "nz_specular", "r_start", "dr", "nr")


@dataclass
class DensityGrids:
    """The (Z, r) grids of densities/ (bohr)."""

    z_bohr: np.ndarray  # (nz_specular,): the specular grid, whose first nz_scattering points are the scattering grid's
    r_bohr: np.ndarray  # (nr,)
    nz_scattering: int
    dz_bohr: float
    dr_bohr: float


def density_grids(manifest: dict) -> DensityGrids:
    """The grids recorded in densities/manifest.json (both Z grids start at z_start with the same dZ)."""
    g = manifest["grids"]
    z = g["z_start"] + g["dz"] * np.arange(g["nz_specular"])
    r = g["r_start"] + g["dr"] * np.arange(g["nr"])
    return DensityGrids(z, r, int(g["nz"]), float(g["dz"]), float(g["dr"]))


def read_densities(directories: list[Path]) -> tuple[Series, DensityGrids, list[float]]:
    """densities/ of a run's segments, joined (readers.segments.read_joined_series), its grids and the segments'
    start times. ReadError when a dataset or a grid field is missing."""
    series, starts = read_joined_series(directories, DIRECTORY)
    grids = series.manifest.get("grids", {})
    missing = [d for d in DATASETS if d not in series.arrays] + [f for f in GRID_FIELDS if f not in grids]
    if missing:
        raise ReadError(
            f"{series.directory}/manifest.json: no {', '.join(missing)}; give densities/ written by spodvr with "
            "writeDensities"
        )
    return series, density_grids(series.manifest), starts
