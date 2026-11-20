"""energy_resolved/ (EnergyResolvedSites, EnergyResolvedEnergies): the run's time sums of the wave packet at chosen
(x, y) grid points and collision energies, an array directory read by readers/arrays.py (docs/output_files.md,
"energy_resolved/"). This module knows the directory's name, its datasets, their axes and shapes, and its manifest
fields.

Shadows: src/report/energy_resolved_files.cpp (energy_resolved_datasets, manifest_fields)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.arrays import read_arrays
from readers.errors import ReadError

DIRECTORY = "energy_resolved"
SITE_AXES = ["site", "energy", "Z", "r", "rot"]
SPECULAR_AXES = ["energy", "Z_specular", "r"]
AXES = {
    "sums": SITE_AXES,
    "specular_sums": SPECULAR_AXES,
    "sums_partial": SITE_AXES,
    "specular_sums_partial": SPECULAR_AXES,
    "partial_time": [],
    "last_time": [],
    "steps_summed": [],
}


@dataclass(frozen=True, kw_only=True)
class SumSite:
    """A site of the sums: its name, (x, y) grid point and the factor with which the specular wave enters the
    scattering grid there (s2s)."""

    name: str
    ix: int
    iy: int
    ixy: int  # iy + ny ix
    specular_to_scattering: complex


@dataclass(frozen=True, kw_only=True)
class EnergyResolved:
    """energy_resolved/ of a run: the sums (bohr^-2 times a.u. of time) and their grids (bohr)."""

    sites: list[SumSite]
    collision_ev: np.ndarray  # (n_energies,)
    total_hartree: np.ndarray  # (n_energies,): E_coll + E(vin, jin)
    sums: np.ndarray  # (n_sites, n_energies, nz, nr, n_rot) complex: the scattering wave's sums
    specular_sums: np.ndarray  # (n_energies, nz_specular, nr) complex: the specular wave's sums
    sums_partial: np.ndarray  # the same up to partial_time_au
    specular_sums_partial: np.ndarray
    partial_time_au: float
    last_time_au: float
    steps_summed: int
    z_bohr: np.ndarray  # (nz,): the scattering grid
    z_specular_bohr: np.ndarray  # (nz_specular,): the specular grid (the same start and spacing)
    r_bohr: np.ndarray  # (nr,)
    analysis_time_step_au: float


def read_energy_resolved(run: Path) -> EnergyResolved:
    """energy_resolved/ of a run directory; ReadError when a dataset or field is missing or does not have the axes
    and shape of the grids.

    Source: <run>/energy_resolved/ (docs/output_files.md, "energy_resolved/").
    """
    directory = run / DIRECTORY
    data = read_arrays(directory)
    m = data.manifest
    source = directory / "manifest.json"
    try:
        g = m["grids"]
        z = g["z_start"] + g["dz"] * np.arange(g["nz"])
        z_specular = g["z_start"] + g["dz"] * np.arange(g["nz_specular"])
        r = g["r_start"] + g["dr"] * np.arange(g["nr"])
        sites = [sum_site(entry) for entry in m["sites"]]
        collision = np.asarray(m["collision_ev"], dtype=float)
        site_shape = (len(sites), len(collision), len(z), len(r), int(g["n_rot"]))
        specular_shape = (len(collision), len(z_specular), len(r))
        check_datasets(source, m["datasets"], site_shape, specular_shape)
        total = np.asarray(m["total_hartree"], dtype=float)
        time_step = float(m["analysis_time_step"])
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{source}: no valid {error}; give energy_resolved/ written by spodvr") from None
    a = data.arrays
    return EnergyResolved(
        sites=sites,
        collision_ev=collision,
        total_hartree=total,
        sums=a["sums"],
        specular_sums=a["specular_sums"],
        sums_partial=a["sums_partial"],
        specular_sums_partial=a["specular_sums_partial"],
        partial_time_au=float(a["partial_time"]),
        last_time_au=float(a["last_time"]),
        steps_summed=int(a["steps_summed"]),
        z_bohr=z,
        z_specular_bohr=z_specular,
        r_bohr=r,
        analysis_time_step_au=time_step,
    )


def sum_site(entry: dict) -> SumSite:
    """A manifest's site entry as a SumSite."""
    real, imaginary = entry["specular_to_scattering"]
    return SumSite(
        name=str(entry["name"]),
        ix=int(entry["ix"]),
        iy=int(entry["iy"]),
        ixy=int(entry["ixy"]),
        specular_to_scattering=complex(real, imaginary),
    )


def check_datasets(source: Path, datasets: dict, site_shape: tuple[int, ...], specular_shape: tuple[int, ...]) -> None:
    """ReadError unless every dataset is there with the axes and the shape of the grids this reader indexes."""
    if set(datasets) != set(AXES):
        raise ReadError(f"{source}: datasets {sorted(datasets)}, not {sorted(AXES)}; give energy_resolved/ of spodvr")
    shapes = {
        "sums": site_shape,
        "sums_partial": site_shape,
        "specular_sums": specular_shape,
        "specular_sums_partial": specular_shape,
    }
    for name, spec in datasets.items():
        if spec["axes"] != AXES[name] or tuple(spec["shape"]) != shapes.get(name, ()):
            raise ReadError(
                f"{source}: {name} has axes {spec['axes']} and shape {spec['shape']}, not {AXES[name]} and "
                f"{list(shapes.get(name, ()))}; update the tool"
            )
