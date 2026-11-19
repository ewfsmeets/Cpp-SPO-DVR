"""potential/ (writePotential): the run's potential on its grid at set-up, an array directory read by
readers/arrays.py (docs/output_files.md, "potential/"): the minimum over the DVR orientations and its orientation
at every grid point of Z < zasymp, and V at every orientation at the sites of PlotSites. This module knows the
directory's name, its datasets, their axes and shapes, and its manifest fields.

Shadows: src/report/potential_files.cpp (potential_datasets, potential_manifest_fields)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.arrays import read_arrays
from readers.errors import ReadError
from readers.sites import Site, sites_of

DIRECTORY = "potential"
AXES = {
    "v_min": ["Z", "r", "x", "y"],
    "v_min_node": ["Z", "r", "x", "y", "node"],
    "v_sites": ["site", "Z", "r", "theta", "phi"],
}


@dataclass
class PotentialLimits:
    """How the run changed the PES (the manifest's limits): switching to the gas phase and the VMin/VMax clamps."""

    use_switching: bool
    z_min_switch_bohr: float
    z_max_switch_bohr: float
    impose_v_min: bool
    v_min_clamp_hartree: float
    impose_v_max: bool
    v_max_clamp_hartree: float


@dataclass
class PotentialData:
    """potential/ of a run (hartree), with its grids (bohr) and DVR nodes (rad)."""

    v_min_hartree: np.ndarray  # (nz_potential, nr, nx, ny): the minimum over the orientations
    v_min_node: np.ndarray  # (nz_potential, nr, nx, ny, 2): (itheta, iphi) of that minimum
    v_sites_hartree: np.ndarray | None  # (n_sites, nz_potential, nr, ntheta, nphi), None without sites
    z_bohr: np.ndarray  # (nz_potential,): the scattering grid's points of Z < zasymp
    r_bohr: np.ndarray  # (nr,)
    nx: int
    ny: int
    theta_rad: np.ndarray  # (ntheta,): the DVR theta nodes (theta <= pi/2 only for a homonuclear molecule)
    phi_rad: np.ndarray  # (nphi,)
    homonuclear: bool
    sites: list[Site]
    limits: PotentialLimits


def read_potential(run: Path) -> PotentialData:
    """potential/ of a run directory; ReadError when a dataset or field is missing or does not have the axes and
    shape of the grids.

    Source: <run>/potential/ (docs/output_files.md, "potential/").
    """
    directory = run / DIRECTORY
    data = read_arrays(directory)
    m = data.manifest
    source = directory / "manifest.json"
    try:
        g = m["grids"]
        z = g["z_start"] + g["dz"] * np.arange(g["nz_potential"])
        r = g["r_start"] + g["dr"] * np.arange(g["nr"])
        nx, ny = int(g["nx"]), int(g["ny"])
        theta, phi = np.asarray(m["theta_nodes_rad"], dtype=float), np.asarray(m["phi_nodes_rad"], dtype=float)
        limits = limits_of(m["limits"])
        homonuclear = bool(m["homonuclear"])
        sites = sites_of(m["sites"], ny, source)
        shapes = {
            "v_min": (len(z), len(r), nx, ny),
            "v_min_node": (len(z), len(r), nx, ny, 2),
            "v_sites": (len(sites), len(z), len(r), len(theta), len(phi)),
        }
        check_datasets(source, m["datasets"], shapes, has_sites=bool(sites))
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{source}: no valid {error}; give potential/ written by spodvr") from None
    arrays = data.arrays
    v_sites = arrays.get("v_sites")
    return PotentialData(
        arrays["v_min"], arrays["v_min_node"], v_sites, z, r, nx, ny, theta, phi, homonuclear, sites, limits
    )


def limits_of(entry: dict) -> PotentialLimits:
    """The manifest's limits as PotentialLimits."""
    return PotentialLimits(
        bool(entry["use_switching"]),
        float(entry["z_min_switch"]),
        float(entry["z_max_switch"]),
        bool(entry["impose_v_min"]),
        float(entry["v_min_hartree"]),
        bool(entry["impose_v_max"]),
        float(entry["v_max_hartree"]),
    )


def check_datasets(source: Path, datasets: dict, shapes: dict[str, tuple[int, ...]], has_sites: bool) -> None:
    """ReadError unless v_min and v_min_node are there, v_sites exactly when there are sites, and each has the axes
    and the shape of the grids this reader indexes."""
    needed = {"v_min", "v_min_node"} | ({"v_sites"} if has_sites else set())
    if set(datasets) != needed:
        raise ReadError(
            f"{source}: datasets {sorted(datasets)}, not {sorted(needed)}; give potential/ written by spodvr"
        )
    for name, spec in datasets.items():
        if spec["axes"] != AXES[name] or tuple(spec["shape"]) != shapes[name]:
            raise ReadError(
                f"{source}: {name} has axes {spec['axes']} and shape {spec['shape']}, not {AXES[name]} and "
                f"{list(shapes[name])}; update the tool"
            )
