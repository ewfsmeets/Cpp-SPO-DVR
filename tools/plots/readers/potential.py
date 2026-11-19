"""potential/ (writePotential): the run's potential on its grid at set-up, an array directory read by
readers/arrays.py (docs/output_files.md, "potential/"): the minimum over the DVR orientations and its orientation
at every grid point of Z < zasymp, and V at every orientation at the sites of PlotSites. This module knows the
directory's name, its datasets and its manifest fields.

Shadows: src/report/potential_files.cpp (potential_datasets, potential_manifest_fields)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.arrays import read_arrays
from readers.errors import ReadError

DIRECTORY = "potential"


@dataclass
class Site:
    """A site of PlotSites: its name, fractional coordinates along the lattice vectors and grid point."""

    name: str
    x: float
    y: float
    ix: int
    iy: int


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
    sites: list[Site]
    limits: PotentialLimits
    hartree_to_ev: float


def read_potential(run: Path) -> PotentialData:
    """potential/ of a run directory; ReadError when a dataset or field is missing.

    Source: <run>/potential/ (docs/output_files.md, "potential/").
    """
    directory = run / DIRECTORY
    data = read_arrays(directory)
    m = data.manifest
    try:
        g = m["grids"]
        z = g["z_start"] + g["dz"] * np.arange(g["nz_potential"])
        r = g["r_start"] + g["dr"] * np.arange(g["nr"])
        sites = [Site(s["name"], s["x"], s["y"], s["ix"], s["iy"]) for s in m["sites"]]
        check_layout(directory, m, sites, int(g["ny"]))
        lim = m["limits"]
        limits = PotentialLimits(
            lim["use_switching"],
            lim["z_min_switch"],
            lim["z_max_switch"],
            lim["impose_v_min"],
            lim["v_min_hartree"],
            lim["impose_v_max"],
            lim["v_max_hartree"],
        )
        return PotentialData(
            data.arrays["v_min"],
            data.arrays["v_min_node"],
            data.arrays.get("v_sites"),
            z,
            r,
            int(g["nx"]),
            int(g["ny"]),
            np.asarray(m["theta_nodes_rad"]),
            np.asarray(m["phi_nodes_rad"]),
            sites,
            limits,
            float(m["hartree_to_ev"]),
        )
    except KeyError as error:
        raise ReadError(f"{directory / 'manifest.json'}: no {error}; give potential/ written by spodvr") from None


def check_layout(directory: Path, manifest: dict, sites: list[Site], ny: int) -> None:
    """ReadError unless the datasets have the axes this reader indexes ((Z, r, x, y) and sites at
    ixy = iy + ny ix): a changed layout would otherwise mix up x and y silently."""
    axes = {
        "v_min": ["Z", "r", "x", "y"],
        "v_min_node": ["Z", "r", "x", "y", "node"],
        "v_sites": ["site", "Z", "r", "theta", "phi"],
    }
    for name, expected in axes.items():
        spec = manifest["datasets"].get(name)
        if spec is not None and spec["axes"] != expected:
            raise ReadError(
                f"{directory / 'manifest.json'}: {name} has axes {spec['axes']}, not {expected}; update the tool"
            )
    for site, entry in zip(sites, manifest["sites"], strict=True):
        if entry["ixy"] != site.iy + ny * site.ix:
            raise ReadError(
                f"{directory / 'manifest.json'}: site {site.name} has ixy {entry['ixy']}, not iy + ny ix; update the tool"
            )
