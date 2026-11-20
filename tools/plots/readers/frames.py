"""frames/ (dtFrames): reductions of the whole wave packet psi_full = psi_scat + s2s chi (chi, the specular wave,
added in the initial rotational state with the specular-to-scattering factor s2s = exp(i k_par.R) unit_cell_norm)
at tStart and every dtFrames, a series directory read by readers/arrays.py (docs/output_files.md, "frames/").
This module knows the directory's name, the datasets it must hold (compute uses them by their C++ names), its
grids, the rotational states, the sites and the optional datasets of FrameExtras with their grids.

Shadows: src/report/frame_files.cpp (frame_datasets, extra_datasets, extras_manifest, frame_manifest_fields)
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readers.arrays import Series
from readers.errors import ReadError
from readers.segments import read_joined_series
from readers.sites import Site, sites_of

DIRECTORY = "frames"
DATASETS = ("prob_zrxy", "prob_tail_zr", "pop_jmj_z", "surface_index")
GRID_FIELDS = ("z_start", "dz", "nz", "nz_specular", "r_start", "dr", "nr", "nx", "ny")
# The datasets of each FrameExtras extra (src/frames/frame_extras.hpp).
EXTRA_DATASETS = {
    "angular": ("prob_angle_z",),
    "diffraction": ("pop_nm_z",),
    "vibration": ("pop_vj_z",),
    "rotation_moments": ("p2_numerator_z", "rot_numerator_z"),
    "momentum": ("prob_kr", "prob_kz_spec"),
}


@dataclass(frozen=True, kw_only=True)
class AngularNodes:
    """The DVR grid of the angular extra: theta_t (Gauss-Legendre nodes; for a homonuclear molecule only those up to
    pi/2, with weights that count the mirror node pi - theta_t too) and the equidistant phi_k."""

    theta_rad: np.ndarray  # (ntheta,)
    theta_weights: np.ndarray  # (ntheta,): the weights w_t; a node's solid angle is w_t 2 pi / nphi
    phi_rad: np.ndarray  # (nphi,)
    homonuclear: bool


@dataclass(frozen=True, kw_only=True)
class FrameExtrasInfo:
    """Which FrameExtras a run wrote, with the grids their datasets need (None for an extra that is off)."""

    names: tuple[str, ...] = ()
    angular: AngularNodes | None = None
    vibration_j_values: np.ndarray | None = None  # (nJ,): the J of pop_vj_z's last axis, ascending


@dataclass(frozen=True, kw_only=True)
class Frames:
    """frames/ of a run (its segments joined): the records and what they are recorded on."""

    series: Series  # times_au and the datasets by name, record first (DATASETS, and EXTRA_DATASETS of the extras)
    z_bohr: np.ndarray  # (nz_specular,): the specular grid, whose first nz_scattering points are the scattering grid's
    r_bohr: np.ndarray  # (nr,)
    nz_scattering: int
    nx: int
    ny: int
    states: np.ndarray  # (nrot, 2) int: (J, mJ) of each rotational state, in the order of pop_jmj_z's last axis
    initial_index: int  # the state of (JIn, MJIn), which holds chi
    sites: list[Site]
    starts_au: list[float]  # the segments' start times
    extras: FrameExtrasInfo = FrameExtrasInfo()  # none in a run without FrameExtras


def read_frames(directories: list[Path]) -> Frames:
    """frames/ of a run's segments, joined (readers.segments.read_joined_series); ReadError when a dataset, a grid
    field, the rotational states, a site's grid point or an extra's dataset or grid is missing or wrong, or a
    dataset's shape does not follow from the grids, the states and the extras' grids."""
    series, starts = read_joined_series(directories, DIRECTORY)
    manifest = series.manifest
    source = series.directory / "manifest.json"
    try:
        g = manifest["grids"]
        missing = [d for d in DATASETS if d not in series.arrays] + [f for f in GRID_FIELDS if f not in g]
        if missing:
            raise ReadError(f"{source}: no {', '.join(missing)}; give frames/ written by spodvr with dtFrames")
        states = np.array(manifest["rotational_states"], dtype=int).reshape(-1, 2)
        sites = sites_of(manifest["sites"], int(g["ny"]), source)
        frames = Frames(
            series=series,
            z_bohr=g["z_start"] + g["dz"] * np.arange(g["nz_specular"]),
            r_bohr=g["r_start"] + g["dr"] * np.arange(g["nr"]),
            nz_scattering=int(g["nz"]),
            nx=int(g["nx"]),
            ny=int(g["ny"]),
            states=states,
            initial_index=int(manifest["initial_index"]),
            sites=sites,
            starts_au=starts,
            extras=extras_of(manifest),
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{source}: no valid {error}; give frames/ written by spodvr with dtFrames") from None
    check_extras(frames, source)
    check_shapes(frames, source)
    return frames


def extras_of(manifest: dict) -> FrameExtrasInfo:
    """The manifest's "extras" (written only with FrameExtras: without it, no extras)."""
    if "extras" not in manifest:
        return FrameExtrasInfo()
    entry = manifest["extras"]
    j_values = entry.get("vibration_j_values")
    return FrameExtrasInfo(
        names=tuple(entry["names"]),
        angular=angular_nodes_of(entry),
        vibration_j_values=None if j_values is None else np.array(j_values, dtype=int),
    )


def angular_nodes_of(entry: dict) -> AngularNodes | None:
    """The DVR grid of the manifest's "extras" (None without the angular extra)."""
    if "angular" not in entry:
        return None
    angular = entry["angular"]
    return AngularNodes(
        theta_rad=np.array(angular["theta_rad"], dtype=float),
        theta_weights=np.array(angular["theta_weights"], dtype=float),
        phi_rad=np.array(angular["phi_rad"], dtype=float),
        homonuclear=bool(angular["homonuclear"]),
    )


def check_extras(frames: Frames, source: Path) -> None:
    """ReadError when an extra is unknown, its datasets are missing, or the angular or vibration extra lacks its
    grid."""
    extras = frames.extras
    for name in extras.names:
        if name not in EXTRA_DATASETS:
            raise ReadError(f"{source}: unknown FrameExtras name {name!r}; update the tool")
        missing = [d for d in EXTRA_DATASETS[name] if d not in frames.series.arrays]
        if missing:
            raise ReadError(
                f"{source}: FrameExtras {name} without {', '.join(missing)}; give frames/ written by spodvr"
            )
    if ("angular" in extras.names) != (extras.angular is not None):
        raise ReadError(
            f"{source}: the angular extra without its DVR grid, or the grid without the extra; give frames/ "
            "written by spodvr"
        )
    if ("vibration" in extras.names) != (extras.vibration_j_values is not None):
        raise ReadError(
            f"{source}: the vibration extra without its J values, or the values without the extra; give "
            "frames/ written by spodvr"
        )


def expected_shapes(frames: Frames) -> dict[str, tuple[int, ...]]:
    """The record shape of every dataset the run wrote, from the grids, the rotational states and the extras'
    grids (after check_extras, an extra's grid is there exactly when its name is)."""
    nz, nr, n_rot = frames.nz_scattering, len(frames.r_bohr), len(frames.states)
    expected = {
        "prob_zrxy": (nz, nr, frames.nx, frames.ny),
        "prob_tail_zr": (len(frames.z_bohr) - nz, nr),
        "pop_jmj_z": (nz, n_rot),
        "surface_index": (),
    }
    extras = frames.extras
    if extras.angular is not None:
        expected["prob_angle_z"] = (nz, len(extras.angular.theta_rad), len(extras.angular.phi_rad))
    if "diffraction" in extras.names:
        expected["pop_nm_z"] = (nz, frames.nx, frames.ny)
    if extras.vibration_j_values is not None:
        expected["pop_vj_z"] = (nz, nr, len(extras.vibration_j_values))
    if "rotation_moments" in extras.names:
        expected["p2_numerator_z"] = expected["rot_numerator_z"] = (nz,)
    if "momentum" in extras.names:
        expected["prob_kr"] = (nr,)
        expected["prob_kz_spec"] = (len(frames.z_bohr),)
    return expected


def check_shapes(frames: Frames, source: Path) -> None:
    """ReadError when a dataset's record shape is not the one the grids, the rotational states and the extras'
    grids give."""
    for name, shape in expected_shapes(frames).items():
        if frames.series.arrays[name].shape[1:] != shape:
            got = frames.series.arrays[name].shape[1:]
            raise ReadError(f"{source}: {name} has records of shape {got}, not {shape}; give frames/ written by spodvr")
