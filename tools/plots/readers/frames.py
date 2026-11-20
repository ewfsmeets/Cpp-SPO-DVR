"""frames/ (dtFrames): reductions of the whole wave packet psi_full = psi_scat + s2s chi (chi, the specular wave,
added in the initial rotational state with the specular-to-scattering factor s2s = exp(i k_par.R) unit_cell_norm)
at tStart and every dtFrames, a series directory read by readers/arrays.py (docs/output_files.md, "frames/").
This module knows the directory's name, the datasets it must hold (compute uses them by their C++ names), its
grids, the rotational states and the sites.

Shadows: src/report/frame_files.cpp (frame_datasets, frame_manifest_fields)
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


@dataclass(frozen=True, kw_only=True)
class Frames:
    """frames/ of a run (its segments joined): the records and what they are recorded on."""

    series: Series  # times_au and the datasets by name, record first (see DATASETS)
    z_bohr: np.ndarray  # (nz_specular,): the specular grid, whose first nz_scattering points are the scattering grid's
    r_bohr: np.ndarray  # (nr,)
    nz_scattering: int
    nx: int
    ny: int
    states: np.ndarray  # (nrot, 2) int: (J, mJ) of each rotational state, in the order of pop_jmj_z's last axis
    initial_index: int  # the state of (JIn, MJIn), which holds chi
    sites: list[Site]
    starts_au: list[float]  # the segments' start times


def read_frames(directories: list[Path]) -> Frames:
    """frames/ of a run's segments, joined (readers.segments.read_joined_series); ReadError when a dataset, a grid
    field, the rotational states or a site's grid point is missing or wrong, or a dataset's shape does not follow
    from the grids and the states."""
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
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ReadError(f"{source}: no valid {error}; give frames/ written by spodvr with dtFrames") from None
    check_shapes(frames, source)
    return frames


def check_shapes(frames: Frames, source: Path) -> None:
    """ReadError when a dataset's record shape is not the one the grids and the rotational states give."""
    nz, nr, n_rot = frames.nz_scattering, len(frames.r_bohr), len(frames.states)
    expected = {
        "prob_zrxy": (nz, nr, frames.nx, frames.ny),
        "prob_tail_zr": (len(frames.z_bohr) - nz, nr),
        "pop_jmj_z": (nz, n_rot),
        "surface_index": (),
    }
    for name, shape in expected.items():
        if frames.series.arrays[name].shape[1:] != shape:
            got = frames.series.arrays[name].shape[1:]
            raise ReadError(f"{source}: {name} has records of shape {got}, not {shape}; give frames/ written by spodvr")
