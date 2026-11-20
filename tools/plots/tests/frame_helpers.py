"""Helpers of the frames/ tests (test_runs_frames, test_runs_frame_extras); not a test module:
- the whole wave packet psi_full = psi_scat + s2s chi of a fixture run, rebuilt in numpy from its checkpoint.bin with
  s2s from its Cartesian formula (the independent input of the numpy checks);
- the count of entries outside a tolerance (NaN counts);
- a copy of a fixture run with a changed frames/manifest.json (the readers' refusals).
Imported as tests.frame_helpers (tests/ has no __init__.py: a namespace package)."""

import json
from pathlib import Path

import numpy as np

from readers.binary import read_checkpoint
from readers.frames import Frames
from readers.run_info import RunInfo


def s2s_of(info: RunInfo, nx: int, ny: int) -> np.ndarray:
    """(nx, ny): the specular-to-scattering factor s2s = exp(i k_par . R) unit_cell_norm from its Cartesian formula:
    k_par of magnitude k_parallel at the azimuth phi_incidence from the first lattice vector, R = (x + y cos(gamma),
    y sin(gamma)) with x = ix dx along the first lattice vector and y = iy dy along the second (independent of the
    program's table, which uses skewed components). Only the oblique runs test its phase (k_par = 0 elsewhere)."""
    wave, grids, lattice = info.data["wave_packet"], info.data["grids"], info.data["lattice"]
    k = wave["k_parallel"] * np.array([np.cos(wave["phi_incidence_rad"]), np.sin(wave["phi_incidence_rad"])])
    gamma = lattice["skew_angle_rad"]
    x, y = np.meshgrid(np.arange(nx) * grids["dx"], np.arange(ny) * grids["dy"], indexing="ij")
    phase = k[0] * (x + y * np.cos(gamma)) + k[1] * y * np.sin(gamma)
    return np.exp(1j * phase) * grids["unit_cell_norm"]


def psi_full_of(run: Path, frames: Frames, s2s: np.ndarray, with_chi: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """psi_full (nz, nr, nx, ny, nrot) on the scattering grid with the given s2s (nx, ny), and chi (nz_specular, nr),
    from the run's checkpoint.bin (bohr^-2 and bohr^-1: |psi|^2 dZ dr dx dy sin(gamma) and |chi|^2 dZ dr are
    probabilities); with_chi false: psi_scat alone instead of psi_full."""
    checkpoint = read_checkpoint(run / "checkpoint.bin")
    nz, nr, nx, ny = frames.nz_scattering, len(frames.r_bohr), frames.nx, frames.ny
    psi = np.array(checkpoint.scattering).reshape(nz, nr, nx * ny, -1)  # izr = ir + nr iz, then ixy, irot
    chi = np.asarray(checkpoint.specular)
    if with_chi:
        psi[:, :, :, frames.initial_index] += s2s.reshape(-1)[None, None, :] * chi[:nz, :, None]
    return psi.reshape(nz, nr, nx, ny, -1), chi


def relative_misses(a: np.ndarray, b: np.ndarray, scale: np.ndarray | float, tolerance: float) -> int:
    """How many entries of a differ from b by more than tolerance x scale (NaN counts)."""
    return int(np.count_nonzero(~(np.abs(a - b) <= tolerance * scale)))


def run_with_frames_manifest(fixture: Path, work: Path, change: dict) -> Path:
    """A copy of a fixture run (its files linked) whose frames/manifest.json has the top-level fields in `change`
    replaced (None: removed), the data files linked."""
    work.mkdir()
    for path in fixture.iterdir():
        if path.name != "frames":
            (work / path.name).symlink_to(path.resolve())
    (work / "frames").mkdir()
    for path in (fixture / "frames").iterdir():
        if path.name != "manifest.json":
            (work / "frames" / path.name).symlink_to(path.resolve())
    manifest = json.loads((fixture / "frames" / "manifest.json").read_text(encoding="utf-8"))
    for key, value in change.items():
        if value is None:
            del manifest[key]
        else:
            manifest[key] = value
    (work / "frames" / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return work
