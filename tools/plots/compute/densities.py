"""The wave packet over (Z, r) from densities/ (writeDensities): the probability of both grids together, its
profiles in Z and r, the mean positions, and the norm per region of the grids against time.

The scattering grid (nz_scattering points) and the specular grid (nz_specular >= nz_scattering points) start at
the same Z with the same dZ, so a cell iz < nz_scattering holds the scattering wave's and the specular wave's
probability; beyond it only the specular wave's (docs/output_files.md, "densities/"). Region edges are the
program's own grid indices from run_info.json (nz_potential, iz_ana, the absorbers' points), so a point on an
edge falls on the same side as in the program.
"""

from dataclasses import dataclass

import numpy as np

from compute.grids import density_per_volume
from readers.arrays import Series
from readers.errors import ReadError
from readers.run_info import RunInfo

DATASET_OF_GRID = {"scattering": "prob_scat_zr", "specular": "prob_spec_zr"}  # an absorber's grid -> its cells
LOW_Z_ABSORBERS = {"ZOptLow", "ZSpOptLow"}  # absorbing from the grid's start up to their end


@dataclass
class DensityGrids:
    """The (Z, r) grids of densities/ (bohr)."""

    z_bohr: np.ndarray  # (nz_specular,): the specular grid, whose first nz_scattering points are the scattering grid's
    r_bohr: np.ndarray  # (nr,)
    nz_scattering: int
    dz_bohr: float
    dr_bohr: float


def density_grids(series: Series) -> DensityGrids:
    """The grids recorded in densities/manifest.json."""
    g = series.manifest["grids"]
    z = g["z_start"] + g["dz"] * np.arange(g["nz_specular"])
    r = g["r_start"] + g["dr"] * np.arange(g["nr"])
    return DensityGrids(z, r, int(g["nz"]), float(g["dz"]), float(g["dr"]))


def usable_records(series: Series) -> tuple[Series, tuple[str, ...]]:
    """The records before the first one with a NaN or inf, and a note when some were left out. ReadError when no
    record is left (the run stopped before its first averages time, or diverged at once)."""
    n = len(series.times_au) if series.first_non_finite is None else series.first_non_finite
    if n == 0:
        raise ReadError(f"{series.directory}: no finite record; the run stopped before its first record or diverged")
    if n == len(series.times_au):
        return series, ()
    arrays = {name: values[:n] for name, values in series.arrays.items()}
    cut = Series(series.directory, series.manifest, series.times_au[:n], arrays, None)
    return cut, (f"records from number {n} on hold a NaN or inf (the propagation diverged): not drawn",)


def combined_probability(series: Series) -> np.ndarray:
    """(n_records, nz_specular, nr): the probability per cell of both waves together (the scattering wave's
    cells padded with zero beyond its grid)."""
    total = np.array(series.arrays["prob_spec_zr"], dtype=float)
    nz_scattering = series.arrays["prob_scat_zr"].shape[1]
    total[:, :nz_scattering, :] += series.arrays["prob_scat_zr"]
    return total


def mean_positions_bohr(total: np.ndarray, grids: DensityGrids) -> tuple[np.ndarray, np.ndarray]:
    """<Z> and <r> per record (bohr) of the combined probability `total`, each cell at its grid point."""
    norm = np.sum(total, axis=(1, 2))
    z_mean = np.sum(total * grids.z_bohr[None, :, None], axis=(1, 2)) / norm
    r_mean = np.sum(total * grids.r_bohr[None, None, :], axis=(1, 2)) / norm
    return z_mean, r_mean


def r_absorber_first_index(info: RunInfo, nr: int) -> int:
    """The first r point of the r absorber (ROpt, its points in run_info.json): r_c, where the outer r region
    starts; nr (no outer region) without one."""
    for absorber in info.data["absorbers"]["list"]:
        if absorber["coordinate"] == "r" and absorber["points"] is not None:
            return int(absorber["points"]["first"])
    return nr


def band_edges(info: RunInfo, grids: DensityGrids) -> list[tuple[str, int, int]]:
    """The Z bands for r < r_c: (label, first index, end index (exclusive)) on the specular grid's indices.
    The edges: nz_potential (the first point at or beyond zasymp), iz_ana and the scattering grid's end."""
    z_asymptotic = info.data["z_asymptotic"]
    iz_asymptotic = int(info.data["grids"]["nz_potential"])
    iz_analysis = int(info.data["analysis"]["iz_ana"])
    return [
        (f"r < r_c, Z < zasymp ({z_asymptotic:g} bohr)", 0, iz_asymptotic),
        ("r < r_c, zasymp <= Z < Zana", iz_asymptotic, iz_analysis),
        ("r < r_c, Zana <= Z < end of the scattering grid", iz_analysis, grids.nz_scattering),
        ("r < r_c, Z beyond the scattering grid", grids.nz_scattering, len(grids.z_bohr)),
    ]


def band_norms(total: np.ndarray, bands: list[tuple[str, int, int]], ir_c: int) -> dict[str, np.ndarray]:
    """The probability of `total` (n_records, nz_specular, nr) in each Z band for r points below ir_c."""
    return {label: np.sum(total[:, first:end, :ir_c], axis=(1, 2)) for label, first, end in bands}


def absorber_contents(series: Series, info: RunInfo, total: np.ndarray) -> dict[str, np.ndarray]:
    """The probability inside each absorber of run_info.json (its points) on its own grid; the r absorber acts on
    both grids, so its content is that of the combined probability `total`."""
    contents = {}
    for absorber in info.data["absorbers"]["list"]:
        points = absorber["points"]
        if points is None:
            continue
        inside = slice(int(points["first"]), int(points["last"]) + 1)
        if absorber["coordinate"] == "r":
            contents[absorber["name"]] = np.sum(total[:, :, inside], axis=(1, 2))
        elif absorber["grid"] in DATASET_OF_GRID:
            cells = series.arrays[DATASET_OF_GRID[absorber["grid"]]]
            contents[absorber["name"]] = np.sum(cells[:, inside, :], axis=(1, 2))
        else:
            raise ReadError(f"run_info.json: absorber {absorber['name']} on grid {absorber['grid']!r}; update the tool")
    return contents


@dataclass
class RegionNorms:
    """The norm per region against time (n_records each): `regions` partition the record's norm, `absorbed` is
    1 - norm, and the absorbers' contents overlap the regions."""

    regions: dict[str, np.ndarray]  # label -> norm in that region
    absorbed: np.ndarray  # 1 - norm: what the absorbers have removed since t = 0
    absorbers: dict[str, np.ndarray]  # absorber name -> the probability inside it on its grid(s)


def region_norms(series: Series, info: RunInfo, grids: DensityGrids) -> RegionNorms:
    """The norm in r >= r_c (both grids; r_c is where the r absorber starts) and, for r < r_c, in the Z bands of
    band_edges; together they make the record's norm. Plus 1 - norm and every absorber's content."""
    total = combined_probability(series)
    ir_c = r_absorber_first_index(info, len(grids.r_bohr))
    regions = {}
    if ir_c < len(grids.r_bohr):
        regions[f"r >= r_c = ROptStart ({grids.r_bohr[ir_c]:g} bohr)"] = np.sum(total[:, :, ir_c:], axis=(1, 2))
    regions.update(band_norms(total, band_edges(info, grids), ir_c))
    absorbed = 1.0 - np.asarray(series.arrays["norm"])
    return RegionNorms(regions, absorbed, absorber_contents(series, info, total))


def z_profile_per_bohr(total: np.ndarray, grids: DensityGrids) -> np.ndarray:
    """(n_records, nz_specular): the probability per bohr of Z, summed over r."""
    return density_per_volume(np.sum(total, axis=2), grids.dz_bohr)


def r_profile_per_bohr(total: np.ndarray, grids: DensityGrids) -> np.ndarray:
    """(n_records, nr): the probability per bohr of r, summed over Z."""
    return density_per_volume(np.sum(total, axis=1), grids.dr_bohr)


def specular_line_bohr(series: Series, grids: DensityGrids) -> np.ndarray:
    """The Z of the specular analysis surface per record (bohr)."""
    return grids.z_bohr[0] + grids.dz_bohr * series.arrays["surface_index"].astype(float)


@dataclass
class Mark:
    """A line a density map marks: its label, position (bohr) and, for an absorber, its grid ("scattering",
    "specular" or "both"; None for the other lines)."""

    label: str
    position_bohr: float
    grid: str | None


@dataclass
class Marks:
    """The lines a density map marks along Z and along r."""

    z_lines: list[Mark]
    r_lines: list[Mark]


def density_marks(info: RunInfo, grids: DensityGrids) -> Marks:
    """zasymp, Zana, the end of the scattering grid and where each absorber begins to absorb (its start; the
    low-Z absorbers, which absorb from the grid's start, at their end)."""
    z_lines = [
        Mark("zasymp", float(info.data["z_asymptotic"]), None),
        Mark("Zana", float(info.data["analysis"]["z_ana"]), None),
        Mark("end of the scattering grid", float(grids.z_bohr[grids.nz_scattering - 1]), None),
    ]
    r_lines = []
    for absorber in info.data["absorbers"]["list"]:
        lines = r_lines if absorber["coordinate"] == "r" else z_lines
        if absorber["name"] in LOW_Z_ABSORBERS:
            lines.append(Mark(f"{absorber['name']} ends", float(absorber["end"]), absorber["grid"]))
        else:
            lines.append(Mark(f"{absorber['name']} starts", float(absorber["start"]), absorber["grid"]))
    return Marks(z_lines, r_lines)


def movie_frames(n_records: int, max_frames: int) -> np.ndarray:
    """At most max_frames (>= 2) record indices, evenly spaced, the first and the last included."""
    if max_frames < 2:
        raise ValueError(f"max_frames {max_frames}: a movie needs at least 2 frames")
    if n_records <= max_frames:
        return np.arange(n_records)
    return np.unique(np.round(np.linspace(0, n_records - 1, max_frames)).astype(int))
