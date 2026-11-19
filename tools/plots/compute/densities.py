"""The wave packet over (Z, r) from densities/ (writeDensities): the probability of both grids together, its
profiles in Z and r, the norm per region of the grids and the absorbers' contents against time.

The scattering grid (nz_scattering points) and the specular grid (nz_specular >= nz_scattering points) start at
the same Z with the same dZ, so a cell iz < nz_scattering holds the scattering wave's and the specular wave's
probability; beyond it only the specular wave's (docs/output_files.md, "densities/"). Region edges are the
program's own grid indices from run_info.json (nz_potential, iz_ana, the absorbers' points), so a point on an
edge falls on the same side as in the program.
"""

from dataclasses import dataclass

import numpy as np

from compute.grids import density_per_volume, r_absorber_first_index
from readers.arrays import Series
from readers.densities import DensityGrids
from readers.errors import ReadError
from readers.run_info import RunInfo

DATASET_OF_GRID = {"scattering": "prob_scat_zr", "specular": "prob_spec_zr"}  # a Z absorber's grid -> its cells
LOW_Z_ABSORBERS = {"ZOptLow", "ZSpOptLow"}  # absorbing on [start, end], strongest at their start (low Z)
WAVE_OF_GRID = {"scattering": "scattering wave", "specular": "specular wave", "both": "both waves"}


def combined_probability(series: Series, records: np.ndarray | None = None) -> np.ndarray:
    """(n_records, nz_specular, nr): the probability per cell of both waves together (the scattering wave's
    cells padded with zero beyond its grid), of the chosen records (all by default; chosen before anything is
    copied)."""
    chosen = slice(None) if records is None else records
    combined = np.array(series.arrays["prob_spec_zr"][chosen], dtype=float)
    scattering = series.arrays["prob_scat_zr"][chosen]
    combined[:, : scattering.shape[1], :] += scattering
    return combined


def band_edges(info: RunInfo, grids: DensityGrids) -> list[tuple[str, int, int]]:
    """The Z bands for r < r_c that hold grid points: (label, first index, end index (exclusive)) on the specular
    grid's indices. The edges: nz_potential (the first point at or beyond zasymp), iz_ana and the scattering
    grid's end (zasymp <= Z < Zana holds none when Zana = zasymp)."""
    z_asymptotic_bohr = info.data["z_asymptotic"]
    iz_asymptotic = int(info.data["grids"]["nz_potential"])
    iz_analysis = int(info.data["analysis"]["iz_ana"])
    bands = [
        (f"r < r_c, Z < zasymp ({z_asymptotic_bohr:g} bohr)", 0, iz_asymptotic),
        ("r < r_c, zasymp <= Z < Zana", iz_asymptotic, iz_analysis),
        ("r < r_c, Zana <= Z < end of the scattering grid", iz_analysis, grids.nz_scattering),
        ("r < r_c, Z beyond the scattering grid", grids.nz_scattering, len(grids.z_bohr)),
    ]
    return [(label, first, end) for label, first, end in bands if first < end]


def band_norms(combined: np.ndarray, bands: list[tuple[str, int, int]], ir_c: int) -> dict[str, np.ndarray]:
    """The probability of `combined` (n_records, nz_specular, nr) in each Z band for r points below ir_c."""
    return {label: np.sum(combined[:, first:end, :ir_c], axis=(1, 2)) for label, first, end in bands}


def absorber_grid(absorber: dict) -> str:
    """The grid an absorber of run_info.json acts on: "scattering", "specular" or "both" (the r absorber);
    ReadError for another."""
    if absorber["grid"] not in WAVE_OF_GRID:
        raise ReadError(f"run_info.json: absorber {absorber['name']} on grid {absorber['grid']!r}; update the tool")
    return absorber["grid"]


def absorber_contents(series: Series, info: RunInfo, combined: np.ndarray) -> dict[str, np.ndarray]:
    """The probability inside each absorber of run_info.json (its points) on its own grid; the r absorber acts on
    both grids, so its content is that of the combined probability `combined`."""
    contents = {}
    for absorber in info.data["absorbers"]["list"]:
        points = absorber["points"]
        if points is None:
            continue
        inside = slice(int(points["first"]), int(points["last"]) + 1)
        grid = absorber_grid(absorber)
        if grid == "both":
            contents[absorber["name"]] = np.sum(combined[:, :, inside], axis=(1, 2))
        else:
            cells = series.arrays[DATASET_OF_GRID[grid]]
            contents[absorber["name"]] = np.sum(cells[:, inside, :], axis=(1, 2))
    return contents


@dataclass
class RegionNorms:
    """The norm per region against time (n_records each): `regions` partition the record's norm, `removed` is
    1 - norm, and the absorbers' contents overlap the regions."""

    regions: dict[str, np.ndarray]  # label -> norm in that region
    removed: np.ndarray  # 1 - norm: removed since t = 0 by the absorbers and the rotational basis cut (JMax, MJMax)
    absorbers: dict[str, np.ndarray]  # absorber name -> the probability inside it on its grid(s)
    absorber_grids: dict[str, str]  # absorber name -> its grid: "scattering", "specular" or "both" (r)


def region_norms(series: Series, info: RunInfo, grids: DensityGrids, combined: np.ndarray) -> RegionNorms:
    """The norm in r >= r_c (both grids; r_c is where the r absorber starts) and, for r < r_c, in the Z bands of
    band_edges; together they make the record's norm. Plus 1 - norm and every absorber's content. `combined`:
    combined_probability(series)."""
    ir_c = r_absorber_first_index(info, len(grids.r_bohr))
    regions = {}
    if ir_c < len(grids.r_bohr):
        regions[f"r >= r_c = ROptStart ({grids.r_bohr[ir_c]:g} bohr)"] = np.sum(combined[:, :, ir_c:], axis=(1, 2))
    regions.update(band_norms(combined, band_edges(info, grids), ir_c))
    removed = 1.0 - np.asarray(series.arrays["norm"])
    contents = absorber_contents(series, info, combined)
    grids_of = {a["name"]: a["grid"] for a in info.data["absorbers"]["list"] if a["name"] in contents}
    return RegionNorms(regions, removed, contents, grids_of)


def z_profile_per_bohr(combined: np.ndarray, grids: DensityGrids) -> np.ndarray:
    """(n_records, nz_specular): the probability per bohr of Z of `combined`, summed over r."""
    return density_per_volume(np.sum(combined, axis=2), grids.dz_bohr)


def r_profile_per_bohr(combined: np.ndarray, grids: DensityGrids) -> np.ndarray:
    """(n_records, nr): the probability per bohr of r of `combined`, summed over Z."""
    return density_per_volume(np.sum(combined, axis=1), grids.dr_bohr)


def specular_inside_line(series: Series, info: RunInfo) -> np.ndarray:
    """The specular wave between Zana and its analysis line per record, not yet analysed: the rest norm of
    averages.dat (cells iz_ana < iz <= the surface index, every r)."""
    iz_analysis = int(info.data["analysis"]["iz_ana"])
    cells = series.arrays["prob_spec_zr"]
    indices = np.asarray(series.arrays["surface_index"])
    return np.array([np.sum(cells[k, iz_analysis + 1 : indices[k] + 1, :]) for k in range(len(indices))])


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
    """zasymp, Zana, the end of the scattering grid and where each absorber begins to absorb, naming the wave it
    acts on (its start; a low-Z absorber, which absorbs below its end, at its end)."""
    z_lines = [
        Mark("zasymp", float(info.data["z_asymptotic"]), None),
        Mark("Zana", float(info.data["analysis"]["z_ana"]), None),
        Mark("end of the scattering grid", float(grids.z_bohr[grids.nz_scattering - 1]), None),
    ]
    r_lines = []
    for absorber in info.data["absorbers"]["list"]:
        lines = r_lines if absorber["coordinate"] == "r" else z_lines
        wave = WAVE_OF_GRID[absorber_grid(absorber)]
        if absorber["name"] in LOW_Z_ABSORBERS:
            mark = Mark(f"{absorber['name']} (absorbs below; {wave})", float(absorber["end"]), absorber["grid"])
        else:
            mark = Mark(f"{absorber['name']} starts ({wave})", float(absorber["start"]), absorber["grid"])
        lines.append(mark)
    return Marks(z_lines, r_lines)
