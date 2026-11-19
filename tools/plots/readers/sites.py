"""The sites of PlotSites, named grid points of the surface cell, as the program records them in a plot-data
manifest (potential/manifest.json "sites"): the fractions along the lattice vectors and the grid point (ix, iy),
with ixy = iy + ny ix. (run_info.json plot_data lists the sites too, without ixy.)

Shadows: src/input/validate_outputs.cpp (plot_site_ixy)
"""

from dataclasses import dataclass
from pathlib import Path

from readers.errors import ReadError


@dataclass
class Site:
    """A site of PlotSites: its name, its fractions along the two lattice vectors and its grid point."""

    name: str
    fraction_x: float
    fraction_y: float
    ix: int
    iy: int


def sites_of(entries: list[dict], ny: int, source: Path) -> list[Site]:
    """The sites of a list of entries (name, x, y, ix, iy, ixy); ReadError when one lacks a field or its ixy is not
    iy + ny ix (a changed layout would mix up x and y silently)."""
    sites = []
    for entry in entries:
        try:
            site = Site(entry["name"], float(entry["x"]), float(entry["y"]), int(entry["ix"]), int(entry["iy"]))
            ixy = int(entry["ixy"])
        except (KeyError, TypeError, ValueError) as error:
            raise ReadError(f"{source}: a site without a valid {error}; give plot data written by spodvr") from None
        if ixy != site.iy + ny * site.ix:
            raise ReadError(f"{source}: site {site.name} has ixy {ixy}, not iy + ny ix; update the tool")
        sites.append(site)
    return sites
