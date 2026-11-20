"""pes_cut_2d.dat (write2DPES) on the fixture runs: the cut at the bridge site and a DVR node equals V there in
dailight_oblique's potential/ (the same PES on the same grid, through the program's potential grid) wherever
neither the switch to the gas phase nor a VMin/VMax clamp changed it. Negative tests: one value moved just beyond
the tolerance fails at exactly that point; V at the neighbouring theta node fails at most compared points."""

from pathlib import Path

import numpy as np

from readers.pes_cut import PesCut, read_pes_cut
from readers.potential import PotentialData, read_potential
from readers.run_info import read_run_info

# Relative to max(|V|, 1 eV), as the C++ cross-check RegressionPesCutChecks.CutAtANodeEqualsThePotentialGrid
# (tests/test_regression_pes.cpp). The two differ by rounding only: theta passes through degrees in the cut
# (rad -> deg -> rad, about 1 ulp), and the potential grid converts eV to hartree with the inverse factor where the
# cut file divides (1 ulp).
TOLERANCE = 1e-13


def allowed_hartree(site_hartree: np.ndarray, to_ev: float) -> np.ndarray:
    """The tolerance at each point of the site's V (hartree)."""
    return TOLERANCE * np.maximum(np.abs(site_hartree) * to_ev, 1.0) / to_ev


def failing_points(values: np.ndarray, site: np.ndarray, compared: np.ndarray, to_ev: float) -> int:
    """The number of compared (Z, r) points where `values` and the site's V differ by more than the tolerance."""
    return int(np.count_nonzero(~(np.abs(values - site) <= allowed_hartree(site, to_ev)) & compared))


def compared_points(potential: PotentialData, site: np.ndarray) -> np.ndarray:
    """(nz_potential, nr) mask: the points where the run's V is the PES's: Z <= ZminSwitch and not a clamp."""
    limits = potential.limits
    unswitched = (potential.z_bohr <= limits.z_min_switch_bohr)[:, None] if limits.use_switching else True
    clamped = (limits.impose_v_max & (site >= limits.v_max_clamp_hartree)) | (
        limits.impose_v_min & (site <= limits.v_min_clamp_hartree)
    )
    return np.broadcast_to(unswitched, site.shape) & ~clamped


def node_of(potential: PotentialData, cut: PesCut) -> tuple[int, int]:
    """(itheta, iphi) of the DVR node the cut was taken at (its angles in degrees)."""
    itheta = int(np.argmin(np.abs(potential.theta_rad - np.radians(cut.theta_deg))))
    iphi = int(np.argmin(np.abs(potential.phi_rad - np.radians(cut.phi_deg))))
    assert abs(np.degrees(potential.theta_rad[itheta]) - cut.theta_deg) < 1e-12 and potential.phi_rad[iphi] == 0.0
    return itheta, iphi


def test_cut_equals_the_site_potential_at_its_node(fixture_runs: Path) -> None:
    """The cut fixture (dailight_oblique's input at the bridge, the last theta node, phi 0) against dailight_oblique's
    V at the bridge site and that node: none of more than 1000 compared points beyond the tolerance."""
    cut = read_pes_cut(fixture_runs / "pes_cut", read_run_info(fixture_runs / "pes_cut"))
    potential = read_potential(fixture_runs / "dailight_oblique")
    to_ev = read_run_info(fixture_runs / "dailight_oblique").hartree_to_ev
    (bridge,) = [k for k, site in enumerate(potential.sites) if site.name == "bridge"]
    itheta, iphi = node_of(potential, cut)
    nz = len(potential.z_bohr)
    assert np.array_equal(cut.z_bohr[:nz], potential.z_bohr) and np.array_equal(cut.r_bohr, potential.r_bohr)
    site = np.asarray(potential.v_sites_hartree[bridge, :, :, itheta, iphi])
    compared = compared_points(potential, site)
    values = cut.v_hartree[:nz]
    n_compared = np.count_nonzero(compared)
    assert n_compared > 1000
    assert failing_points(values, site, compared, to_ev) == 0

    # One compared value moved by twice its tolerance: exactly that point fails.
    points = np.argwhere(compared)
    iz, ir = points[len(points) // 2]
    moved = values.copy()
    moved[iz, ir] += 2.0 * allowed_hartree(site[iz, ir], to_ev)
    assert failing_points(moved, site, compared, to_ev) == 1
    # V at the neighbouring theta node (a cut at the wrong node): most compared points fail.
    neighbour = np.asarray(potential.v_sites_hartree[bridge, :, :, itheta - 1, iphi])
    assert failing_points(neighbour, site, compared, to_ev) > n_compared // 2
