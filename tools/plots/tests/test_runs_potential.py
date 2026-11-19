"""potential/ (writePotential) on the fixture runs, each check with its negative test: V_min at each site is the
minimum of the site's V over the orientations and its node holds it; the hard wall's V equals its formula, computed
here in numpy; the hard wall's minimax barrier is its gas-phase V at r_c just below zasymp."""

from pathlib import Path

import numpy as np
import pytest

from compute.densities import r_absorber_first_index
from compute.potential import barrier_on
from readers.potential import read_potential
from readers.run_info import read_run_info

RUNS = ["hardwall_tiny", "dailight_oblique", "potential_only"]
HARD_WALL_D_EV = 4.8286  # the hard wall's gas-phase well depth (src/potential/hard_wall_potential.cpp)
FORMULA_TOLERANCE = 1e-15  # relative to d + V_z: numpy's exp and pow against glibc's, and the eV -> hartree factor


def hard_wall_ev(z_bohr: np.ndarray, r_bohr: np.ndarray, switched: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """The hard wall's V(Z, r) in eV and its wall term V_z (Shadows: src/potential/hard_wall_potential.cpp
    HardWallPotential::value, gas_phase), on the (Z, r) grid; switched=False leaves out the wall's cosine switch
    (a porting slip)."""
    z, r = np.meshgrid(z_bohr, r_bohr, indexing="ij")
    a1, a2, a3, a4, re = 2.282, 1.555, 0.7533, 2.23, 1.40
    x = r - re
    v_gas = -HARD_WALL_D_EV * (1.0 + a1 * x + a2 * x**2.0 + a3 * x**3.0) * np.exp(-a4 * x) + HARD_WALL_D_EV
    v_z = 1.0 / (z * z * z * z + 0.001)
    z_max, z_min = 7.0, 0.1778
    switch = 1.0 - (0.5 - 0.5 * np.cos(np.pi * (z - z_min) / (z_max - z_min)))
    factor = np.where(z < z_min, 1.0, switch) if switched else 1.0
    return np.where(z < z_max, v_gas + factor * v_z, v_gas), v_z


@pytest.mark.parametrize("run", RUNS)
def test_sites_hold_the_minimum(fixture_runs: Path, run: str) -> None:
    """At each site's grid point, V_min is the minimum of the site's V over the DVR orientations (bitwise) and V at
    the stored node (itheta, iphi) equals it; negative tests: V_min of one point x (1 + 1e-12), or NaN, fails exactly
    there, and the node read as (iphi, itheta) (a swap) fails where the two orders give different V."""
    p = read_potential(fixture_runs / run)
    assert p.sites and p.v_sites_hartree is not None
    for k, site in enumerate(p.sites):
        v_min = np.array(p.v_min_hartree[:, :, site.ix, site.iy])
        sites = np.asarray(p.v_sites_hartree[k])
        lowest = sites.min(axis=(2, 3))
        assert np.count_nonzero(v_min != lowest) == 0
        node = np.asarray(p.v_min_node[:, :, site.ix, site.iy])
        iz, ir = np.indices(v_min.shape)
        assert np.count_nonzero(sites[iz, ir, node[..., 0], node[..., 1]] != v_min) == 0
        for planted in (v_min[3, 4] * (1 + 1e-12), np.nan):
            wrong = v_min.copy()
            wrong[3, 4] = planted
            assert np.count_nonzero(~(wrong == lowest)) == 1
        n_theta, n_phi = sites.shape[2:]
        swapped_theta, swapped_phi = node[..., 1] % n_theta, node[..., 0] % n_phi
        swapped = sites[iz, ir, swapped_theta, swapped_phi]
        differs = sites[iz, ir, node[..., 0], node[..., 1]] != swapped
        assert np.count_nonzero(swapped != v_min) == np.count_nonzero(differs)


def test_hard_wall_equals_its_formula(fixture_runs: Path) -> None:
    """hardwall_tiny (no switching to the gas phase, no VMin/VMax clamp): V_min and every site's V at every
    orientation equal the formula within 1e-15 of d + V_z (eV); negative test: the formula without the wall's cosine
    switch fails exactly at the points 0.1778 <= Z < 7 where the switch changes V by more than that."""
    run = fixture_runs / "hardwall_tiny"
    p = read_potential(run)
    assert not (p.limits.use_switching or p.limits.impose_v_min or p.limits.impose_v_max)
    expected, v_z = hard_wall_ev(p.z_bohr, p.r_bohr)
    tolerance = FORMULA_TOLERANCE * (HARD_WALL_D_EV + np.abs(v_z))
    got_min = np.asarray(p.v_min_hartree) * p.hartree_to_ev
    assert np.all(np.abs(got_min - expected[:, :, None, None]) <= tolerance[:, :, None, None])
    got_sites = np.asarray(p.v_sites_hartree) * p.hartree_to_ev
    assert np.all(np.abs(got_sites - expected[None, :, :, None, None]) <= tolerance[None, :, :, None, None])
    unswitched, _ = hard_wall_ev(p.z_bohr, p.r_bohr, switched=False)
    failing = ~(np.abs(got_min[:, :, 0, 0] - unswitched) <= tolerance)
    changed = np.abs(unswitched - expected) > tolerance
    assert changed.any() and np.array_equal(failing, changed)


def test_hard_wall_barrier_is_the_gas_phase_at_r_c(fixture_runs: Path) -> None:
    """The hard wall's wall term falls with Z, so the lowest way to r >= r_c runs along the last Z row: the barrier
    is V there at the r absorber's first point (the program's index, run_info.json), a value of the file exactly;
    negative test: a start at the first Z row (the wrong side of the grid) gives a higher barrier."""
    run = fixture_runs / "hardwall_tiny"
    p = read_potential(run)
    ir_c = r_absorber_first_index(read_run_info(run), len(p.r_bohr))
    v = np.asarray(p.v_min_hartree[:, :, 0, 0])
    barrier = barrier_on(v, ir_c, p.limits)
    assert barrier.height_hartree == v[-1, ir_c] and (barrier.iz, barrier.ir) == (v.shape[0] - 1, ir_c)
    wrong_side = barrier_on(v[::-1], ir_c, p.limits)
    assert wrong_side.height_hartree > barrier.height_hartree
