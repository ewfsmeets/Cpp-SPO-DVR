"""Diffraction on the fixture runs: the Cartesian parallel momenta against the program's channel energies, the
incident direction, the plane of incidence at oblique incidence, and the specular angle."""

from pathlib import Path

import numpy as np
import pytest

from compute.analysis import analysis_energies
from compute.diffraction import diffraction, out_of_plane_fraction, scattering_angles
from compute.grids import parallel_momentum_of
from readers.binary import read_bk_coefficients
from readers.run_info import read_run_info
from readers.tables import read_gas_phase_states, read_s_matrix


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique"])
def test_parallel_energy_is_the_channel_energy_above_the_level(fixture_runs: Path, run: str) -> None:
    """For every analysed channel (v, state): |K|^2 / 2M = channel energy (BK file) - E(v, j)
    (gas_phase_states.dat), to 1e-12 hartree: K is the momentum the program's channel energies use."""
    info = read_run_info(fixture_runs / run)
    bk = read_bk_coefficients(fixture_runs / run / "bk_coefficients.bin")
    states = read_gas_phase_states(fixture_runs / run / "gas_phase_states.dat")
    level = {(int(v), int(j)): e for v, j, e in zip(states.v, states.j, states.energy_hartree, strict=True)}
    j, _, n, m, _ = bk.labels.T
    kx, ky = parallel_momentum_of(info, n, m)
    parallel = (kx * kx + ky * ky) / (2.0 * info.data["molecule"]["total_mass"])
    for v in range(bk.channel_energy_hartree.shape[0]):
        internal = np.array([level[(v, int(final_j))] for final_j in j])
        assert np.max(np.abs(bk.channel_energy_hartree[v] - internal - parallel)) <= 1e-12


def test_incident_parallel_momentum_points_along_phi(fixture_runs: Path) -> None:
    """K(0, 0) = k_par (cos phi, sin phi) with k_par and PhiParallel of run_info.json (17 degrees), to 1e-12
    relative; negative: PhiParallel used as radians misses by far more."""
    info = read_run_info(fixture_runs / "dailight_oblique")
    kx, ky = parallel_momentum_of(info, np.array([0]), np.array([0]))
    k_par, phi_deg = info.data["wave_packet"]["k_parallel"], info.data["wave_packet"]["phi_parallel_deg"]
    for phi_rad, close in ((np.radians(phi_deg), True), (phi_deg, False)):
        error = np.hypot(kx[0] - k_par * np.cos(phi_rad), ky[0] - k_par * np.sin(phi_rad)) / k_par
        assert (error <= 1e-12) == close, (phi_rad, error)


def test_only_the_specular_channel_is_in_the_plane_at_17_degrees(fixture_runs: Path) -> None:
    """At an azimuth off every mirror line of the 60 degree cell, no other channel lies in the plane of
    incidence, so the out-of-plane fraction is 1 - P(0, 0) / P_scattered; an --azimuth is ignored there."""
    run = fixture_runs / "dailight_oblique"
    info = read_run_info(run)
    result = diffraction(read_s_matrix(run / "s_matrix.tsv"), info, analysis_energies(info).collision_ev, None)
    assert [c for c, inside in zip(result.channels, result.in_plane, strict=True) if inside] == [(0, 0)]
    assert result.plane == "the plane of incidence"
    p_specular = result.p_channel[:, result.channels.index((0, 0))]
    expected = 1 - p_specular / np.sum(result.p_channel, axis=1)
    assert np.max(np.abs(out_of_plane_fraction(result) - expected)) <= 1e-15
    turned = diffraction(read_s_matrix(run / "s_matrix.tsv"), info, analysis_energies(info).collision_ev, 0.5)
    assert turned.plane == "the plane of incidence"


def test_specular_channel_leaves_at_the_incidence_angle(fixture_runs: Path) -> None:
    """The elastic, specular final state (v' = 0, j' = 0, (0, 0)) leaves at theta_f = atan2(k_par, k0), the
    angle of incidence of that energy, and along phi = PhiParallel, to 1e-12 rad (rounding of the arctangents);
    the probabilities are s_matrix.tsv's summed over mJ'."""
    run = fixture_runs / "dailight_oblique"
    info = read_run_info(run)
    table, bk = read_s_matrix(run / "s_matrix.tsv"), read_bk_coefficients(run / "bk_coefficients.bin")
    energies = analysis_energies(info)
    for energy in (0, len(energies.collision_ev) // 2, len(energies.collision_ev) - 1):
        angles = scattering_angles(table, bk, info, energy)
        elastic = (angles.v == 0) & (angles.j == 0) & (angles.n == 0) & (angles.m == 0)
        expected = np.arctan2(info.data["wave_packet"]["k_parallel"], energies.k0_per_bohr[energy])
        assert abs(angles.theta_rad[elastic][0] - expected) <= 1e-12
        assert abs(angles.phi_rad[elastic][0] - np.radians(info.data["wave_packet"]["phi_parallel_deg"])) <= 1e-12
        at_energy = (
            table["E_coll"] == table["E_coll"][np.argmin(np.abs(table["E_coll"] - energies.collision_ev[energy]))]
        )
        assert abs(np.sum(angles.probability) - np.sum(table["P"][at_energy])) <= 1e-13


def in_plane_channels(fixture_runs: Path, azimuth_deg: float | None) -> list[tuple[int, int]]:
    """The in-plane channels of the hard-wall run (normal incidence, a rectangular 2 x 2 cell) at an azimuth."""
    run = fixture_runs / "hardwall_tiny"
    info = read_run_info(run)
    azimuth_rad = None if azimuth_deg is None else np.radians(azimuth_deg)
    result = diffraction(read_s_matrix(run / "s_matrix.tsv"), info, analysis_energies(info).collision_ev, azimuth_rad)
    return [c for c, inside in zip(result.channels, result.in_plane, strict=True) if inside]


def test_plane_at_normal_incidence_follows_the_azimuth(fixture_runs: Path) -> None:
    """At normal incidence the plane is the chosen azimuth: 0 (default) holds the channels (n, 0), 90 deg the
    channels (0, m) of the rectangular cell."""
    assert in_plane_channels(fixture_runs, None) == [(0, 0), (1, 0)]
    assert in_plane_channels(fixture_runs, 90.0) == [(0, 0), (0, 1)]
