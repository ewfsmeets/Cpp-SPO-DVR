"""The S matrix rebuilt from bk_coefficients.bin (copies/s_matrix.py through compute/analysis.py, as the
figures use it) equals the program's s_matrix.tsv on the fixture runs: the check that the BK reader, the
time sums, the prefactor and the phases (each channel's referred to the line where it is measured) are right."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from compute.analysis import analysis_energies, analysis_time_step_au, s_matrix_factors
from compute.time_series import analysis_lines_bohr
from copies.s_matrix import normal_momentum, time_sums
from readers.binary import read_bk_coefficients, read_checkpoint
from readers.run_info import read_run_info
from readers.tables import read_s_matrix

TOLERANCE = 1e-12  # of the largest |S| at each energy: the sums run in another order than the C++'s


def rebuilt_s_matrix(run: Path, slip: str = "") -> tuple[np.ndarray, np.ndarray]:
    """S (n_energies, n_v, n_states) of a run from its BK file, and the BK file's state labels. The slips of
    the negative tests: "collision energy" puts E_coll instead of E_tot into the time phase; "cartesian"
    takes the skewed reciprocal-axis components kx, ky as Cartesian in k0; "zana" refers every phase to Zana
    (the Fortran's convention: the initial state is measured at the specular line)."""
    info = read_run_info(run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    energies = analysis_energies(info)
    if slip == "cartesian":
        packet, mass = info.data["wave_packet"], info.data["molecule"]["total_mass"]
        kx, ky = packet["kx"], packet["ky"]
        k0 = normal_momentum(
            energies.collision_ev,
            total_mass_me=mass,
            hartree_to_ev=info.hartree_to_ev,
            kx_per_bohr=kx,
            ky_per_bohr=ky,
            skew_angle_rad=np.pi / 2,
        )
        energies = replace(energies, k0_per_bohr=k0)
    use_total_energy = slip != "collision energy"
    phase_energy = energies.total_hartree if use_total_energy else energies.collision_ev / info.hartree_to_ev
    n_times, n_v, n_states = bk.coefficients.shape
    sums = time_sums(
        bk.times_au, bk.coefficients.reshape(n_times, n_v * n_states), phase_energy, analysis_time_step_au(info)
    )
    checkpoint = read_checkpoint(run / "checkpoint.bin")
    lines = analysis_lines_bohr(info, bk, checkpoint.surface_index, checkpoint.surface_fixed and slip != "zana")
    s = sums * s_matrix_factors(info, energies, bk, lines)
    return s.reshape(len(energies.collision_ev), n_v, n_states), bk.labels


def worst_relative_difference(run: Path, slip: str = "") -> np.ndarray:
    """Per energy: max over the open channels of s_matrix.tsv of |S_rebuilt - S_file| / max |S_file|."""
    rebuilt, labels = rebuilt_s_matrix(run, slip)
    table = read_s_matrix(run / "s_matrix.tsv")
    collision = analysis_energies(read_run_info(run)).collision_ev
    state_of = {tuple(label[:4]): i for i, label in enumerate(labels)}  # (j, mJ, n, m) -> i_state
    worst = np.zeros(len(collision))
    largest = np.zeros(len(collision))
    for row in range(len(table["E_coll"])):
        energy = int(np.argmin(np.abs(collision - table["E_coll"][row])))
        key = tuple(int(table[c][row]) for c in ("j", "mJ", "n", "m"))
        s_file = complex(table["Re_S"][row], table["Im_S"][row])
        s_here = rebuilt[energy, int(table["v"][row]), state_of[key]]
        worst[energy] = np.max([worst[energy], abs(s_here - s_file)])  # np.max keeps a NaN
        largest[energy] = np.max([largest[energy], abs(s_file)])
    return worst / largest


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique", "dailight_oblique_split/second"])
def test_s_matrix_from_the_bk_file(fixture_runs: Path, run: str) -> None:
    """Every open channel at every energy, |dS| <= 1e-12 of that energy's largest |S| (the second segment of a
    split run holds every record of the run)."""
    relative = worst_relative_difference(fixture_runs / run)
    assert np.max(relative) <= TOLERANCE, f"worst {np.max(relative):.3e} at energy index {np.argmax(relative)}"


@pytest.mark.parametrize("slip", ["collision energy", "cartesian", "zana"])
def test_a_porting_slip_fails_at_every_energy(fixture_runs: Path, slip: str) -> None:
    """Negative tests, three realistic slips: E_coll instead of E_tot in exp(i E t); the skewed kx, ky taken as
    Cartesian in k0 (the oblique fixture has kx 5.7 and ky 4.4 bohr^-1 on the 60 degree cell); the initial state's
    phase referred to Zana as in the Fortran (the fixture's specular line is fixed; measured: all 41 energies fail,
    the smallest by 0.022 of the largest |S|)."""
    relative = worst_relative_difference(fixture_runs / "dailight_oblique", slip)
    assert np.sum(relative > TOLERANCE) == len(relative)
