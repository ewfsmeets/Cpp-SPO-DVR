"""The results on the fixture runs: the channel probabilities of s_matrix.tsv sum to P_scattered; the
computed results, the time series of P_scattered and the gas-phase levels agree with the program's files."""

from pathlib import Path

import numpy as np
import pytest

from compute.gas_phase import levels
from compute.results import results
from compute.time_series import convergence_ratio, p_scattered_stopped_at, scattering_in_time
from readers.binary import read_bk_coefficients
from readers.run_info import read_run_info
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_s_matrix, read_scattering_summary

SUM_TOLERANCE = 1e-13  # relative: the C++ and numpy add the same few hundred terms in other orders
S_TOLERANCE = 1e-12  # relative: as tests/test_runs_s_matrix.py


def summed_probabilities(run: Path, double_count: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """(sum of P over the open channels of each energy, P_scattered); double_count adds the first channel twice
    (the slip of the negative test)."""
    s = read_s_matrix(run / "s_matrix.tsv")
    summary = read_scattering_summary(run / "scattering_summary.tsv")
    sums = np.array([np.sum(s["P"][s["E_coll"] == e]) for e in summary["E_coll"]])
    if double_count:
        firsts = np.array([s["P"][s["E_coll"] == e][0] for e in summary["E_coll"]])
        sums = sums + firsts
    return sums, summary["P_scattered"]


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique"])
def test_channels_sum_to_p_scattered(fixture_runs: Path, run: str) -> None:
    """Sum over channels = P_scattered at every energy, to SUM_TOLERANCE relative."""
    sums, p = summed_probabilities(fixture_runs / run)
    assert np.max(np.abs(sums - p) / p) <= SUM_TOLERANCE


def test_a_channel_counted_twice_fails(fixture_runs: Path) -> None:
    """Negative test: the first channel counted twice fails at every energy where it carries probability."""
    sums, p = summed_probabilities(fixture_runs / "dailight_oblique", double_count=True)
    assert np.sum(np.abs(sums - p) / p > SUM_TOLERANCE) == len(p)


@pytest.mark.parametrize("run", ["hardwall_tiny", "dailight_oblique"])
def test_p_scattered_in_time_ends_at_the_programs(fixture_runs: Path, run: str) -> None:
    """P_scattered(E; tEnd) of the time series (sums grown record by record) and of the sums stopped at the
    last record both equal scattering_summary.tsv's; the convergence ratio is 1 at tEnd."""
    info = read_run_info(fixture_runs / run)
    bk = read_bk_coefficients(fixture_runs / run / "bk_coefficients.bin")
    expected = read_scattering_summary(fixture_runs / run / "scattering_summary.tsv")["P_scattered"]
    series = scattering_in_time(info, bk)
    assert np.max(np.abs(series.p_scattered[-1] - expected) / expected) <= S_TOLERANCE
    stopped = p_scattered_stopped_at(info, bk, len(bk.times_au) - 1)
    assert np.max(np.abs(stopped - expected) / expected) <= S_TOLERANCE
    assert np.max(np.abs(convergence_ratio(series)[:, -1] - 1.0)) == 0.0


def test_results_of_the_oblique_run(fixture_runs: Path) -> None:
    """dailight_oblique: the unitarity error is -min P_reaction; two levels are scattered into (v' = 1 opens
    at the highest energies)."""
    summary = read_scattering_summary(fixture_runs / "dailight_oblique" / "scattering_summary.tsv")
    r = results(summary)
    assert r.unitarity_error == -np.min(summary["P_reaction"]) > 0
    assert sorted(r.p_vibration) == [0, 1]


def test_gas_phase_levels_bound_below_the_last_potential_point(fixture_runs: Path) -> None:
    """A level is bound when its energy is below V_asy at the last r point; the ground state is bound and its
    energy is the initial state's E(v_in, j_in) of run_info.json."""
    run = fixture_runs / "hardwall_tiny"
    info = read_run_info(run)
    potential = read_asymptotic_potential(run / "asymptotic_potential.dat")
    gas = levels(read_gas_phase_states(run / "gas_phase_states.dat"), potential, info.hartree_to_ev)
    assert list(gas.bound) == list(gas.energy_ev < potential.v_ev[-1])
    ground = int(np.argmin(gas.energy_ev))
    assert gas.bound[ground] and (gas.v[ground], gas.j[ground]) == (0, 0)
    initial_ev = info.data["initial_state"]["energy_ev"]
    assert abs(gas.energy_ev[ground] - initial_ev) <= 1e-12 * abs(initial_ev)
