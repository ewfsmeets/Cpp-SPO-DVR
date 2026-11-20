"""energy_resolved/ (EnergyResolvedSites, EnergyResolvedEnergies) on the fixture runs, each check with its negative
tests:
- the reader: shapes, the steps summed and the times of the final and partial sums;
- the sums (dailight_rho_e: all 16 grid points of the oblique 4 x 4 cell, the analysis energies 0.45 and 0.55 eV):
  projected at the analysis line onto the channels exactly as the program projects the wave packet
  (src/analysis/balint_kurti.cpp), they give back s_matrix.tsv's S of every channel but the initial one (analysed on
  the specular grid); the sums are linear in psi, so this checks the gather, the phases and the bookkeeping;
- the stationary state at oblique incidence (the same run): the plot tool's psi_E projected at every Z between zasymp
  and the valid region's top gives S_f(Z) = s_matrix.tsv's S_f for every open non-initial channel;
- the stationary state of the hard wall (hardwall_long): the initial channel's F(Z) against an exact Numerov solution
  exp(-i k' Z) - S_N exp(i k' Z) (S_N the Numerov reflection coefficient, not s_matrix.tsv's convention) of the wall
  V_z at the normal energy the split operator propagates, E_n - delta (delta: the effective-energy shift of the
  Strang r propagator for phi_v; k' = sqrt(2 M (E_n - delta))), phase included;
- the figure's content.
"""

from pathlib import Path

import numpy as np
import pytest

from compute.analysis import analysis_energies, analysis_time_step_au, s_matrix_factors
from compute.energy_resolved import in_valid_z, stationary_state, valid_region
from copies.s_matrix import time_sums
from figures import default_options, figure_energy_resolved
from readers.binary import read_bk_coefficients
from readers.energy_resolved import read_energy_resolved
from readers.run_info import RunInfo, read_run_info
from readers.tables import read_asymptotic_potential, read_gas_phase_states, read_s_matrix
from tests.absorber_profiles import absorber_factors, absorbers_on

S_TOLERANCE = 1e-13  # of the largest |S| at each energy: the same sums in another order (measured 1.3e-15)


def test_reader_reads_the_sums(fixture_runs: Path) -> None:
    """hardwall_long: one site, three energies; 2000 analysis steps summed (tEnd 20000, dtAna 10), the final sums
    at 20000 and the partial ones at the first step at or after 0.8 tEnd, 16000 a.u.; the partial sums differ."""
    sums = read_energy_resolved(fixture_runs / "hardwall_long")
    assert [s.name for s in sums.sites] == ["origin"] and list(sums.collision_ev) == [0.3, 0.4, 0.5]
    assert sums.sums.shape[:2] == (1, 3) and sums.specular_sums.shape[0] == 3
    assert sums.steps_summed == 2000 and sums.last_time_au == 20000.0 and sums.partial_time_au == 16000.0
    assert not np.array_equal(np.asarray(sums.sums), np.asarray(sums.sums_partial))


# --- the S-matrix gate -------------------------------------------------------------------------------------------


def channel_sums(run: Path, slip: str = "") -> tuple[np.ndarray, list[int]]:
    """(n_energies_of_the_sums, n_v * n_states): the BK projection of the sums at the analysis line (Zana) for every
    analysed channel (index i_state + n_states v, as bk_coefficients.bin), and the indices of the sums' energies
    among the analysis energies. Slips of the negative tests: "no shift" leaves the incident parallel momentum in
    (shift_forward = 1); "transposed" takes the sites' grid points as (iy, ix)."""
    info, sums = read_run_info(run), read_energy_resolved(run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    gas = read_gas_phase_states(run / "gas_phase_states.dat")
    g, nx, ny = info.data["grids"], info.data["grids"]["nx"], info.data["grids"]["ny"]
    iz = int(round((info.data["analysis"]["z_ana"] - g["z_start"]) / g["dz"]))
    plane = np.zeros((len(sums.collision_ev), nx, ny, len(sums.r_bohr), sums.sums.shape[4]), dtype=complex)
    for k, site in enumerate(sums.sites):
        shift = np.conj(site.specular_to_scattering) / g["unit_cell_norm"]  # shift_forward = conj(s2s) / norm
        ix, iy = (site.iy, site.ix) if slip == "transposed" else (site.ix, site.iy)
        plane[:, ix, iy] = (1.0 if slip == "no shift" else shift) * np.asarray(sums.sums[k, :, iz])
    states = [tuple(s) for s in info.data["rotational_basis"]["states"]]
    xy_scale = np.sqrt(g["dx"] * g["dy"] * np.sin(info.data["lattice"]["skew_angle_rad"]) / (nx * ny))
    n_v, n_states = bk.channel_energy_hartree.shape
    out = np.zeros((len(sums.collision_ev), n_v * n_states), dtype=complex)
    for i_state, (j, m_j, n, m, _) in enumerate(bk.labels):
        irot = states.index((j, m_j))
        for v in range(n_v):
            phi = gas.wave_functions[np.flatnonzero((gas.v == v) & (gas.j == j))[0]]
            overlap = np.einsum("exyr,r->exy", plane[..., irot], phi) * g["dr"]
            momentum = np.fft.fft2(overlap, axes=(1, 2))  # exponent sign minus, unscaled (FftSign::Minus, scale 1)
            out[:, i_state + n_states * v] = momentum[:, n % nx, m % ny] * xy_scale
    collision = analysis_energies(info).collision_ev
    indices = [int(np.argmin(np.abs(collision - e))) for e in sums.collision_ev]
    return out, indices


def s_from_the_sums(run: Path, slip: str = "") -> np.ndarray:
    """(n_energies_of_the_sums, n_v * n_states): S of every channel from the sums (channel_sums times the S matrix's
    prefactor, phases and momenta, compute/analysis.py s_matrix_factors)."""
    info = read_run_info(run)
    projected, indices = channel_sums(run, "" if slip == "missing step" else slip)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    energies = analysis_energies(info)
    if slip == "missing step":  # the sums without the analysis step in the middle of the run
        middle = len(bk.times_au) // 2
        n_times, n_v, n_states = bk.coefficients.shape
        one_step = time_sums(
            bk.times_au[middle : middle + 1],
            bk.coefficients[middle : middle + 1].reshape(1, n_v * n_states),
            energies.total_hartree,
            analysis_time_step_au(info),
        )
        projected = projected - one_step[indices]
    return projected * s_matrix_factors(info, energies, bk)[indices]


def s_differences(run: Path, slip: str = "") -> np.ndarray:
    """Per energy of the sums: max over s_matrix.tsv's open channels but the initial one of |S_sums - S_file|, over
    that energy's largest |S_file|."""
    info = read_run_info(run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    rebuilt = s_from_the_sums(run, slip)
    sums = read_energy_resolved(run)
    table = read_s_matrix(run / "s_matrix.tsv")
    state_of = {tuple(label[:4]): i for i, label in enumerate(bk.labels)}
    initial = (info.data["initial_state"]["j"], info.data["initial_state"]["mJ"], 0, 0)
    n_states = bk.labels.shape[0]
    worst, largest = np.zeros(len(sums.collision_ev)), np.zeros(len(sums.collision_ev))
    for row in range(len(table["E_coll"])):
        matches = np.flatnonzero(np.abs(sums.collision_ev - table["E_coll"][row]) < 1e-9)
        key = tuple(int(table[c][row]) for c in ("j", "mJ", "n", "m"))
        if len(matches) == 0 or key == initial:
            continue
        e = int(matches[0])
        s_file = complex(table["Re_S"][row], table["Im_S"][row])
        s_here = rebuilt[e, state_of[key] + n_states * int(table["v"][row])]
        worst[e] = np.max([worst[e], abs(s_here - s_file)])  # np.max keeps a NaN
        largest[e] = np.max([largest[e], abs(s_file)])
    return worst / largest


def test_sums_give_back_the_s_matrix(fixture_runs: Path) -> None:
    """dailight_rho_e at both energies: every non-initial channel of s_matrix.tsv to S_TOLERANCE of the largest |S|."""
    relative = s_differences(fixture_runs / "dailight_rho_e")
    assert len(relative) == 2 and np.max(relative) <= S_TOLERANCE, relative


@pytest.mark.parametrize("slip", ["no shift", "transposed", "missing step"])
def test_a_slip_in_the_projection_fails_at_both_energies(fixture_runs: Path, slip: str) -> None:
    """Negative tests: the incident parallel momentum left in (oblique incidence: kx 5.7, ky 4.4 bohr^-1), the sites'
    grid points transposed (the 60 degree cell is not symmetric under the swap), and one analysis step (of 1400)
    missing from the sums."""
    relative = s_differences(fixture_runs / "dailight_rho_e", slip)
    assert np.sum(relative > S_TOLERANCE) == 2


# --- the stationary state's channels at every Z (oblique incidence) -----------------------------------------------


def stationary_s(run: Path, slip: str = "", partial: bool = False) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """S_f(Z) of every channel from the plot tool's psi_E (compute/energy_resolved.py stationary_state at all 16 sites)
    projected at each Z point of [zasymp, the valid region's top) as the program projects psi (balint_kurti.cpp):
    S_f(Z) = -Phi_f(Z) exp(-i k_f Z) exp(-i (k0 - kz) Z0) sqrt(k_f / k0), which equals S_f at every such Z for the
    stationary state (s_matrix.cpp: arg = -k_f Zana + (kz - k0) Z0 - pi). Returns the Z points, S (n_energies, n_z,
    n_v * n_states) and the indices of the sums' energies among the analysis energies. Slips: "no weight" leaves out
    psi_E's 1 / sqrt(weight); "k0 for kf" puts k0 into the phase exp(-i k_f Z); "into ZOpt" extends the Z range to
    the Z absorber's end. partial: from the partial sums (up to 0.8 tEnd)."""
    info, sums = read_run_info(run), read_energy_resolved(run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    gas = read_gas_phase_states(run / "gas_phase_states.dat")
    g, nx, ny = info.data["grids"], info.data["grids"]["nx"], info.data["grids"]["ny"]
    packet, mass = info.data["wave_packet"], info.data["molecule"]["total_mass"]
    region = valid_region(info)
    z_top = (
        next(a["end"] for a in info.data["absorbers"]["list"] if a["name"] == "ZOpt")
        if slip == "into ZOpt"
        else (region.z_high_bohr)
    )
    rows = np.flatnonzero((sums.z_bohr >= info.data["z_asymptotic"]) & (sums.z_bohr < z_top))
    energies = analysis_energies(info)
    indices = [int(np.argmin(np.abs(energies.collision_ev - e))) for e in sums.collision_ev]
    k0 = energies.k0_per_bohr[indices]
    weight = np.exp(-((packet["kz"] - k0) ** 2) / (2 * packet["alpha"]))
    states = [tuple(x) for x in info.data["rotational_basis"]["states"]]
    xy_scale = np.sqrt(g["dx"] * g["dy"] * np.sin(info.data["lattice"]["skew_angle_rad"]) / (nx * ny))
    n_v, n_states = bk.channel_energy_hartree.shape
    out = np.zeros((len(indices), len(rows), n_v * n_states), dtype=complex)
    for e, ie in enumerate(indices):
        plane = np.zeros((len(rows), nx, ny, len(sums.r_bohr), sums.sums.shape[4]), dtype=complex)
        for k, site in enumerate(sums.sites):
            psi = stationary_state(sums, info, k, e, partial)[rows] * (
                np.sqrt(weight[e]) if slip == "no weight" else 1.0
            )
            plane[:, site.ix, site.iy] = np.conj(site.specular_to_scattering) / g["unit_cell_norm"] * psi
        k_f = np.sqrt(np.maximum(2 * mass * (energies.total_hartree[ie] - bk.channel_energy_hartree), 0.0))
        for i_state, (j, m_j, n, m, _) in enumerate(bk.labels):
            irot = states.index((j, m_j))
            for v in range(n_v):
                phi = gas.wave_functions[np.flatnonzero((gas.v == v) & (gas.j == j))[0]]
                overlap = np.einsum("zxyr,r->zxy", plane[..., irot], phi) * g["dr"]
                amplitude = np.fft.fft2(overlap, axes=(1, 2))[:, n % nx, m % ny] * xy_scale
                kf = k0[e] if slip == "k0 for kf" else k_f[v, i_state]
                phase = np.exp(-1j * kf * sums.z_bohr[rows]) * np.exp(-1j * (k0[e] - packet["kz"]) * packet["z0"])
                out[e, :, i_state + n_states * v] = -amplitude * phase * np.sqrt(k_f[v, i_state] / k0[e])
    return sums.z_bohr[rows], out, indices


def stationary_s_differences(run: Path, slip: str = "") -> np.ndarray:
    """Per energy: max over Z and over s_matrix.tsv's open channels but the initial one of |S_f(Z) - S_tsv|, over that
    energy's largest |S_tsv|."""
    info = read_run_info(run)
    bk = read_bk_coefficients(run / "bk_coefficients.bin")
    _, s_of_z, indices = stationary_s(run, slip)
    table = read_s_matrix(run / "s_matrix.tsv")
    collision = analysis_energies(info).collision_ev
    state_of = {tuple(label[:4]): i for i, label in enumerate(bk.labels)}
    initial = (info.data["initial_state"]["j"], info.data["initial_state"]["mJ"], 0, 0)
    n_states = bk.labels.shape[0]
    worst, largest = np.zeros(len(indices)), np.zeros(len(indices))
    for row in range(len(table["E_coll"])):
        ie = int(np.argmin(np.abs(collision - table["E_coll"][row])))
        key = tuple(int(table[c][row]) for c in ("j", "mJ", "n", "m"))
        if ie not in indices or key == initial:
            continue
        e = indices.index(ie)
        s_file = complex(table["Re_S"][row], table["Im_S"][row])
        s_here = s_of_z[e, :, state_of[key] + n_states * int(table["v"][row])]
        worst[e] = np.max([worst[e], np.max(np.abs(s_here - s_file))])  # np.max keeps a NaN
        largest[e] = np.max([largest[e], abs(s_file)])
    return worst / largest


# The phase drift of the split operator's propagated momentum against k_f over the checked Z range: delta / v =
# 6.65e-7 hartree x M / k ~ 2.3e-4 bohr^-1 (dtStep 2) over at most 1.2 bohr of Z, rounded up.
SPLIT_OPERATOR_DRIFT = 3e-4


def test_stationary_channels_hold_the_s_matrix_at_every_z(fixture_runs: Path) -> None:
    """dailight_rho_e at both energies: S_f(Z) of the plot tool's psi_E equals s_matrix.tsv's S_f at every Z between
    zasymp (8.0) and the valid region's top (9.2 bohr) for every open non-initial channel, within the truncation
    indicator of the same quantity (the partial sums against the final ones: 2.7e-3, 1.4e-3 of max |S|) plus
    SPLIT_OPERATOR_DRIFT; measured 1.1e-3, 7.1e-4 (2026-11-20). The only check of psi_E at oblique incidence (k0 from
    the parallel momentum, the weight) and of the valid region's top."""
    run = fixture_runs / "dailight_rho_e"
    relative = stationary_s_differences(run)
    _, final, _ = stationary_s(run)
    _, partial, _ = stationary_s(run, partial=True)
    largest = np.max(np.abs(final), axis=(1, 2))
    tolerance = np.max(np.abs(final - partial), axis=(1, 2)) / largest + SPLIT_OPERATOR_DRIFT
    assert len(relative) == 2 and np.all(relative <= tolerance), (relative, tolerance)


@pytest.mark.parametrize("slip", ["no weight", "k0 for kf", "into ZOpt"])
def test_a_slip_in_the_stationary_channels_fails_at_both_energies(fixture_runs: Path, slip: str) -> None:
    """Negative tests: psi_E without its 1 / sqrt(weight) (weight 0.77 at 0.55 eV), k0 in place of k_f in the outgoing
    phase, and the Z range extended into the Z absorber: each exceeds 4e-3 of max |S| (above both tolerances) at both
    energies."""
    assert np.sum(stationary_s_differences(fixture_runs / "dailight_rho_e", slip) > 4e-3) == 2


# --- the hard-wall gate -----------------------------------------------------------------------------------------

WALL_Z_MAX, WALL_Z_MIN = 7.0, 0.1778  # the hard wall's cosine switch (src/potential/hard_wall_potential.cpp)
NUMEROV_START_BOHR = 0.3  # inside the wall (V_z ~ 120 eV): the solution starts at 0 there
NUMEROV_REFINE = 32  # Numerov steps per scattering-grid step (its error ~ h^4: negligible)


def wall_hartree(z_bohr: np.ndarray, info: RunInfo) -> np.ndarray:
    """The hard wall's Z term V_z (hartree): 1 / (Z^4 + 0.001) eV with its cosine switch to 0 at Z = 7 bohr
    (Shadows: src/potential/hard_wall_potential.cpp HardWallPotential::value; tests/test_runs_potential.py)."""
    v_z = 1.0 / (z_bohr * z_bohr * z_bohr * z_bohr + 0.001)
    switch = 1.0 - (0.5 - 0.5 * np.cos(np.pi * (z_bohr - WALL_Z_MIN) / (WALL_Z_MAX - WALL_Z_MIN)))
    factor = np.where(z_bohr < WALL_Z_MIN, 1.0, switch)
    return np.where(z_bohr < WALL_Z_MAX, factor * v_z, 0.0) / info.hartree_to_ev


def strang_shift_hartree(run: Path, info: RunInfo) -> float:
    """delta: the effective energy of phi_v under one Strang step of the r motion, exp(-i T dt/2) [ROpt]
    exp(-i V_asy dt) exp(-i T dt/2) on numpy's FFT grid, minus its FGH energy: -arg <phi|U|phi> / dt - E_v."""
    gas = read_gas_phase_states(run / "gas_phase_states.dat")
    state = np.flatnonzero((gas.v == info.data["initial_state"]["v"]) & (gas.j == info.data["initial_state"]["j"]))[0]
    phi, r = gas.wave_functions[state], gas.r_bohr
    dt, mu = info.data["times"]["dt_step"], info.data["molecule"]["reduced_mass"]
    k = 2 * np.pi * np.fft.fftfreq(len(r), r[1] - r[0])
    half_kinetic = np.exp(-1j * k**2 / (2 * mu) * dt / 2)
    potential = absorber_factors(info, absorbers_on(info, "r", "scattering"), r) * np.exp(
        -1j * read_asymptotic_potential(run / "asymptotic_potential.dat").v_hartree * dt
    )
    stepped = np.fft.ifft(half_kinetic * np.fft.fft(potential * np.fft.ifft(half_kinetic * np.fft.fft(phi))))
    amplitude = np.sum(phi * stepped) / np.sum(phi * phi)
    return float(-np.angle(amplitude) / dt - gas.energy_hartree[state])


def numerov_wall_state(z_bohr: np.ndarray, energy_hartree: float, mass: float, info: RunInfo) -> np.ndarray:
    """The stationary state of the wall at a normal energy, exp(-i k Z) + S exp(i k Z) where V_z = 0 (Z >= 7), on the
    points z_bohr: Numerov outward from 0 deep inside the wall on a grid NUMEROV_REFINE times finer, the incoming
    amplitude fitted at two points beyond the wall and divided out."""
    dz = (z_bohr[1] - z_bohr[0]) / NUMEROV_REFINE
    fine = np.arange(NUMEROV_START_BOHR, z_bohr[-1] + dz, dz)
    q = 2 * mass * (energy_hartree - wall_hartree(fine, info))
    f = np.zeros(len(fine))
    f[1] = 1e-30
    c = dz * dz / 12
    for n in range(1, len(fine) - 1):
        f[n + 1] = (2 * f[n] * (1 - 5 * c * q[n]) - f[n - 1] * (1 + c * q[n - 1])) / (1 + c * q[n + 1])
    k = np.sqrt(2 * mass * energy_hartree)
    a, b = len(fine) - 1, len(fine) - 1 - NUMEROV_REFINE  # two points beyond the wall
    basis = np.array(
        [[np.exp(-1j * k * fine[a]), np.exp(1j * k * fine[a])], [np.exp(-1j * k * fine[b]), np.exp(1j * k * fine[b])]]
    )
    incoming, _ = np.linalg.solve(basis, np.array([f[a], f[b]], dtype=complex))
    state = f / incoming
    return np.where(
        z_bohr < NUMEROV_START_BOHR, 0.0, np.interp(z_bohr, fine, state.real) + 1j * np.interp(z_bohr, fine, state.imag)
    )


def wall_amplitudes(run: Path, energy: int, slip: str = "") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """On the valid region's Z points above NUMEROV_START_BOHR: the points, F(Z) of the initial channel at one energy,
    the same of the partial sums, and the Numerov state. F = exp(-i (k' - kz) Z0) sqrt(A) sum_r phi_v(r)
    psi_E(Z, r, initial) dr: the packet's phase at Z0 divided out with the momentum k' = sqrt(2 M (E_n - delta))
    with which its component at E_tot moves (k_par = 0 here). Slips of the negative tests: "no k0" leaves out
    psi_E's k0 factor, "total" solves at E_tot - delta (E(vin, jin) not subtracted), "no delta" leaves delta out,
    "scale" makes F 1.001 times too large (a normalisation slip). F is rescaled by (k'/k0) sqrt(weight(k0) /
    weight(k')) to the momentum k' of the state at E_tot (psi_E's prefactor uses k0, as the S matrix)."""
    info, sums = read_run_info(run), read_energy_resolved(run)
    gas = read_gas_phase_states(run / "gas_phase_states.dat")
    state = np.flatnonzero((gas.v == info.data["initial_state"]["v"]) & (gas.j == info.data["initial_state"]["j"]))[0]
    packet, mass = info.data["wave_packet"], info.data["molecule"]["total_mass"]
    initial = int(info.data["rotational_basis"]["initial_index"])
    collision_hartree = sums.collision_ev[energy] / info.hartree_to_ev
    delta = 0.0 if slip == "no delta" else strang_shift_hartree(run, info)
    energy_hartree = collision_hartree - delta
    k0, k_prime = np.sqrt(2 * mass * collision_hartree), np.sqrt(2 * mass * energy_hartree)
    phase = np.exp(-1j * (k_prime - packet["kz"]) * packet["z0"])

    # psi_E's prefactor takes k0 and weight(k0) (as the S matrix); the state at E_tot moves with k': rescale to it.
    def packet_weight_at(k: float) -> float:
        """exp(-(kz - k)^2 / (2 alpha))."""
        return float(np.exp(-((packet["kz"] - k) ** 2) / (2 * packet["alpha"])))

    consistent = k_prime / k0 * np.sqrt(packet_weight_at(k0) / packet_weight_at(k_prime))
    weight = np.sqrt(info.data["grids"]["unit_cell_area"]) * info.data["grids"]["dr"] * consistent
    if slip == "scale":
        weight = weight * 1.001
    amplitudes = []
    for partial in (False, True):
        psi = stationary_state(sums, info, 0, energy, partial)[:, :, initial]
        if slip == "no k0":
            psi = psi / np.sqrt(2 * mass * collision_hartree)
        amplitudes.append(phase * weight * (psi @ gas.wave_functions[state]))
    if slip == "total":
        energy_hartree = collision_hartree + gas.energy_hartree[state] - delta
    exact = numerov_wall_state(sums.z_bohr, energy_hartree, mass, info)
    region = valid_region(info)
    shown = (sums.z_bohr > max(region.z_low_bohr, NUMEROV_START_BOHR)) & (sums.z_bohr < region.z_high_bohr)
    return sums.z_bohr[shown], amplitudes[0][shown], amplitudes[1][shown], exact[shown]


def amplitudes_beyond_the_wall(z_bohr: np.ndarray, f: np.ndarray, k_per_bohr: float) -> tuple[complex, complex]:
    """The incoming and outgoing amplitudes A, B of F = A exp(-i k Z) + B exp(i k Z) fitted where V_z = 0 (Z >= 7)."""
    beyond = z_bohr >= WALL_Z_MAX
    waves = np.stack([np.exp(-1j * k_per_bohr * z_bohr[beyond]), np.exp(1j * k_per_bohr * z_bohr[beyond])], axis=1)
    (incoming, outgoing), *_ = np.linalg.lstsq(waves, f[beyond], rcond=None)
    return complex(incoming), complex(outgoing)


def wall_check(run: Path, energy: int, slip: str = "") -> tuple[float, float, float, float]:
    """(max |F - F_exact|, its tolerance, ||A| - 1|, its tolerance). The truncation indicator tau = max |F - F_partial|
    (what the time sums still gained after 0.8 tEnd: an indicator, not a bound) enters both. F's tolerance adds the
    absorbers' non-unitarity |F|max ||B| / |A| - 1| (independent of F's scale); |A| - 1 (the normalisation, 1 for the
    exact state) is allowed tau / |F|max (the indicator relative to the state) plus the |B|/|A| term, both measured on
    the correct F (the slips change F, not the tolerances)."""
    z, f, f_partial, exact = wall_amplitudes(run, energy, slip)
    z, f_right, f_right_partial, _ = wall_amplitudes(run, energy)
    info, sums = read_run_info(run), read_energy_resolved(run)
    e_n = sums.collision_ev[energy] / info.hartree_to_ev - strang_shift_hartree(run, info)
    k = np.sqrt(2 * info.data["molecule"]["total_mass"] * e_n)
    largest = float(np.max(np.abs(f_right)))
    tau = float(np.max(np.abs(f_right - f_right_partial)))
    incoming, outgoing = amplitudes_beyond_the_wall(z, f_right, k)
    unitarity = abs(abs(outgoing) / abs(incoming) - 1)
    measured_incoming, _ = amplitudes_beyond_the_wall(z, f, k)
    return (
        float(np.max(np.abs(f - exact))),
        tau + largest * unitarity,
        abs(abs(measured_incoming) - 1),
        tau / largest + unitarity,
    )


@pytest.mark.parametrize("energy", [0, 1, 2])
def test_hard_wall_stationary_state_is_the_exact_one(fixture_runs: Path, energy: int) -> None:
    """At 0.3, 0.4 and 0.5 eV, phase included: max |F - F_exact| and ||A| - 1| within their derived tolerances
    (wall_check), which stay below 1e-2 (beyond that the comparison would say nothing at |F| up to 2.4). delta, the
    Strang shift at dtStep 2 on this r grid, is 6.65e-7 hartree."""
    run = fixture_runs / "hardwall_long"
    assert strang_shift_hartree(run, read_run_info(run)) == pytest.approx(6.652e-7, rel=1e-3)
    error, tolerance, incoming, incoming_tolerance = wall_check(run, energy)
    assert tolerance < 1e-2 and error <= tolerance, (error, tolerance)
    assert incoming <= incoming_tolerance, (incoming, incoming_tolerance)


@pytest.mark.parametrize("slip", ["no k0", "total", "no delta", "scale"])
def test_a_slip_in_the_stationary_state_fails_at_every_energy(fixture_runs: Path, slip: str) -> None:
    """Negative tests: psi_E without its k0 factor, the Numerov energy E_tot instead of E_n, the Strang shift delta
    left out (its phase (k0 - k') Z0 = 2.8e-3 rad gives 5e-3 at |F| = 2.4), and F 1.001 times too large each fail at
    all three energies (in F or in |A|)."""
    run = fixture_runs / "hardwall_long"
    checks = [wall_check(run, e, slip) for e in range(3)]
    failing = [error > tolerance or incoming > bound for error, tolerance, incoming, bound in checks]
    assert sum(failing) == 3, checks


def test_figure_content(fixture_runs: Path) -> None:
    """dailight_rho_e: one row per energy (2); left the map above the first site with its colour bar, drawn only in the
    valid region (NaN at Z outside it and at r beyond the r absorber's start); right the profiles of the first six of
    the 16 sites (NaN outside the valid Z) with the valid region in the legend and, in the title, the change after the
    partial sums' time (0.8 tEnd = 11200 a.u.); the notes say that it is a density and that 6 of 16 sites are drawn."""
    run = fixture_runs / "dailight_rho_e"
    info = read_run_info(run)
    (stem, fig), *_ = figure_energy_resolved([run], info, default_options())
    assert stem == "energy_resolved"
    panels = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
    assert len(panels) == 4 and sum(ax.get_label() == "<colorbar>" for ax in fig.axes) == 2
    titles = [ax.get_title() for ax in panels[0::2]]
    assert titles == ["E_coll 0.45 eV: A rho_E above g00", "E_coll 0.55 eV: A rho_E above g00"]
    sums, region = read_energy_resolved(run), valid_region(info)
    inside = in_valid_z(sums.z_bohr, region)[:, None] & (sums.r_bohr < region.r_high_bohr)[None, :]
    for ax in panels[0::2]:
        drawn = np.ma.filled(ax.collections[0].get_array(), np.nan).reshape(len(sums.r_bohr), len(sums.z_bohr)).T
        assert np.array_equal(np.isfinite(drawn), inside)
    for ax in panels[1::2]:
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert labels == ["g00", "g01", "g02", "g03", "g10", "g11", "valid region"]
        assert "changed after t = 11200 a.u. by up to" in ax.get_title()
        for line in ax.get_lines()[:6]:
            assert np.array_equal(np.isfinite(line.get_ydata()), in_valid_z(sums.z_bohr, region))
    info_text = fig.texts[0].get_text()
    assert "not a reaction probability" in info_text and "profiles of the first 6 of 16 sites" in info_text
