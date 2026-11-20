"""The FrameExtras datasets of frames/ on the fixture runs, each check with its negative test:
- at every record they sum to the frame's probability (the momenta: to the grids' norms);
- on the checkpoints (at the end of a run, and mid-run while the packet is at the surface) each equals an independent
  numpy computation from psi_full: the angular grid and the normalised associated Legendre functions from numpy's
  Legendre series, the parallel momenta from numpy's FFT of exp(-i k_par.R) psi_full, the vibrational overlaps from
  gas_phase_states.dat, the momenta from numpy's FFTs, the P2 elements by a Gauss-Legendre quadrature;
- <P2> agrees with the DVR quadrature of the angular extra within its known error;
- on the hard wall (no x, y or angle coupling; V separable in Z and r) nothing leaves the (0, 0) order, and the
  vibrational populations follow a numpy split-operator propagation of the r motion alone;
- a split run joins to the unsplit run; the reader refuses manifests it cannot index."""

from math import factorial
from pathlib import Path

import numpy as np
import pytest

from readers.densities import read_densities
from readers.errors import ReadError
from readers.frames import DATASETS, EXTRA_DATASETS, Frames, read_frames
from readers.run_info import RunInfo, read_run_info
from readers.tables import read_asymptotic_potential, read_gas_phase_states
from tests.absorber_profiles import absorber_factors, absorbers_on
from tests.frame_helpers import psi_full_of, relative_misses, run_with_frames_manifest, s2s_of

RUNS = ["hardwall_tiny", "dailight_oblique"]
CHECKPOINTS = ["hardwall_tiny", "dailight_oblique", "dailight_oblique_split/first"]
SUM_TOLERANCE = 1e-13  # relative to the scale named at each use: sums in another order than the program's
# relative to the record's largest value: numpy's FFTs and einsum sum over more terms, in another order
NUMPY_TOLERANCE = 1e-12
QUADRATURE_POINTS = 64  # Gauss-Legendre points of the P2 elements: exact up to degree 127, i.e. for JMax <= 62
MANY = 10  # a realistic slip in a negative test fails at more than this many points (the fixtures' grids give
#            hundreds; one failing point would be an accident)


def last_record(frames: Frames, name: str) -> np.ndarray:
    """The last record of a dataset: the one at the checkpoint's time."""
    return np.asarray(frames.series.arrays[name][-1])


def cell_volume_bohr4(info: RunInfo) -> float:
    """dZ dr dx dy sin(gamma): the volume of a scattering-grid cell."""
    g = info.data["grids"]
    return g["dz"] * g["dr"] * g["xy_volume_element"]


@pytest.mark.parametrize("run", RUNS)
def test_extras_sum_to_the_frames(fixture_runs: Path, run: str) -> None:
    """Every record: at each Z the angular, diffraction and vibration populations sum to pop_jmj_z's sum over the
    states (1e-13 of the largest Z sum of the run); prob_kr sums to prob_zrxy's total, prob_kz_spec to the densities'
    specular cells (1e-13 of their largest; dtAve = dtFrames). Negative (one vibration value off by 10 times the
    tolerance, or NaN) fails exactly its Z point."""
    frames = read_frames([fixture_runs / run])
    arrays = frames.series.arrays
    per_z = np.asarray(arrays["pop_jmj_z"]).sum(axis=2)  # (records, nz)
    scale = float(per_z.max())
    for name in ("prob_angle_z", "pop_nm_z", "pop_vj_z"):
        summed = np.asarray(arrays[name]).sum(axis=(2, 3))
        assert relative_misses(summed, per_z, scale, SUM_TOLERANCE) == 0, name
    total = np.asarray(arrays["prob_zrxy"]).sum(axis=(1, 2, 3, 4))
    assert relative_misses(np.asarray(arrays["prob_kr"]).sum(axis=1), total, total.max(), SUM_TOLERANCE) == 0
    densities, _, _ = read_densities([fixture_runs / run])
    specular_norm = np.asarray(densities.arrays["prob_spec_zr"]).sum(axis=(1, 2))
    by_k_z = np.asarray(arrays["prob_kz_spec"]).sum(axis=1)
    assert relative_misses(by_k_z, specular_norm, specular_norm.max(), SUM_TOLERANCE) == 0
    vibration = np.array(arrays["pop_vj_z"])
    for planted in (vibration[3, 20, 0, 0] + 10 * SUM_TOLERANCE * scale, np.nan):
        wrong = vibration.copy()
        wrong[3, 20, 0, 0] = planted
        assert relative_misses(wrong.sum(axis=(2, 3)), per_z, scale, SUM_TOLERANCE) == 1


def normalised_legendre(j: int, m: int, x: np.ndarray) -> np.ndarray:
    """The program's normalised associated Legendre function of x = cos(theta): sqrt((2J + 1) / 2 (J - |m|)! /
    (J + |m|)!) (1 - x^2)^(|m|/2) d^|m| P_J / dx^|m|, times (-1)^m for m > 0 (the Condon-Shortley phase only there),
    from numpy's Legendre series (independent of the program's recursion)."""
    abs_m = abs(m)
    series = np.zeros(j + 1)
    series[j] = 1.0
    derivative = np.polynomial.legendre.legder(series, abs_m)
    norm = np.sqrt((2 * j + 1) / 2 * factorial(j - abs_m) / factorial(j + abs_m))
    sign = (-1) ** abs_m if m > 0 else 1
    return sign * norm * (1 - x**2) ** (abs_m / 2) * np.polynomial.legendre.legval(x, derivative)


def angular_nodes(j_max: int, homonuclear: bool) -> tuple[np.ndarray, np.ndarray]:
    """theta (rad) and weights of the DVR: the j_max + 1 Gauss-Legendre roots in cos(theta), from theta = 0; for a
    homonuclear molecule those with cos(theta) >= 0 only, the weight doubled where the mirror node is another one."""
    x, w = np.polynomial.legendre.leggauss(j_max + 1)
    x, w = x[::-1], w[::-1]
    if homonuclear:
        keep = x > -1e-14
        x, w = x[keep], np.where(np.abs(x[keep]) < 1e-14, w[keep], 2 * w[keep])
    return np.arccos(x), w


def dvr_matrix(frames: Frames, theta_rad: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """(ntheta, nphi, nrot): psi at the node (theta_t, phi_k) times sqrt(w_t 2 pi / nphi) per FBR state (J, mJ):
    sqrt(w_t / nphi) Pbar_J^mJ(cos theta_t) exp(i mJ phi_k)."""
    phi_rad = frames.extras.angular.phi_rad
    columns = []
    for j, m in frames.states:
        legendre = normalised_legendre(int(j), int(m), np.cos(theta_rad)) * np.sqrt(weights / len(phi_rad))
        columns.append(legendre[:, None] * np.exp(1j * m * phi_rad)[None, :])
    return np.stack(columns, axis=2)


def node_probabilities(psi: np.ndarray, matrix: np.ndarray, volume_bohr4: float) -> np.ndarray:
    """(nz, ntheta, nphi): |psi at the nodes|^2 summed over r, x, y, x the cell volume."""
    return (np.abs(np.einsum("zrxys,tks->zrxytk", psi, matrix)) ** 2).sum(axis=(1, 2, 3)) * volume_bohr4


@pytest.mark.parametrize("run", CHECKPOINTS)
def test_angular_equals_numpy(fixture_runs: Path, run: str) -> None:
    """The manifest's DVR grid is the Gauss-Legendre grid (1e-14), phi_k = 2 pi k / nphi; the last record equals the
    node probabilities of psi_full (1e-12 of its largest). Negative (mid-run, where psi_scat and s2s chi overlap):
    their node probabilities added without the cross terms fail at more than MANY nodes."""
    directory = fixture_runs / run
    frames, info = read_frames([directory]), read_run_info(directory)
    nodes = frames.extras.angular
    theta_rad, weights = angular_nodes(int(frames.states[:, 0].max()), nodes.homonuclear)
    assert np.max(np.abs(theta_rad - nodes.theta_rad)) <= 1e-14
    assert np.max(np.abs(weights - nodes.theta_weights)) <= 1e-14
    assert np.max(np.abs(nodes.phi_rad - 2 * np.pi * np.arange(len(nodes.phi_rad)) / len(nodes.phi_rad))) <= 1e-14
    matrix = dvr_matrix(frames, theta_rad, weights)
    s2s = s2s_of(info, frames.nx, frames.ny)
    psi, _ = psi_full_of(directory, frames, s2s)
    expected = node_probabilities(psi, matrix, cell_volume_bohr4(info))
    got = last_record(frames, "prob_angle_z")
    assert relative_misses(got, expected, got.max(), NUMPY_TOLERANCE) == 0
    if run.endswith("split/first"):
        scat, _ = psi_full_of(directory, frames, s2s, with_chi=False)
        apart = node_probabilities(scat, matrix, cell_volume_bohr4(info))
        apart += node_probabilities(psi - scat, matrix, cell_volume_bohr4(info))
        assert relative_misses(got, apart, got.max(), NUMPY_TOLERANCE) > MANY


@pytest.mark.parametrize("run", CHECKPOINTS)
def test_diffraction_equals_numpy(fixture_runs: Path, run: str) -> None:
    """The last record equals |FFT over (x, y) of exp(-i k_par.R) psi_full|^2 summed over r and the states, x dZ dr
    dx dy sin(gamma) / (nX nY) (1e-12 of its largest). Negative on the oblique runs: without the shift the
    distribution moves (k_par is no multiple of the reciprocal steps) and fails at more than MANY points."""
    directory = fixture_runs / run
    frames, info = read_frames([directory]), read_run_info(directory)
    s2s = s2s_of(info, frames.nx, frames.ny)
    psi, _ = psi_full_of(directory, frames, s2s)
    shift = np.conj(s2s) / info.data["grids"]["unit_cell_norm"]  # exp(-i k_par . R)
    scale = cell_volume_bohr4(info) / (frames.nx * frames.ny)
    expected = (np.abs(np.fft.fft2(psi * shift[None, None, :, :, None], axes=(2, 3))) ** 2).sum(axis=(1, 4)) * scale
    got = last_record(frames, "pop_nm_z")
    assert relative_misses(got, expected, got.max(), NUMPY_TOLERANCE) == 0
    if info.data["wave_packet"]["k_parallel"] > 0:
        unshifted = (np.abs(np.fft.fft2(psi, axes=(2, 3))) ** 2).sum(axis=(1, 4)) * scale
        assert relative_misses(got, unshifted, got.max(), NUMPY_TOLERANCE) > MANY


def test_hard_wall_stays_in_the_zero_order(fixture_runs: Path) -> None:
    """The hard wall couples no (x, y) motion: in every record the orders other than (0, 0) hold at most 1e-14 of
    the record's probability. Negative: the same check on the oblique Dai-Light run, which diffracts, fails."""
    for run, stays in (("hardwall_tiny", True), ("dailight_oblique", False)):
        pop = np.asarray(read_frames([fixture_runs / run]).series.arrays["pop_nm_z"]).sum(axis=1)  # (records, nx, ny)
        total = pop.sum(axis=(1, 2))
        assert (np.max((total - pop[:, 0, 0]) / total) <= 1e-14) == stays, run


def vibrational_states(directory: Path, j_values: np.ndarray, nr: int) -> np.ndarray:
    """(nJ, nr, nr): phi_vJ(r) of gas_phase_states.dat (bohr^-1/2), [J index, v, r point]."""
    states = read_gas_phase_states(directory / "gas_phase_states.dat")
    table = np.zeros((len(j_values), nr, nr))
    for v, j, wave in zip(states.v, states.j, states.wave_functions, strict=True):
        if j in j_values:
            table[int(np.searchsorted(j_values, j)), v] = wave
    return table


def vibration_populations(psi: np.ndarray, frames: Frames, info: RunInfo, phi: np.ndarray) -> np.ndarray:
    """(nz, nr, nJ): |sum_r phi_vJ(r) psi_s(r) dr|^2 summed over x, y and the states s of each J, x dZ dx dy
    sin(gamma); phi: (nJ, v, r)."""
    g = info.data["grids"]
    j_values = frames.extras.vibration_j_values
    populations = np.zeros((frames.nz_scattering, len(frames.r_bohr), len(j_values)))
    for s, (j, _) in enumerate(frames.states):
        j_index = int(np.searchsorted(j_values, j))
        overlaps = np.einsum("vr,zrxy->zvxy", phi[j_index], psi[..., s]) * g["dr"]
        populations[:, :, j_index] += (np.abs(overlaps) ** 2).sum(axis=(2, 3)) * g["dz"] * g["xy_volume_element"]
    return populations


@pytest.mark.parametrize("run", CHECKPOINTS)
def test_vibration_equals_numpy(fixture_runs: Path, run: str) -> None:
    """The last record equals the vibrational populations of psi_full with phi_vJ from gas_phase_states.dat (1e-12
    of its largest). Negative on the Dai-Light runs, where J = 2 is populated: phi_v0 (J = 0) used for every J fails
    at more than MANY points."""
    directory = fixture_runs / run
    frames, info = read_frames([directory]), read_run_info(directory)
    psi, _ = psi_full_of(directory, frames, s2s_of(info, frames.nx, frames.ny))
    phi = vibrational_states(directory, frames.extras.vibration_j_values, len(frames.r_bohr))
    got = last_record(frames, "pop_vj_z")
    assert relative_misses(got, vibration_populations(psi, frames, info, phi), got.max(), NUMPY_TOLERANCE) == 0
    if run.startswith("dailight"):
        wrong = vibration_populations(psi, frames, info, np.repeat(phi[:1], len(phi), axis=0))
        assert relative_misses(got, wrong, got.max(), NUMPY_TOLERANCE) > MANY


@pytest.mark.parametrize("run", CHECKPOINTS)
def test_momentum_equals_numpy(fixture_runs: Path, run: str) -> None:
    """The last record: prob_kr = |FFT over r of psi_full|^2 summed over Z, x, y and the states x the cell volume /
    nR; prob_kz_spec = |FFT over Z of chi|^2 summed over r x dZ dr / nZSp (1e-12 of the largest). Negative: prob_kr
    of psi_scat alone (chi left out) fails at more than MANY points."""
    directory = fixture_runs / run
    frames, info = read_frames([directory]), read_run_info(directory)
    s2s = s2s_of(info, frames.nx, frames.ny)
    psi, chi = psi_full_of(directory, frames, s2s)
    g = info.data["grids"]
    per_k_r = cell_volume_bohr4(info) / len(frames.r_bohr)
    expected_r = (np.abs(np.fft.fft(psi, axis=1)) ** 2).sum(axis=(0, 2, 3, 4)) * per_k_r
    expected_z = (np.abs(np.fft.fft(chi, axis=0)) ** 2).sum(axis=1) * g["dz"] * g["dr"] / len(frames.z_bohr)
    for name, expected in (("prob_kr", expected_r), ("prob_kz_spec", expected_z)):
        got = last_record(frames, name)
        assert relative_misses(got, expected, got.max(), NUMPY_TOLERANCE) == 0, name
    scat, _ = psi_full_of(directory, frames, s2s, with_chi=False)
    wrong = (np.abs(np.fft.fft(scat, axis=1)) ** 2).sum(axis=(0, 2, 3, 4)) * per_k_r
    got = last_record(frames, "prob_kr")
    assert relative_misses(got, wrong, got.max(), NUMPY_TOLERANCE) > MANY


def p2_elements(states: np.ndarray) -> np.ndarray:
    """(nrot, nrot): <J m| P2(cos theta) |J' m'> by a QUADRATURE_POINTS-point Gauss-Legendre quadrature of the
    normalised associated Legendre functions (zero for m != m')."""
    x, w = np.polynomial.legendre.leggauss(QUADRATURE_POINTS)
    p2 = (3 * x**2 - 1) / 2
    elements = np.zeros((len(states), len(states)))
    for a, (j, m) in enumerate(states):
        for b, (jp, mp) in enumerate(states):
            if m == mp:
                elements[a, b] = np.sum(w * p2 * normalised_legendre(j, m, x) * normalised_legendre(jp, m, x))
    return elements


def moment(psi: np.ndarray, operator: np.ndarray, volume_bohr4: float) -> np.ndarray:
    """(nz,): sum over r, x, y of psi^+ A psi x the cell volume, for a matrix A over the rotational states."""
    return np.einsum("zrxys,st,zrxyt->z", np.conj(psi), operator, psi).real * volume_bohr4


@pytest.mark.parametrize("run", CHECKPOINTS)
def test_rotation_moments_equal_numpy(fixture_runs: Path, run: str) -> None:
    """The last record: p2_numerator_z = the moment of P2 (elements by quadrature) and rot_numerator_z = sum of
    |psi_s|^2 J(J+1) / (2 mu r^2) x the cell volume, to 1e-12 of the largest probability at a Z point (times the
    largest J(J+1) / (2 mu r^2) for the energy). Negatives where the packet is (more than 1e-6 of that probability):
    P2 replaced by cos^2, and J(J+1) / (mu r^2) (the 2 dropped), fail at every such Z point."""
    directory = fixture_runs / run
    frames, info = read_frames([directory]), read_run_info(directory)
    psi, _ = psi_full_of(directory, frames, s2s_of(info, frames.nx, frames.ny))
    volume = cell_volume_bohr4(info)
    elements = p2_elements(frames.states)
    j = frames.states[:, 0]
    mu_me = info.data["molecule"]["reduced_mass"]
    energy_hartree = j * (j + 1) / (2 * mu_me * frames.r_bohr[:, None] ** 2)  # (nr, nrot)
    rot = np.einsum("zrxys,rs->z", np.abs(psi) ** 2, energy_hartree) * volume
    probability = last_record(frames, "pop_jmj_z").sum(axis=1)
    scale = float(probability.max())
    energy_scale = scale * float(energy_hartree.max())
    got_p2, got_rot = last_record(frames, "p2_numerator_z"), last_record(frames, "rot_numerator_z")
    assert relative_misses(got_p2, moment(psi, elements, volume), scale, NUMPY_TOLERANCE) == 0
    assert relative_misses(got_rot, rot, energy_scale, NUMPY_TOLERANCE) == 0
    occupied = np.count_nonzero(probability > 1e-6 * scale)
    cos2 = (2 * elements + np.eye(len(frames.states))) / 3  # <cos^2> = (2 P2 + 1) / 3
    assert relative_misses(got_p2, moment(psi, cos2, volume), scale, NUMPY_TOLERANCE) >= occupied
    if j.max() > 0 and run.startswith("dailight"):  # the hard wall stays in J = 0: no rotational energy to scale
        assert relative_misses(got_rot, 2 * rot, energy_scale, NUMPY_TOLERANCE) > MANY


@pytest.mark.parametrize("run", RUNS)
def test_p2_against_the_dvr_quadrature(fixture_runs: Path, run: str) -> None:
    """Every record and Z: <P2> times the probability (analytic elements) equals the DVR quadrature sum over the
    nodes of P2(cos theta_t) prob_angle_z within 1.5 P(J = JMax) at that Z (the quadrature is exact unless both
    states are at JMax; -1/2 <= P2 <= 1, so the two differ by at most 1.5 there) plus 1e-13 of the record's
    probability. The hard wall keeps J = JMax empty: there the two agree to rounding. Negative: cos^2 at the nodes
    instead of P2 fails."""
    frames = read_frames([fixture_runs / run])
    arrays = frames.series.arrays
    cos_theta = np.cos(frames.extras.angular.theta_rad)
    pop = np.asarray(arrays["pop_jmj_z"])
    at_top = pop[:, :, frames.states[:, 0] == frames.states[:, 0].max()].sum(axis=2)
    bound = 1.5 * at_top + SUM_TOLERANCE * pop.sum(axis=(1, 2))[:, None]
    for nodes, right in (((3 * cos_theta**2 - 1) / 2, True), (cos_theta**2, False)):
        quadrature = np.einsum("nztk,t->nz", np.asarray(arrays["prob_angle_z"]), nodes)
        assert np.all(np.abs(np.asarray(arrays["p2_numerator_z"]) - quadrature) <= bound) == right


def strang_r_populations(directory: Path, frames: Frames, info: RunInfo, kinetic_share: float = 0.5) -> np.ndarray:
    """(records, nr): |<phi_v0 | g(t)>|^2 / sum over v, for the r motion alone at the frame times: g starts as
    phi_00 (gas_phase_states.dat) and takes time steps exp(-i T dt/2) [ROpt factor] exp(-i V_asy dt) exp(-i T dt/2),
    T = k^2 / (2 mu) on numpy's FFT grid, V_asy from asymptotic_potential.dat, the ROpt profile from run_info.json
    (kinetic_share: the kinetic steps' share of dt, 0.5 but in the negative test)."""
    r_bohr = frames.r_bohr
    dt_au, mu_me = info.data["times"]["dt_step"], info.data["molecule"]["reduced_mass"]
    phi = vibrational_states(directory, np.array([0]), len(r_bohr))[0]  # (v, r) of J = 0
    k_per_bohr = 2 * np.pi * np.fft.fftfreq(len(r_bohr), r_bohr[1] - r_bohr[0])
    kinetic_part = np.exp(-1j * k_per_bohr**2 / (2 * mu_me) * dt_au * kinetic_share)
    v_hartree = read_asymptotic_potential(directory / "asymptotic_potential.dat").v_hartree
    step = absorber_factors(info, absorbers_on(info, "r", "specular"), r_bohr) * np.exp(-1j * v_hartree * dt_au)
    wave = phi[0].astype(complex)
    times_au = frames.series.times_au - info.data["times"]["t_start"]
    populations, done = [], 0
    for steps in np.round(times_au / dt_au).astype(int):
        for _ in range(steps - done):
            wave = np.fft.ifft(kinetic_part * np.fft.fft(step * np.fft.ifft(kinetic_part * np.fft.fft(wave))))
        done = steps
        c2 = np.abs(phi @ wave) ** 2
        populations.append(c2 / c2.sum())
    return np.array(populations)


def test_hard_wall_vibration_follows_the_r_motion_alone(fixture_runs: Path) -> None:
    """The hard wall is separable in Z and r, so psi_full = f(Z) g(r) and pop_vj_z(Z) / its sum over v is the same at
    every Z: |c_v(t)|^2 / sum |c|^2 of the r motion alone, propagated here in numpy with the run's grid, time step,
    V_asy and r absorber (to 1e-12, where the Z point holds more than 1e-6 of the record's largest; the split
    operator moves ~1e-8 out of v = 0). Its premises are asserted: no VMax clamp (it would couple Z and r in the
    wall), J = 0 (no centrifugal term), and T dt < pi at every r momentum (the program's half step is the principal
    square root of exp(-i T dt), numpy's exp(-i T dt / 2): equal only then). Negative: full kinetic steps fail."""
    directory = fixture_runs / "hardwall_tiny"
    frames, info = read_frames([directory]), read_run_info(directory)
    assert info.data["potential"]["impose_v_max"] is False and info.data["initial_state"]["j"] == 0
    dr_bohr = frames.r_bohr[1] - frames.r_bohr[0]
    assert (np.pi / dr_bohr) ** 2 / (2 * info.data["molecule"]["reduced_mass"]) * info.data["times"]["dt_step"] < np.pi
    pop = np.asarray(frames.series.arrays["pop_vj_z"])[:, :, :, 0]  # (records, nz, v) of J = 0
    per_z = pop.sum(axis=2)
    shown = per_z > 1e-6 * per_z.max(axis=1)[:, None]
    ratios = pop / np.where(shown, per_z, 1.0)[:, :, None]
    expected = strang_r_populations(directory, frames, info)
    assert np.max(np.where(shown[:, :, None], np.abs(ratios - expected[:, None, :]), 0.0)) <= 1e-12
    assert np.max(expected[:, 1:].sum(axis=1)) > 1e-9  # the split operator does move some out of v = 0
    wrong = strang_r_populations(directory, frames, info, kinetic_share=1.0)
    assert np.max(np.where(shown[-1][:, None], np.abs(ratios[-1] - wrong[-1][None, :]), 0.0)) > 1e-12


@pytest.mark.parametrize("run", RUNS)
def test_split_run_joins_to_the_unsplit_extras(fixture_runs: Path, run: str) -> None:
    """The two segments of a split run, joined, hold the unsplit run's extras bitwise."""
    split = fixture_runs / f"{run}_split"
    joined = read_frames([split / "first", split / "second"])
    whole = read_frames([fixture_runs / run])
    for names in EXTRA_DATASETS.values():
        for name in names:
            assert np.array_equal(np.asarray(joined.series.arrays[name]), np.asarray(whole.series.arrays[name])), name


def test_reader_refuses_extras_it_cannot_index(fixture_runs: Path, tmp_path: Path) -> None:
    """frames/ whose manifest names an unknown extra, lists an extra without its dataset, the angular extra without
    its DVR grid, the vibration extra without its J values, or a DVR grid that does not give prob_angle_z's shape is
    a ReadError naming what is wrong; a manifest without "extras" and datasets (frames/ of an older program) reads as
    no extras."""
    fixture = fixture_runs / "dailight_oblique"
    manifest = run_with_frames_manifest(fixture, tmp_path / "same", {}) / "frames" / "manifest.json"
    extras = read_frames([manifest.parent.parent]).series.manifest["extras"]
    datasets = read_frames([manifest.parent.parent]).series.manifest["datasets"]
    fewer_nodes = {**extras["angular"], "theta_rad": extras["angular"]["theta_rad"][:-1]}
    cases = {
        "unknown FrameExtras name 'spin'": {"extras": {**extras, "names": [*extras["names"], "spin"]}},
        "FrameExtras momentum without prob_kr": {"datasets": {k: v for k, v in datasets.items() if k != "prob_kr"}},
        "the angular extra without its DVR grid": {"extras": {k: v for k, v in extras.items() if k != "angular"}},
        "the vibration extra without its J values": {
            "extras": {k: v for k, v in extras.items() if k != "vibration_j_values"}
        },
        "prob_angle_z has records of shape": {"extras": {**extras, "angular": fewer_nodes}},
    }
    for k, (message, change) in enumerate(cases.items()):
        with pytest.raises(ReadError, match=message):
            read_frames([run_with_frames_manifest(fixture, tmp_path / f"case{k}", change)])
    older = {"extras": None, "datasets": {k: v for k, v in datasets.items() if k in DATASETS}}
    assert read_frames([run_with_frames_manifest(fixture, tmp_path / "older", older)]).extras.names == ()
