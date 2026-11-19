"""The S matrix from the Balint-Kurti coefficients, as the program computes it: the time sums, and the
prefactor and phase of each open channel.

S_f(E) = sum_t c_f(t) exp(i E_tot t) dt_ana * (alpha / 2 pi)^(1/4) / M * exp(i arg_f) * sqrt(k_f k0 / weight)
with arg_f = -k_f Zana + (kz - k0) Z0 - pi, weight = exp(-(kz - k0)^2 / (2 alpha)), k0 the normal momentum
of the collision energy and k_f that of channel f. The initial packet is (2 alpha / pi)^(1/4)
exp(-alpha (Z - Z0)^2 - i kz Z): alpha (bohr^-2) sets its width, kz (bohr^-1) its mean normal momentum.

Shadows: src/analysis/s_matrix.cpp (SMatrixAccumulator::add, SMatrixAccumulator::result_at)
"""

import numpy as np


def time_sums(
    times_au: np.ndarray, coefficients: np.ndarray, total_energy_hartree: np.ndarray, analysis_time_step_au: float
) -> np.ndarray:
    """sum_t c(t) exp(i E_tot t) dt_ana for every energy and channel.

    times_au: (n_times,); coefficients: (n_times, n_channels) complex, in bohr^-1/2; total_energy_hartree:
    (n_energies,). Returns: (n_energies, n_channels) complex. A matrix product with dt_ana applied once at
    the end, where the C++ adds c dt_ana exp(i E t) record by record: equal to rounding (the S matrices agree
    to 1e-12 of the largest |S|, tests/test_runs_s_matrix.py).
    """
    phases = np.exp(1j * np.outer(total_energy_hartree, times_au))  # (n_energies, n_times)
    return (phases @ coefficients) * analysis_time_step_au


def normal_momentum(
    collision_energy_ev: np.ndarray,
    *,
    total_mass_me: float,
    hartree_to_ev: float,
    kx_per_bohr: float,
    ky_per_bohr: float,
    skew_angle_rad: float,
) -> np.ndarray:
    """k0 (bohr^-1): the momentum normal to the surface at the collision energies.

    k_parallel^2 = kx^2 - 2 cos(skew) kx ky + ky^2, with kx, ky the components along the skewed reciprocal
    axes (run_info.json wave_packet.kx, .ky).
    """
    kx, ky = kx_per_bohr, ky_per_bohr
    k_parallel_sq = kx * kx - (2.0 * np.cos(skew_angle_rad) * kx * ky) + ky * ky
    k_total = np.sqrt(2.0 * total_mass_me * collision_energy_ev / hartree_to_ev)
    return np.sqrt(k_total * k_total - k_parallel_sq)


def packet_weight(k0_per_bohr: np.ndarray, kz_per_bohr: float, alpha_per_bohr2: float) -> np.ndarray:
    """exp(-(kz - k0)^2 / (2 alpha)): the weight of the normal momentum k0 in the initial packet, 1 at its
    centre kz (dimensionless)."""
    dk = kz_per_bohr - k0_per_bohr
    return np.exp(-(dk * dk) / (2.0 * alpha_per_bohr2))


def s_prefactor(alpha_per_bohr2: float, total_mass_me: float) -> float:
    """(alpha / 2 pi)^(1/4) / M, the prefactor of every S-matrix element (atomic units)."""
    return float(np.power(0.5 * alpha_per_bohr2 / np.pi, 0.25) / total_mass_me)


def channel_factors(
    total_energy_ev: np.ndarray,
    channel_energy_hartree: np.ndarray,
    k0_per_bohr: np.ndarray,
    *,
    kz_per_bohr: float,
    alpha_per_bohr2: float,
    total_mass_me: float,
    z_analysis_bohr: float,
    z0_bohr: float,
    hartree_to_ev: float,
) -> np.ndarray:
    """prefactor * exp(i arg_f) * sqrt(k_f k0 / weight) for every energy and channel, 0 for a closed channel
    (k_f^2 <= 0): S = time sums * these factors.

    total_energy_ev, k0_per_bohr: (n_energies,); channel_energy_hartree: (n_channels,).
    Returns: (n_energies, n_channels) complex, in bohr^1/2 (the time sums are in bohr^-1/2).
    """
    k0 = k0_per_bohr[:, None]
    weight = packet_weight(k0, kz_per_bohr, alpha_per_bohr2)
    kf_sq = 2.0 * total_mass_me * (total_energy_ev[:, None] / hartree_to_ev - channel_energy_hartree[None, :])
    open_channel = kf_sq > 0.0
    kf = np.sqrt(np.where(open_channel, kf_sq, 1.0))  # 1.0: any positive value; closed channels are zeroed below
    arg = -kf * z_analysis_bohr + (kz_per_bohr - k0) * z0_bohr - np.pi
    factors = s_prefactor(alpha_per_bohr2, total_mass_me) * np.exp(1j * arg) * np.sqrt(kf * k0 / weight)
    return np.where(open_channel, factors, 0.0)
