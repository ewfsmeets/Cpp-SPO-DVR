"""Reflection R and transmission T of a Z absorber in one dimension, and the amplitude it returns to the
analysis line (inputs/README.md sections 7 and 8; the verdict's "absorber reflection" finding).

A plane wave of energy E meets V = -i A (3/2) ((q - start)/L)^2 on [start, start + L]: R is what comes back, T
what passes the end. On the periodic Z grid the transmitted part wraps round to ZStart, is reflected by the
repulsive wall and crosses the absorber again, so it returns with amplitude about T: the returned amplitude is
a = sqrt(R) + T. The method is the wave-packet designer's (tools/wave_packet_designer/design/transmission.py,
reflection_transmission: a transfer matrix through thin slices of constant V, integrated backwards from the
transmitted wave), vectorised over the energies; tests/test_absorber_reflection.py compares the two.
"""

import numpy as np

SLICE_BOHR = 0.002  # the designer's SLICE (design/transmission.py); tests/test_absorber_reflection.py compares with it
RESCALE_ABOVE = 1e100  # psi and psi' are divided by their largest modulus above this (no inf/inf in a long absorber)


def reflection_transmission(
    energies_hartree: np.ndarray, *, length_bohr: float, strength_hartree: float, mass_me: float
) -> tuple[np.ndarray, np.ndarray]:
    """R and T at the kinetic energies (hartree, > 0) for an absorber of the given length and strength A. Each
    energy is computed on its own: the result of a row does not depend on the other rows."""
    energies = np.asarray(energies_hartree, dtype=float)
    k_free_per_bohr = np.sqrt(2.0 * mass_me * energies)
    n_slices = max(1, round(length_bohr / SLICE_BOHR))
    h_bohr = length_bohr / n_slices
    # Beyond the end only the transmitted wave exp(i k q), amplitude 1; go back slice by slice to the start.
    psi = np.ones(energies.shape, dtype=complex)
    dpsi_per_bohr = 1j * k_free_per_bohr
    log_scale = np.zeros(energies.shape)  # log of what psi and psi' were divided by
    for i in reversed(range(n_slices)):
        x = (i + 0.5) / n_slices  # slice centre, relative to the absorber
        potential_hartree = -1j * strength_hartree * 1.5 * x * x
        k_per_bohr = np.sqrt(2.0 * mass_me * (energies - potential_hartree))
        c, s = np.cos(k_per_bohr * h_bohr), np.sin(k_per_bohr * h_bohr)
        psi, dpsi_per_bohr = c * psi - s / k_per_bohr * dpsi_per_bohr, k_per_bohr * s * psi + c * dpsi_per_bohr
        psi, dpsi_per_bohr, log_scale = rescaled(psi, dpsi_per_bohr, log_scale)
    # Split the wave at the start into incoming exp(i k q) and reflected exp(-i k q) parts.
    incoming = (psi + dpsi_per_bohr / (1j * k_free_per_bohr)) / 2.0
    reflected = (psi - dpsi_per_bohr / (1j * k_free_per_bohr)) / 2.0
    transmission = np.exp(-2.0 * log_scale) / np.abs(incoming) ** 2  # the factor is 1 without a rescaling
    return np.abs(reflected / incoming) ** 2, transmission


def rescaled(
    psi: np.ndarray, dpsi_per_bohr: np.ndarray, log_scale: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """psi and psi' divided by their largest modulus where it exceeds RESCALE_ABOVE (the ratios R needs are
    unchanged; log_scale records the factor for T); unchanged elsewhere."""
    largest = np.maximum(np.abs(psi), np.abs(dpsi_per_bohr))
    large = largest > RESCALE_ABOVE
    if not large.any():
        return psi, dpsi_per_bohr, log_scale
    divisor = np.where(large, largest, 1.0)
    return psi / divisor, dpsi_per_bohr / divisor, log_scale + np.log(divisor)


def returned_amplitude(
    energies_hartree: np.ndarray, *, length_bohr: float, strength_hartree: float, mass_me: float
) -> tuple[np.ndarray, np.ndarray]:
    """Per energy (hartree, > 0): R and the amplitude a = sqrt(R) + T the absorber returns to the analysis line.
    Computed once per distinct energy (the rows of an S matrix repeat channel energies)."""
    distinct, row_to_distinct = np.unique(np.asarray(energies_hartree, dtype=float), return_inverse=True)
    reflection, transmission = reflection_transmission(
        distinct, length_bohr=length_bohr, strength_hartree=strength_hartree, mass_me=mass_me
    )
    amplitude = np.sqrt(reflection) + transmission
    return reflection[row_to_distinct], amplitude[row_to_distinct]
