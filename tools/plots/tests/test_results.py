"""compute/results.py's packet spectrum on hand-made run_info data: normalised over the energies it holds, peaked at
the packet's mean normal energy, and exp(-2) of the peak per momentum at the ends of the packet's energy range."""

from pathlib import Path

import numpy as np

from compute.results import packet_spectrum_per_ev
from readers.run_info import RunInfo


def packet_info(kz: float, alpha: float) -> RunInfo:
    """A run_info of H2 (M = 3674 me) at normal incidence with the packet (kz, alpha)."""
    data = {
        "wave_packet": {"kx": 0.0, "ky": 0.0, "kz": kz, "alpha": alpha},
        "lattice": {"skew_angle_rad": np.pi / 2},
        "molecule": {"total_mass": 3674.0},
    }
    return RunInfo(Path("."), data, 27.2114, 0.0, 1.0)


def test_spectrum_is_normalised_and_peaks_at_the_mean_energy() -> None:
    """kz = 7 bohr^-1 (E_n = 49 / (2 x 3674) hartree = 0.181 eV), alpha = 0.5: the spectrum integrates to 1 over
    the energies (dE/dk = k / M, so the momentum distribution's integral is 1). Without the Jacobian M / k0 (a slip)
    the integral is about 1.9e-3, so the check catches it; the peak (moved by the Jacobian by about 4 meV) lies
    within 10 meV of 0.181 eV."""
    info = packet_info(7.0, 0.5)
    e = np.linspace(1e-6, 2.0, 400001)
    spectrum = packet_spectrum_per_ev(info, e)
    assert abs(np.trapezoid(spectrum, e) - 1.0) < 1e-6
    peak = 49.0 / (2 * 3674.0) * 27.2114
    assert abs(e[np.argmax(spectrum)] - peak) < 0.01


def test_weight_at_the_ends_of_the_energy_range() -> None:
    """A packet made as the program makes it from EMinWP = 0.1 and EMaxWP = 0.9 eV (wave_packet_params.cpp: kz the
    mean of k_min and k_max, alpha = (kz - k_min)^2 / 4): at both ends the weight per momentum (the spectrum times
    k0 / M) is exp(-2) of its peak, what the results explanation states. A slip alpha = (kz - k_min)^2 / 2 gives
    exp(-1)."""
    mass, to_ev = 3674.0, 27.2114
    k_min, k_max = (np.sqrt(2.0 * mass * e / to_ev) for e in (0.1, 0.9))
    kz = 0.5 * (k_min + k_max)
    info = packet_info(kz, (kz - k_min) ** 2 / 4.0)
    ends = np.array([0.1, 0.9, kz * kz / (2.0 * mass) * to_ev])
    per_momentum = packet_spectrum_per_ev(info, ends) * np.sqrt(2.0 * mass * ends / to_ev)
    assert np.allclose(per_momentum[:2] / per_momentum[2], np.exp(-2.0), rtol=1e-12, atol=0)
