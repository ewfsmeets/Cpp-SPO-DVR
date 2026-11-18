"""The rovibrational states of the free molecule and the open channels.

Shadows: src/gasphase/fourier_grid_hamiltonian.cpp (FourierGridHamiltonian, fix_sign, normalise),
src/gasphase/gas_phase.cpp (solve_rovibrational_states, count_open_vibrational_levels,
highest_open_rotational_level, compute_gas_phase)

The program diagonalises the Fourier grid Hamiltonian (FGH) of V(r) + j(j+1)/(2 mu r^2) for j = 0 .. JMax,
with V(r) the potential's gas-phase curve on the input's r grid. The designer evaluates no potential: it
gets a curve on its own grid (GasPhaseCurve, e.g. dz_gas_phase.py) and solves on that grid. The energies
come from numpy's LAPACK, the program's from oneMKL or the system LAPACK: equal to rounding.
"""

from dataclasses import dataclass

import numpy as np

from program.angular_grid import AngularGrid
from program.config import DerivedParameters, InputParameters
from program.fft_sizes import nearest_fft_size
from program.unit_conversion import PI, ev_to_hartree

FORTRAN_MACHINE_EPSILON = 2.0**-54  # Constants_EPS of the Fortran (fortran_machine_epsilon)


@dataclass(frozen=True)
class GasPhaseCurve:
    """V(r) of the free molecule on an equidistant r grid (atomic units)."""

    r_start: float  # bohr
    dr: float  # bohr
    potential: tuple[float, ...]  # hartree, at r_start + i * dr
    system_id: str  # the SystemID whose curve it is
    source: str  # where the curve comes from, for the reports


@dataclass(frozen=True)
class FourierGridHamiltonian:
    """The FGH grid and its kinetic energy matrix (FourierGridHamiltonian)."""

    r: np.ndarray  # bohr: lower_bound + l * dx, l = 0 .. n-1
    dx: float  # bohr
    kinetic: np.ndarray  # hartree, n x n


@dataclass(frozen=True)
class GasPhase:
    """The rovibrational states and the open channels (GasPhase)."""

    r: np.ndarray  # bohr: the FGH grid
    energies: np.ndarray  # hartree, [v, j]
    wave_functions: np.ndarray  # [v, j, ir], sum psi^2 dr = 1
    v_max_open: int  # number of open vibrational levels. Fortran: ivmaxAnalysis
    j_max_open: int  # highest open j. Fortran: jmaxAnalysis
    n_v_analysis: int  # number of vibrational levels analysed (VAna, or v_max_open). Fortran: nviba
    j_max_analysis: int  # highest j analysed (JAna, or j_max_open). Fortran: jana
    warnings: tuple[str, ...]
    curve: GasPhaseCurve  # the V(r) the states were computed for


def compute_gas_phase(p: InputParameters, d: DerivedParameters, angular: AngularGrid, curve: GasPhaseCurve) -> GasPhase:
    """The states for the input's masses and JMax on the curve's grid, and the channels open up to EMaxWP plus
    the initial state's internal and parallel energy (compute_gas_phase). VIn must be a level of the curve's
    grid (VIn < its number of points)."""
    n = len(curve.potential)
    if nearest_fft_size(n) != n:  # the program's nR is a valid FFT size; the FGH needs one
        raise ValueError(f"the gas-phase curve has {n} points, not a valid FFT size")
    fgh = fourier_grid_hamiltonian(n, curve.r_start, curve.r_start + n * curve.dr, d.reduced_mass)
    energies, wave_functions = solve_rovibrational_states(fgh, curve.potential, p.j_max, d.reduced_mass)

    e00 = energies[0, 0]
    e_initial = energies[p.v_in, p.j_in]
    available_energy = ev_to_hartree(p.e_max_wave_packet) + (e_initial - e00) + ev_to_hartree(p.e_parallel)
    warnings: list[str] = []
    v_max_open = count_open_vibrational_levels(energies, available_energy, warnings)
    j_max_open = highest_open_rotational_level(energies, angular.delta_j, available_energy, warnings)
    # VAna < 1 (or absent) and absent JAna mean: determined by the program.
    n_v_analysis = p.v_analysis if p.v_analysis is not None and p.v_analysis >= 1 else v_max_open
    j_max_analysis = p.j_analysis if p.j_analysis is not None and p.j_analysis >= 0 else j_max_open
    warnings += analysis_warnings(p, d, angular, energies, e_initial, n_v_analysis, j_max_analysis)
    return GasPhase(
        fgh.r,
        energies,
        wave_functions,
        v_max_open,
        j_max_open,
        n_v_analysis,
        j_max_analysis,
        tuple(warnings),
        curve,
    )


def analysis_warnings(
    p: InputParameters,
    d: DerivedParameters,
    angular: AngularGrid,
    energies: np.ndarray,
    e_initial: float,
    n_v_analysis: int,
    j_max_analysis: int,
) -> list[str]:
    """Warnings if a channel within the packet's 95% energy interval is not analysed (JAna or VAna too small);
    they do not change the analysis (compute_gas_phase)."""
    warnings = []
    e_check = e_initial + ev_to_hartree(d.confidence_e_max)
    next_j = j_max_analysis + angular.delta_j
    j_limit = p.j_max - 1 if p.homonuclear else p.j_max
    if j_max_analysis < j_limit and e_check > energies[0, next_j]:
        warnings.append("Not all open j channels analysed!")
    if n_v_analysis < energies.shape[0] and e_check > energies[n_v_analysis, 0]:
        warnings.append("Not all open v channels analysed!")
    return warnings


def fourier_grid_hamiltonian(
    n_approx: int, lower_bound: float, upper_bound: float, reduced_mass: float
) -> FourierGridHamiltonian:
    """The FGH on [lower_bound, upper_bound): its grid and kinetic energy matrix (FourierGridHamiltonian).

    T applies the kinetic energy (2 pi l / (n dx))^2 / (2 mu) in momentum space to each column of the identity.
    The program's FFTs have the Fortran sign convention; for a real kinetic energy T is the same.
    """
    n = nearest_fft_size(n_approx, 1)
    dx = (upper_bound - lower_bound) / n
    r = np.array([l * dx + lower_bound for l in range(n)])
    kinetic_diagonal = np.array([2.0 / reduced_mass * (l * PI / n / dx) ** 2.0 for l in range(n // 2 + 1)])
    kinetic = np.fft.irfft(kinetic_diagonal[:, np.newaxis] * np.fft.rfft(np.eye(n), axis=0), n=n, axis=0)
    return FourierGridHamiltonian(r, dx, kinetic)


def solve_rovibrational_states(
    fgh: FourierGridHamiltonian, potential: tuple[float, ...], j_max: int, reduced_mass: float
) -> tuple[np.ndarray, np.ndarray]:
    """Energies [v, j] and wave functions [v, j, ir] for j = 0 .. JMax: V(r) plus the centrifugal term
    j(j+1)/(2 mu r^2) (solve_rovibrational_states)."""
    n = len(fgh.r)
    energies = np.zeros((n, j_max + 1))
    wave_functions = np.zeros((n, j_max + 1, n))
    for j in range(j_max + 1):
        total_potential = list(potential)
        if j > 0:
            for i in range(n):
                r = fgh.r[i] if fgh.r[i] != 0.0 else FORTRAN_MACHINE_EPSILON  # prevent division by zero
                total_potential[i] = potential[i] + (j * (j + 1)) / 2.0 / reduced_mass / (r * r)
        values, vectors = np.linalg.eigh(fgh.kinetic + np.diag(total_potential))  # ascending energies
        for v in range(n):
            energies[v, j] = values[v]
            wave_functions[v, j] = normalise(fix_sign(vectors[:, v]), fgh.dx)
    return energies, wave_functions


def fix_sign(psi: np.ndarray) -> np.ndarray:
    """psi with its first local maximum of |psi|, seen from the last grid point, positive (fix_sign)."""
    i = len(psi) - 1
    while i > 0 and abs(psi[i - 1]) >= abs(psi[i]):
        i -= 1
    return -1.0 * psi if psi[i] < 0.0 else psi


def normalise(psi: np.ndarray, dx: float) -> np.ndarray:
    """psi scaled to sum psi^2 dx = 1 (normalise)."""
    norm = 0.0
    for y in psi:
        norm = norm + y * y
    return psi * (1.0 / np.sqrt(norm * dx))


def count_open_vibrational_levels(energies: np.ndarray, available_energy: float, warnings: list[str]) -> int:
    """Open vibrational levels: count v (j = 0) until E(v, 0) - E(0, 0) exceeds the available energy
    (count_open_vibrational_levels; 1-based like the Fortran loop)."""
    e00 = energies[0, 0]
    n_v = energies.shape[0]
    iv = 1
    while energies[iv - 1, 0] - e00 <= available_energy:
        iv = iv + 1
        if iv == n_v:
            warnings.append("Vibrational basis too small to describe all open channels!")
            break
    return iv - 1


def highest_open_rotational_level(
    energies: np.ndarray, delta_j: int, available_energy: float, warnings: list[str]
) -> int:
    """The highest open j (v = 0), in steps of delta_j (highest_open_rotational_level)."""
    e00 = energies[0, 0]
    j_max = energies.shape[1] - 1
    j = 0
    while energies[0, j] - e00 <= available_energy:
        j = j + delta_j
        if j > j_max:
            warnings.append("Rotational basis too small to describe all open channels!")
            break
    return j - delta_j
