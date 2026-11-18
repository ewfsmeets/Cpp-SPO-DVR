"""The initial wave packet in Z on the specular grid, its momentum distribution, and the TolSp check.

Shadows: src/propagation/initial_wave_packet.cpp (make_initial_wave_packet, last_point_inside, norm_sum)

psi(Z) = (2 alpha / pi)^(1/4) exp(-alpha (Z - Z0)^2 - i kz Z): a Gaussian moving towards the surface. The
program multiplies it by the rovibrational state phi_vj(r); the Z part alone decides the TolSp check, so the
designer computes only that. The norm sums run in the C++ order (a plain loop, not numpy's pairwise sum),
so the TolSp decision is the program's even at the boundary. The momentum distribution uses numpy's FFT
(the program FFTW or oneMKL): equal to rounding, which matters only for values near zero.
"""

import cmath
import math
from dataclasses import dataclass

import numpy as np

from program.config import DerivedParameters, InputParameters
from program.findings import Finding, FindingList
from program.real_compare import rge
from program.spatial_grids import SpatialGrids, first_index
from program.unit_conversion import PI


@dataclass(frozen=True)
class InitialWavePacketInfo:
    """The quantities the program prints in "ROTATIONAL BASIS AND INITIAL WAVE PACKET" (Z part)."""

    norm: float  # sum |psi|^2 dZ over the whole specular grid
    k_space_maximum: float  # max |psi(k)|
    k_space_at_zero: float  # |psi(k = 0)|
    k_space_at_edge: float  # |psi(k)| at the largest momentum of the grid (index n/2)
    i_low: int  # first specular point of the allowed region (the analysis line)
    i_high: int  # last specular point of the allowed region (before the absorber or ZAnaSp)
    norm_inside: float  # sum |psi|^2 dZ over the allowed region


@dataclass(frozen=True)
class InitialWavePacket:
    """psi(Z) on the specular grid and the printed quantities (the Z part of the C++ InitialWavePacket)."""

    psi: np.ndarray  # complex, one value per specular point
    info: InitialWavePacketInfo


def make_initial_wave_packet(
    p: InputParameters, d: DerivedParameters, grids: SpatialGrids
) -> tuple[InitialWavePacket, list[Finding]]:
    """The packet and the TolSp error the program gives when too much norm lies outside the allowed region."""
    checks = FindingList("make_initial_wave_packet")
    z = grids.z_specular
    alpha = d.wave_packet_alpha
    amplitude = (2.0 * alpha / PI) ** 0.25
    psi = []
    for iz in range(z.size):
        dz = z.point(iz) - p.z0
        psi.append(amplitude * cmath.exp(complex(-(alpha * (dz * dz)), -(d.kz * z.point(iz)))))

    # Momentum distribution, scaled like the program's k_spectrum (dZ / sqrt(2 pi)).
    psi_k = np.fft.fft(np.array(psi)) * (p.dz / math.sqrt(2.0 * PI))
    magnitudes = np.abs(psi_k)

    # The allowed region: from the analysis line to the last point before the specular absorber (or ZAnaSp).
    i_low = first_index(z, lambda point: rge(point, grids.z_analysis))
    assert i_low is not None  # Zana is a point of the scattering grid, which the specular grid contains
    i_high = last_point_inside(p, grids)
    info = InitialWavePacketInfo(
        norm=norm_sum(psi, 0, z.size - 1) * p.dz,
        k_space_maximum=float(magnitudes.max()),
        k_space_at_zero=float(magnitudes[0]),
        k_space_at_edge=float(magnitudes[z.size // 2]),
        i_low=i_low,
        i_high=i_high,
        norm_inside=norm_sum(psi, i_low, i_high) * p.dz,
    )
    outside = 1.0 - info.norm_inside
    checks.error_if(
        outside > p.tol_specular,
        "TolSp",
        f"too much norm of the initial wave packet outside the allowed region (Z = {z.point(i_low)} .. "
        f"{z.point(i_high)}): 1 - norm = {outside} > tolsp; move Z0, widen the energy range or "
        f"enlarge the specular grid",
    )
    return InitialWavePacket(np.array(psi), info), checks.findings


def norm_sum(psi: list[complex], first: int, last: int) -> float:
    """sum |psi|^2 over first .. last, in order, as the Fortran's conjg(psi)*psi sum (norm_sum)."""
    total = 0.0
    for i in range(first, last + 1):
        total = total + (psi[i].real * psi[i].real + psi[i].imag * psi[i].imag)
    return total


def last_point_inside(p: InputParameters, grids: SpatialGrids) -> int:
    """The last specular point allowed to carry the packet: before ZAnaSp if that is fixed, else before the absorber."""
    fixed_specular_analysis = not p.dynamical_zsp_analysis
    if fixed_specular_analysis and p.z_analysis_specular is not None and p.z_analysis_specular > 0.0:
        first_beyond = first_index(grids.z_specular, lambda point: rge(point, p.z_analysis_specular))
        if first_beyond is not None:
            return first_beyond - 1
    return grids.z_specular_optical.first - 1  # Fortran: nschmi
