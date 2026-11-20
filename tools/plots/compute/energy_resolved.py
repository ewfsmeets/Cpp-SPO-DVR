"""The stationary scattering state at one collision energy from the run's energy-resolved sums (energy_resolved/):
psi_E = (s + s2s sigma in the initial channel) * (alpha / 2 pi)^(1/4) / M * k0 / sqrt(weight), the S matrix's own
prefactor (src/analysis/s_matrix.cpp) without its channel factors, so that the incident wave has unit modulus per
unit-cell area: psi_E = exp(i (k0 - kz) Z0) (exp(-i k0 Z) phi_v(r) + scattered waves) exp(i K_par R) / sqrt(A), and
A |psi_E|^2 integrated over r is 1 for the incident wave alone (A the unit cell's area). The constant phase
exp(i (k0 - kz) Z0) (the packet's at Z0) is left in: densities do not see it. k0 is the S matrix's (without the
split operator's O(dt^2) shift of the propagated energy, as s_matrix.tsv and the Fortran). Valid only where the
packet's own start and the absorbers do not enter (valid_region).

Units: psi_E in bohr^-3/2 (per sqrt(area) and sqrt(r)); A rho_E per bohr of r.
"""

from dataclasses import dataclass

import numpy as np

from copies.s_matrix import normal_momentum, packet_weight, s_prefactor
from readers.energy_resolved import EnergyResolved
from readers.run_info import RunInfo

# Derived: 6 sigma_Z below Z0 the starting packet's amplitude tail is exp(-9) = 1.2e-4; psi_E is linear in psi, and a
# free-particle model of the time sum gives a 1e-5 error of the incident wave there (2-6e-4 at 5 sigma_Z).
PACKET_WIDTHS = 6.0


def energy_factor(info: RunInfo, collision_ev: float) -> complex:
    """(alpha / 2 pi)^(1/4) / M * k0 / sqrt(weight) at one collision energy (k0 from the packet's parallel momentum,
    weight = exp(-(kz - k0)^2 / 2 alpha), copies/s_matrix.py)."""
    packet, lattice = info.data["wave_packet"], info.data["lattice"]
    mass = info.data["molecule"]["total_mass"]
    k0 = float(
        normal_momentum(
            np.array([collision_ev]),
            total_mass_me=mass,
            hartree_to_ev=info.hartree_to_ev,
            kx_per_bohr=packet["kx"],
            ky_per_bohr=packet["ky"],
            skew_angle_rad=lattice["skew_angle_rad"],
        )[0]
    )
    weight = float(packet_weight(np.array([k0]), packet["kz"], packet["alpha"])[0])
    return s_prefactor(packet["alpha"], mass) * k0 / np.sqrt(weight)


def whole_wave_sums(sums: EnergyResolved, info: RunInfo, site: int, energy: int, partial: bool = False) -> np.ndarray:
    """(nz, nr, n_rot) complex: the sums of the whole wave packet at one site and energy on the scattering grid, the
    scattering wave's plus s2s times the specular wave's in the initial rotational channel (the specular wave is
    the initial channel kept apart; docs/output_files.md "energy_resolved/"). partial: the sums up to 0.8 tEnd."""
    scattering = sums.sums_partial if partial else sums.sums
    specular = sums.specular_sums_partial if partial else sums.specular_sums
    whole = np.array(scattering[site, energy])
    initial = int(info.data["rotational_basis"]["initial_index"])
    nz = len(sums.z_bohr)
    whole[:, :, initial] += sums.sites[site].specular_to_scattering * np.asarray(specular[energy, :nz])
    return whole


def stationary_state(sums: EnergyResolved, info: RunInfo, site: int, energy: int, partial: bool = False) -> np.ndarray:
    """psi_E (nz, nr, n_rot) at one site and energy (bohr^-3/2): whole_wave_sums times energy_factor."""
    return whole_wave_sums(sums, info, site, energy, partial) * energy_factor(info, float(sums.collision_ev[energy]))


def area_density(state: np.ndarray, info: RunInfo) -> np.ndarray:
    """A rho_E(Z, r) = A sum over the rotational states of |psi_E|^2 (nz, nr; per bohr of r), A the unit cell's
    area: 1 integrated over r for the incident wave alone."""
    return info.data["grids"]["unit_cell_area"] * np.sum(np.abs(state) ** 2, axis=2)


@dataclass(frozen=True, kw_only=True)
class ValidRegion:
    """Where psi_E is the stationary state: above the low-Z absorbers' end (of both grids: the specular wave is its
    initial channel), below the Z absorbers' start (both grids) and the starting packet (Z0 - PACKET_WIDTHS sigma_Z),
    and below the r absorber's start."""

    z_low_bohr: float
    z_high_bohr: float
    r_high_bohr: float


def valid_region(info: RunInfo) -> ValidRegion:
    """The region of the scattering grid where the stationary state holds (the absorbers from run_info.json)."""
    absorbers = info.data["absorbers"]["list"]
    z_low = info.data["grids"]["z_start"]
    z_high = info.data["wave_packet"]["z0"] - PACKET_WIDTHS * info.data["wave_packet"]["sigma_z"]
    r_high = np.inf
    for a in absorbers:
        if a["name"] in ("ZOpt", "ZSpOpt"):
            z_high = min(z_high, a["start"])
        elif a["name"] in ("ZOptLow", "ZSpOptLow"):
            z_low = max(z_low, a["end"])
        elif a["name"] == "ROpt":
            r_high = a["start"]
    return ValidRegion(z_low_bohr=float(z_low), z_high_bohr=float(z_high), r_high_bohr=float(r_high))


def z_profile(density: np.ndarray, r_bohr: np.ndarray, region: ValidRegion) -> np.ndarray:
    """(nz,): A rho_E integrated over r below the r absorber (dr times the sum; dimensionless: 1 for the incident wave
    alone), at every Z (the caller masks the Z outside the valid region)."""
    inside_r = r_bohr < region.r_high_bohr
    return np.sum(density[:, inside_r], axis=1) * (r_bohr[1] - r_bohr[0])


def profile_change(sums: EnergyResolved, info: RunInfo, site: int, energy: int, region: ValidRegion) -> float:
    """The largest change of a site's profile (z_profile) in the valid Z range between the partial sums (up to
    sums.partial_time_au) and the final ones: an indicator of what ending the time sums at tEnd leaves out (where it
    is large, a longer tEnd changes the result; it does not bound the remainder)."""
    whole = z_profile(area_density(stationary_state(sums, info, site, energy), info), sums.r_bohr, region)
    partial = z_profile(area_density(stationary_state(sums, info, site, energy, True), info), sums.r_bohr, region)
    return float(np.max(np.abs(whole - partial)[in_valid_z(sums.z_bohr, region)]))


def in_valid_z(z_bohr: np.ndarray, region: ValidRegion) -> np.ndarray:
    """(nz,) bool: the Z points of the valid region."""
    return (z_bohr > region.z_low_bohr) & (z_bohr < region.z_high_bohr)
