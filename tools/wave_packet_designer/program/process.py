"""The derived parameters the designer needs: Z_asymptotic, masses, absorber strengths, packet width and momentum.

Shadows: src/input/process.cpp (build_config, process_input, miscellaneous_setup, strength_for_energy,
low_z_absorber_strength), src/report/sections.cpp (absorber_table_warning, reduced_length_text)

Only the scalar quantities of the Z and r design are copied (no lattice, no x/y momenta). The order of
operations is the C++'s, so the values are bitwise the program's.
"""

import math

from program.config import AbsorberStrength, DerivedParameters, InputParameters, low_z_absorber_enabled
from program.findings import Finding, FindingList
from program.point_group import PHI_POINTS_MULTIPLE
from program.potential_systems import SystemParameters
from program.unit_conversion import PI, amu_to_electron_masses, ev_to_hartree, hartree_to_ev, rad_to_deg
from program.wave_packet_params import (
    AT_FIRST_POINT,
    BELOW,
    INSIDE,
    absorber_reduced_length,
    optical_potential_strength,
    wave_packet_alpha,
    wave_packet_average_energy,
)


def effective_z_asymptotic(p: InputParameters, system: SystemParameters) -> float:
    """Z_asymptotic of the run: ZmaxSwitch with UseSwitching, else the potential's (build_config)."""
    return p.z_max_switch if p.use_switching else system.z_asymptotic


def build_derived(p: InputParameters, system: SystemParameters) -> DerivedParameters:
    """The derived parameters of a validated input (build_config: process_input + miscellaneous_setup)."""
    total_mass = amu_to_electron_masses(p.mass1 + p.mass2)
    reduced_mass = amu_to_electron_masses((p.mass1 * p.mass2) / (p.mass1 + p.mass2))

    e_min = ev_to_hartree(p.e_min_wave_packet)
    e_max = ev_to_hartree(p.e_max_wave_packet)
    translational_energy = p.e_parallel + hartree_to_ev(wave_packet_average_energy(e_min, e_max, total_mass))
    # The angle of incidence in degrees first, then radians, as the Fortran does (same rounding).
    theta_degrees = rad_to_deg(math.asin(math.sqrt(p.e_parallel / translational_energy)))
    theta = PI / 180.0 * theta_degrees
    k_total = math.sqrt(ev_to_hartree(2.0 * total_mass * translational_energy))
    kz = k_total * math.cos(theta)  # = k_total for EParallel = 0
    alpha = wave_packet_alpha(e_min, e_max, total_mass)
    k_max = kz + 2.0 * math.sqrt(alpha)  # the top of the packet's 95% interval in momentum

    return DerivedParameters(
        z_asymptotic=effective_z_asymptotic(p, system),
        point_group=system.point_group,
        run_mode=run_mode(p),
        total_mass=total_mass,
        reduced_mass=reduced_mass,
        mj_max_basis=p.j_max if p.mj_max < 1 else p.mj_max,
        phi_points_multiple=PHI_POINTS_MULTIPLE[system.point_group] if p.use_symmetry and p.phis_optimal else 1,
        wave_packet_alpha=alpha,
        kz=kz,
        confidence_e_max=0.5 / total_mass * (k_max * k_max) * hartree_to_ev(1.0),
        absorbers=absorber_strengths(p, total_mass, reduced_mass),
    )


def run_mode(p: InputParameters) -> int:
    """Fortran nread: 0 new run, 1 restart, 2 analysis only, 6 potential only."""
    mode = 0 if p.new_calculation else 1
    if p.only_analysis:
        mode = 2
    if p.only_potential:
        mode = 6
    return mode


def absorber_strengths(p: InputParameters, total_mass: float, reduced_mass: float) -> tuple[AbsorberStrength, ...]:
    """The absorbers in the order of the Fortran's ProcessInput; Z absorbers use the total mass, r the reduced."""
    absorbers = [absorber_strength("ZOpt", p.z_opt_end - p.z_opt_start, p.z_opt_energy, p.z_opt_factor, total_mass)]
    if low_z_absorber_enabled(p.z_opt_low_start, p.z_opt_low_end):
        absorbers.append(
            absorber_strength(
                "ZOptLow", p.z_opt_low_end - p.z_opt_low_start, p.z_opt_low_energy, p.z_opt_low_factor, total_mass
            )
        )
    absorbers.append(
        absorber_strength("ROpt", p.r_opt_end - p.r_opt_start, p.r_opt_energy, p.r_opt_factor, reduced_mass)
    )
    absorbers.append(
        absorber_strength("ZSpOpt", p.zsp_opt_end - p.zsp_opt_start, p.zsp_opt_energy, p.zsp_opt_factor, total_mass)
    )
    if low_z_absorber_enabled(p.zsp_opt_low_start, p.zsp_opt_low_end):
        absorbers.append(
            absorber_strength(
                "ZSpOptLow",
                p.zsp_opt_low_end - p.zsp_opt_low_start,
                p.zsp_opt_low_energy,
                p.zsp_opt_low_factor,
                total_mass,
            )
        )
    return tuple(absorbers)


def absorber_strength(
    name: str, length: float, energy: float | None, factor: float | None, mass: float
) -> AbsorberStrength:
    """Strength from the target energy if it is > 0, else the factor (0 when absent).

    The input checks allow an energy of 0 only for ZOptLowEnergy; it gives strength 0, no absorption, as the
    Fortran's low-Z default (docs/fortran_deviations.md, absorber D1).
    """
    if energy is not None and energy > 0.0:
        kinetic_energy = ev_to_hartree(energy)
        reduced = absorber_reduced_length(kinetic_energy, length, mass)
        strength = optical_potential_strength(kinetic_energy, reduced)
        return AbsorberStrength(name, energy, reduced.value, reduced.range, strength)
    return AbsorberStrength(name, None, None, None, factor if factor is not None else 0.0)


def absorber_table_warnings(absorbers: tuple[AbsorberStrength, ...]) -> list[Finding]:
    """The set-up warning for every energy-tuned absorber whose L/lambda is not inside the strength table."""
    see = '(inputs/README.md, section "Absorbers")'
    warnings = FindingList("process_input")
    for a in absorbers:
        if a.table_range is None or a.table_range == INSIDE:
            continue
        key = f"{a.name}Energy"
        if a.table_range == AT_FIRST_POINT:
            warnings.warning(
                key,
                "at this energy the absorber is exactly 1 wavelength long (L/lambda = 1), the first point of the "
                "Vibok-Balint-Kurti table, which no interval of the table includes: its strength is 0 and the "
                "absorber is off (the Fortran's strength is undefined here: docs/fortran_deviations.md, absorber D2) "
                + see,
            )
            continue
        below = a.table_range == BELOW
        warnings.warning(
            key,
            f"at this energy the absorber is {reduced_length_text(a.reduced_length)} wavelengths long (L/lambda), "
            f"{'below' if below else 'above'} the Vibok-Balint-Kurti table (L/lambda = 1 .. 15), so its strength "
            f"is the table's value at L/lambda = {1 if below else 15} {see}",
        )
    return warnings.findings


def reduced_length_text(reduced_length: float) -> str:
    """L/lambda with 6 significant digits, or all digits where 6 would print a table end (1 or 15) for a value
    that is not that end."""
    text = f"{reduced_length:.6g}"
    looks_like_an_end = (text == "1" and reduced_length != 1.0) or (text == "15" and reduced_length != 15.0)
    return repr(reduced_length) if looks_like_an_end else text
