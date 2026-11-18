"""The recommendations of inputs/README.md, each applied to the input by one small function.

Each rule returns a DesignCheck with the input's value, the recommended limit and the README section, or None
where it does not apply. To add a rule: write the function, add it to RULES, and add a failing row to
FAILING in tests/test_design.py. Values are compared with the limits by the program's tolerant comparisons
(absolute 1e-12), so a limit the README states as inclusive holds at its boundary. Atomic units inside; the
checks carry the units the README uses.
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

from program.config import DerivedParameters, InputParameters
from program.gas_phase import GasPhase
from program.real_compare import req, rge, rle
from program.spatial_grids import SpatialGrids
from program.unit_conversion import hartree_to_ev

TURNING_POINT = 2.0  # bohr: the classical turning point Z_turn of section 8
BOUND_STATE_TAIL = 1e-6  # the designer's limit for the norm of an analysed bound state beyond ROptStart


@dataclass(frozen=True)
class DesignInput:
    """What the rules look at: the input, its derived parameters and grids, and the gas phase (if computed)."""

    p: InputParameters
    d: DerivedParameters
    grids: SpatialGrids
    gas: GasPhase | None


@dataclass(frozen=True)
class DesignCheck:
    """One recommendation applied to the input."""

    name: str  # short name, e.g. "dZ: (k_grid - k0) / sigma_k"
    passed: bool | None  # None: information only, no limit
    value: float
    relation: str | None  # how the value must compare with the limit: ">=", "<=" or "="; None: information
    limit: float | None
    unit: str  # of value and limit, "" if none
    advice: str  # what the rule protects, or what to change when it fails
    section: str  # the inputs/README.md section, e.g. "4. dZ"
    designer_limit: bool = False  # the limit is the designer's own choice, not a number of the README


def meets(value: float, relation: str, limit: float) -> bool:
    """value >= limit, <= limit or = limit, within the program's absolute tolerance 1e-12."""
    compare = {">=": rge, "<=": rle, "=": req}[relation]
    return compare(value, limit)


def checked(
    name: str,
    value: float,
    relation: str,
    limit: float,
    unit: str,
    advice: str,
    section: str,
    designer_limit: bool = False,
) -> DesignCheck:
    """A DesignCheck whose `passed` follows from value, relation and limit."""
    passed = meets(value, relation, limit)
    return DesignCheck(name, passed, value, relation, limit, unit, advice, section, designer_limit)


def gaussian_tail(n_sigma: float) -> float:
    """The one-sided tail of a Gaussian beyond n_sigma standard deviations."""
    return 0.5 * math.erfc(n_sigma / math.sqrt(2.0))


def packet_width_k(c: DesignInput) -> float:
    """sigma_k = sqrt(alpha) (1/bohr), the width of |psi(k)|^2."""
    return math.sqrt(c.d.wave_packet_alpha)


def energy_range(c: DesignInput) -> DesignCheck:
    """k0 / sigma_k >= 6, i.e. EMaxWP <= 4 EMinWP: the packet itself has (almost) no outward momentum."""
    n_sigma = c.d.kz / packet_width_k(c)
    advice = (
        f"outward-moving part of the packet {gaussian_tail(n_sigma):.1e}; it reaches the specular absorber and "
        "stops the specular analysis line: narrow the energy range"
    )
    return checked("energy range: k0 / sigma_k", n_sigma, ">=", 6.0, "", advice, "3. Wave packet")


def dz_holds_the_momenta(c: DesignInput) -> DesignCheck:
    """(k_grid - k0) / sigma_k >= 6 with k_grid = pi/dZ: hardly any momentum wraps around to outward motion."""
    n_sigma = (math.pi / c.p.dz - c.d.kz) / packet_width_k(c)
    advice = (
        f"wrapped to outward motion {gaussian_tail(n_sigma):.1e} (must be << TolSp): a smaller dZ or a "
        "narrower energy range"
    )
    return checked("dZ: (k_grid - k0) / sigma_k", n_sigma, ">=", 6.0, "", advice, "4. dZ")


def analysis_line_on_the_first_point(c: DesignInput) -> DesignCheck:
    """Zana is the first grid point >= Z_asymptotic: any further out only lengthens the grid and the run."""
    first = c.grids.z_scattering.point(c.grids.nz_potential)
    excess = c.grids.z_analysis - first
    advice = f"the first grid point >= Z_asymptotic is {first:.4f} bohr"
    return checked("Zana above the first possible point", excess, "<=", 0.0, "bohr", advice, "2. Recipe")


def z0_margin(c: DesignInput) -> DesignCheck:
    """Z0 >= Zana + 8 sigma_Z: the packet keeps away from the region where the two grids are coupled."""
    sigma_z = 1.0 / (2.0 * packet_width_k(c))
    n_sigma = (c.p.z0 - c.grids.z_analysis) / sigma_z
    advice = f"Z0 >= {c.grids.z_analysis + 8.0 * sigma_z:.2f} bohr"
    return checked("Z0 - Zana in sigma_Z", n_sigma, ">=", 8.0, "", advice, "5. Z0")


def outer_tail_gap(c: DesignInput) -> DesignCheck:
    """The packet's outer tail (Z0 + 5.2 sigma_Z, the 1e-7 point) lies >= ~2 bohr below ZSpOptStart."""
    sigma_z = 1.0 / (2.0 * packet_width_k(c))
    gap = c.p.zsp_opt_start - (c.p.z0 + 5.2 * sigma_z)
    advice = f"ZSpOptStart >= {c.p.z0 + 5.2 * sigma_z + 2.0:.2f} bohr"
    return checked("ZSpOptStart - outer tail of the packet", gap, ">=", 2.0, "bohr", advice, "5. Z0")


def scattering_absorber_length(c: DesignInput) -> DesignCheck:
    """The scattering absorber is >= ~4 bohr long (H2): it must absorb every channel down to ~0.05 eV."""
    length = c.p.z_opt_end - c.p.z_opt_start
    return checked(
        "scattering absorber length", length, ">=", 4.0, "bohr", "a longer scattering grid (nZ)", "7. Absorbers"
    )


def absorber_ends(c: DesignInput) -> list[DesignCheck]:
    """Each absorber ends on the last grid point: whatever passes the end wraps around to the grid's start."""
    g = c.grids
    ends = [("ZOptEnd", c.p.z_opt_end, g.z_scattering.last), ("ZSpOptEnd", c.p.zsp_opt_end, g.z_specular.last)]
    if c.p.nr > 1:
        ends.append(("ROptEnd", c.p.r_opt_end, g.r.last))
    return [
        checked(
            f"{key} on the last grid point",
            last - end,
            "=",
            0.0,
            "bohr",
            f"{key} = {last:.10g}, the last grid point",
            "7. Absorbers",
        )
        for key, end, last in ends
    ]


def r_absorber_beyond_bound_states(c: DesignInput) -> DesignCheck | None:
    """ROptStart lies beyond the analysed bound states: their norm beyond it is small. The README asks for
    "beyond the extent"; the number BOUND_STATE_TAIL is the designer's. States above the curve's last value
    are not bound on its grid (continuum): a larger ROptStart cannot contain them, so they are left out."""
    gas = c.gas
    if gas is None or c.p.nr == 1:
        return None
    dr = gas.r[1] - gas.r[0]
    largest = 0.0
    for v in range(gas.n_v_analysis):
        for j in range(gas.j_max_analysis + 1):
            if gas.energies[v, j] > gas.curve.potential[-1]:
                continue
            psi = gas.wave_functions[v, j]
            beyond = sum(psi[i] * psi[i] for i in range(len(psi)) if gas.r[i] >= c.p.r_opt_start) * dr
            largest = max(largest, float(beyond))
    return checked(
        "norm of analysed bound states beyond ROptStart",
        largest,
        "<=",
        BOUND_STATE_TAIL,
        "",
        "a larger ROptStart (and r grid)",
        "7. Absorbers",
        designer_limit=True,
    )


def lowest_converged_energy(c: DesignInput) -> DesignCheck:
    """From tEnd: T = (Z0 - Z_turn) M / k_min + (Zana - Z_turn) M / k_n gives the lowest outgoing normal energy
    E_n that has crossed Zana by tEnd (information: slower flux is counted as reaction)."""
    mass = c.d.total_mass
    k_min = c.d.kz - 2.0 * packet_width_k(c)  # the slow end of the packet
    time_left = c.p.t_end - (c.p.z0 - TURNING_POINT) * mass / k_min
    if time_left <= 0.0:
        energy = math.inf
    else:
        k_n = (c.grids.z_analysis - TURNING_POINT) * mass / time_left
        energy = hartree_to_ev(k_n * k_n / (2.0 * mass))
    advice = "outgoing normal energies below it are not converged; check the norm left at tEnd"
    return DesignCheck(
        "lowest converged outgoing normal energy", None, energy, None, None, "eV", advice, "8. Total time"
    )


def energy_resolution(c: DesignInput) -> DesignCheck:
    """The energy resolution of the analysis, about 2 pi / tEnd (information)."""
    resolution = hartree_to_ev(2.0 * math.pi / c.p.t_end)
    advice = "a finer NrOfEs spacing only interpolates"
    return DesignCheck("energy resolution 2 pi / tEnd", None, resolution, None, None, "eV", advice, "8. Total time")


# The rules in report order: each returns a check, a list of checks, or None where it does not apply.
RULES: list[Callable[[DesignInput], DesignCheck | list[DesignCheck] | None]] = [
    energy_range,
    dz_holds_the_momenta,
    analysis_line_on_the_first_point,
    z0_margin,
    outer_tail_gap,
    scattering_absorber_length,
    absorber_ends,
    r_absorber_beyond_bound_states,
    lowest_converged_energy,
    energy_resolution,
]


def design_checks(c: DesignInput) -> tuple[DesignCheck, ...]:
    """Every rule applied to the input, in the order of RULES."""
    checks: list[DesignCheck] = []
    for rule in RULES:
        result = rule(c)
        if isinstance(result, list):
            checks += result
        elif result is not None:
            checks.append(result)
    return tuple(checks)
