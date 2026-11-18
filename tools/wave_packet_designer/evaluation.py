"""evaluate(): everything the designer knows about one input, computed once.

Shadows: src/app/run.cpp (run_program)

The model layer (README, "Architecture"): the front ends call `evaluate` and only draw or print its
`Evaluation`; they compute nothing themselves. evaluate() runs the set-up stages of run_program that
follow the reading of the file: validate_input, the derived parameters, the grids, the flux refusal and
the initial wave packet. A stage runs only when the ones before it passed, so `stopped_at` names the stage
where spodvr would stop, and `findings` lists that stage's errors (all of them, where spodvr shows the
first).

Two things happen outside the program's stages:
- Reading the file (program/json_input.py, parse_input) comes before evaluate(): a front end reads the
  file and calls evaluate() only when spodvr would get through the reading.
- A CRP data set brings its own Z_asymptotic and point group, which the designer gets from its settings
  (settings.py). Without them evaluate() stops at the stage SETTINGS_STAGE, which is not a program stage.
"""

import math
from dataclasses import dataclass

from program.config import DerivedParameters, InputParameters
from program.findings import Finding, Severity, has_errors
from program.initial_wave_packet import InitialWavePacket, make_initial_wave_packet
from program.potential_systems import SystemParameters, is_data_set_system, system_parameters
from program.process import absorber_table_warnings, build_derived, effective_z_asymptotic
from program.spatial_grids import IndexRange, SpatialGrids, make_spatial_grids
from program.validate import validate_input
from settings import Settings

SETTINGS_STAGE = "settings"  # the designer lacks a value that is not an input key; spodvr itself is not involved


@dataclass(frozen=True)
class AbsorberSummary:
    """One absorber with its grid points, as the reports show it."""

    name: str  # the key prefix, e.g. "ZSpOpt"
    grid: str  # "Z scattering", "Z specular" or "r"
    start: float  # bohr
    end: float  # bohr
    indices: IndexRange  # grid points inside
    target_energy: float | None  # eV, when tuned by its …Energy key
    reduced_length: float | None  # L / lambda at the target energy
    strength: float  # hartree (0: off)


@dataclass(frozen=True)
class PacketSummary:
    """The initial wave packet in the terms of inputs/README.md, sections 3 and 4 (atomic units)."""

    sigma_k: float  # 1/bohr: sqrt(alpha), the width of |psi(k)|^2
    sigma_z: float  # bohr: 1 / (2 sigma_k), the width of |psi(Z)|^2
    grid_momentum_limit: float  # 1/bohr: k_grid = pi / dZ, the largest |k| the Z grids hold
    # From the sampled packet; None in the run modes that build none (onlyAnalysis, Potential).
    norm_outside: float | None  # norm - norm inside: the part of the packet outside the allowed region
    tolsp_measure: float | None  # 1 - norm inside: what the program compares with TolSp
    edge_to_maximum: float | None  # |psi(k)| at the grid edge / its maximum


@dataclass(frozen=True)
class Evaluation:
    """The result of evaluate(); a field keeps its default (None or empty) when its stage did not run."""

    parameters: InputParameters
    settings: Settings
    system: SystemParameters | None = None
    z_asymptotic: float | None = None  # bohr: the one the run uses (ZmaxSwitch with UseSwitching)
    derived: DerivedParameters | None = None
    grids: SpatialGrids | None = None
    absorbers: tuple[AbsorberSummary, ...] = ()
    packet: InitialWavePacket | None = None
    packet_summary: PacketSummary | None = None
    findings: tuple[Finding, ...] = ()  # the program's errors and warnings, stage by stage
    stopped_at: str | None = None  # the stage where evaluation stops; None if the set-up passes

    @property
    def accepted(self) -> bool:
        """True if the program gets through the set-up (and builds the initial wave packet, if it propagates)."""
        return self.stopped_at is None


def evaluate(p: InputParameters, settings: Settings) -> Evaluation:
    """Runs the program's set-up stages on `p` as far as they go."""
    system = system_parameters(p.system_id, settings.z_asymptotic, settings.point_group)
    if system is None:
        return missing_system_parameters(p, settings)

    z_asymptotic = effective_z_asymptotic(p, system)
    findings = validate_input(p, system, z_asymptotic)
    if has_errors(findings):
        return Evaluation(p, settings, system, z_asymptotic, findings=tuple(findings), stopped_at="validate_input")

    derived = build_derived(p, system)
    findings += absorber_table_warnings(derived.absorbers)
    grids, grid_findings = make_spatial_grids(p, derived)
    findings += grid_findings
    if grids is None:
        return Evaluation(
            p, settings, system, z_asymptotic, derived, findings=tuple(findings), stopped_at="make_spatial_grids"
        )
    absorbers = absorber_summaries(p, derived, grids)

    # Potential-only and analysis-only runs build no wave packet; a propagating run refuses the flux analysis.
    propagates = derived.run_mode in (0, 1)
    packet, stopped_at = None, None
    if p.do_flux_analysis and propagates:
        message = "doFluxAnalysis = true: the flux analysis is not implemented yet"
        findings.append(Finding(Severity.ERROR, "doFluxAnalysis", message, "run_program"))
        stopped_at = "run_program"
    elif propagates:
        packet, packet_findings = make_initial_wave_packet(p, derived, grids)
        findings += packet_findings
        if has_errors(packet_findings):
            stopped_at = "make_initial_wave_packet"
    return Evaluation(
        p,
        settings,
        system,
        z_asymptotic,
        derived,
        grids,
        absorbers,
        packet=packet,
        packet_summary=packet_summary(p, derived, packet),
        findings=tuple(findings),
        stopped_at=stopped_at,
    )


def missing_system_parameters(p: InputParameters, settings: Settings) -> Evaluation:
    """Without Z_asymptotic and the point group nothing else can be computed: an unknown SystemID (spodvr
    stops in validate_input), or a CRP data set whose values are not in the settings."""
    if not is_data_set_system(p.system_id):
        finding = Finding(Severity.ERROR, "SystemID", "Unknown system ID", "validate_input")
        return Evaluation(p, settings, findings=(finding,), stopped_at="validate_input")
    message = (
        f"SystemID {p.system_id}: Z_asymptotic and the point group come from the CRP data set and must be given "
        "to the designer (the program's set-up output prints both)"
    )
    finding = Finding(Severity.ERROR, "SystemID", message, "designer settings")
    return Evaluation(p, settings, findings=(finding,), stopped_at=SETTINGS_STAGE)


def absorber_summaries(p: InputParameters, d: DerivedParameters, g: SpatialGrids) -> tuple[AbsorberSummary, ...]:
    """Every absorber of the run with its grid points: ZOpt, ROpt and ZSpOpt always, the low-Z ones when on."""
    # Where each absorber acts: grid, start, end and the grid points inside (make_spatial_grids).
    placement = {
        "ZOpt": ("Z scattering", p.z_opt_start, p.z_opt_end, g.z_optical),
        "ZOptLow": ("Z scattering", p.z_opt_low_start, p.z_opt_low_end, g.z_optical_low),
        "ROpt": ("r", p.r_opt_start, p.r_opt_end, g.r_optical),
        "ZSpOpt": ("Z specular", p.zsp_opt_start, p.zsp_opt_end, g.z_specular_optical),
        "ZSpOptLow": ("Z specular", p.zsp_opt_low_start, p.zsp_opt_low_end, g.z_specular_optical_low),
    }
    summaries = []
    for a in d.absorbers:  # ZOpt, ROpt and ZSpOpt always; the low-Z ones only when given
        grid, start, end, indices = placement[a.name]
        summaries.append(
            AbsorberSummary(
                name=a.name,
                grid=grid,
                start=start,
                end=end,
                indices=indices,
                target_energy=a.target_energy,
                reduced_length=a.reduced_length,
                strength=a.strength,
            )
        )
    return tuple(summaries)


def packet_summary(p: InputParameters, d: DerivedParameters, packet: InitialWavePacket | None) -> PacketSummary:
    """The packet's widths and momentum limit, and with a sampled packet its norms and |psi(k)| at the edge."""
    sigma_k = math.sqrt(d.wave_packet_alpha)
    sigma_z = 1.0 / (2.0 * sigma_k)
    grid_momentum_limit = math.pi / p.dz
    if packet is None:
        return PacketSummary(sigma_k, sigma_z, grid_momentum_limit, None, None, None)
    info = packet.info
    return PacketSummary(
        sigma_k,
        sigma_z,
        grid_momentum_limit,
        norm_outside=info.norm - info.norm_inside,
        tolsp_measure=1.0 - info.norm_inside,
        edge_to_maximum=info.k_space_at_edge / info.k_space_maximum,
    )
