"""evaluate(): everything the designer knows about one input, computed once.

Shadows: src/app/run.cpp (run_program)

The model layer (README, "Architecture"): the front ends call `evaluate` and only draw or print its
`Evaluation`; they compute nothing themselves. evaluate() runs the set-up stages of run_program that
follow the reading of the file: validate_input, the derived parameters, the grids, the flux refusal and
the initial wave packet. A stage runs only when the ones before it passed, so `stopped_at` names the stage
where spodvr would stop, and `findings` lists that stage's errors (all of them, where spodvr shows the
first). Once the grids exist it also computes what the program prints later in its set-up (the memory
estimates and the gas-phase channels) and the advice of inputs/README.md (design/: the design checks and the
absorbers' R + T). suggest() applies the README's recipe to the Z keys and evaluates the changed input.

Two things happen outside the program's stages:
- Reading the file (program/json_input.py, parse_input) comes before evaluate(): a front end reads the
  file and calls evaluate() only when spodvr would get through the reading.
- A CRP data set brings its own Z_asymptotic and point group, which the designer gets from its settings
  (settings.py). Without them evaluate() stops at the stage SETTINGS_STAGE, which is not a program stage.
"""

import math
from dataclasses import dataclass, replace

import numpy as np

from design.rules import (
    INFORMATION,
    LOWEST_OUTGOING_ENERGY,
    TURNING_POINT,
    DesignCheck,
    DesignInput,
    analysed_bound_states,
    design_checks,
    memory_fits,
    norm_beyond,
)
from design.suggest import SUGGESTED_FIELDS, fit_r_grid, fit_time, fit_z_grids, suggested_keys
from design.timing import TimeBudget, time_budget
from design.transmission import AbsorberResponse, absorber_mass, absorber_response, response_window
from dz_gas_phase import dz_gas_phase_curve
from program.angular_grid import AngularGrid, make_angular_grid, rotational_basis_size
from program.config import DerivedParameters, InputParameters
from program.fft_sizes import nearest_fft_size
from program.findings import Finding, FindingList, Severity, has_errors
from program.gas_phase import GasPhase, GasPhaseCurve, compute_gas_phase
from program.gas_phase_potentials import gas_phase_potential
from program.initial_wave_packet import InitialWavePacket, make_initial_wave_packet
from program.json_input import key_of_field
from program.memory_estimate import MemoryEstimate, SymmetryCounts, estimate_gpu_memory, estimate_memory
from program.point_group import GROUP_ORDER
from program.potential_systems import SystemParameters, is_data_set_system, system_parameters
from program.process import absorber_table_warnings, build_derived, effective_z_asymptotic
from program.spatial_grids import IndexRange, SpatialGrids, UniformGrid, make_spatial_grids
from program.unit_conversion import hartree_to_ev
from program.validate import validate_input
from settings import Settings

REFERENCE_REFINEMENT = 4  # the converged gas phase uses this many times the input's r points, over the same range
SETTINGS_STAGE = "settings"  # the designer lacks a value that is not an input key; spodvr itself is not involved
# Where the surface-symmetry counts of the memory estimate come from (Resources.symmetry_source).
SYMMETRY_OFF = "off"  # UseSymmetry is false: every (x, y) point is its own representative
SYMMETRY_GIVEN = "given"  # the representatives were given in the settings
SYMMETRY_NOT_GIVEN = "not given"  # UseSymmetry without counts: no reduction assumed (an upper bound)


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
class PacketDensities:
    """The initial packet as the figure draws it: |psi(Z)|^2 and |psi(k)|^2 (ascending k)."""

    z: np.ndarray  # bohr: the specular grid's points
    z_density: np.ndarray  # 1/bohr: |psi(Z)|^2
    k: np.ndarray  # 1/bohr: the grid's momenta, ascending
    k_density: np.ndarray  # bohr: |psi(k)|^2 with the program's scale dZ / sqrt(2 pi)
    outward_fraction: float  # the part of the norm at k > 0: moving outwards (own tail plus wrap-around)


@dataclass(frozen=True)
class Resources:
    """The sizes that set the memory of a run, and its CPU and GPU memory."""

    angular: AngularGrid
    rotational_basis_size: int  # (J, mJ) functions per (x, y, Z, r) point of the wave function
    symmetry: SymmetryCounts
    symmetry_source: str  # SYMMETRY_OFF, SYMMETRY_GIVEN or SYMMETRY_NOT_GIVEN
    threads: int
    memory: MemoryEstimate | None  # None for an analysis-only run (its fine-analysis estimate is not copied)
    gpu_memory: MemoryEstimate | None  # only for a propagating run


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
    resources: Resources | None = None
    gas_phase: GasPhase | None = None  # as the program computes it; None if VIn lies beyond the curve's levels
    gas_reference: GasPhase | None = None  # the same potential on a finer r grid: the converged states
    # The curve is another potential's, or the input's r grid is not the curve's: then the program's states
    # differ, and these are only indicative.
    gas_phase_indicative: bool = False
    # Per analysed v: the norm beyond each r of its most extended analysed bound j, on the gas phase's r grid.
    bound_state_tails: tuple[np.ndarray, ...] = ()
    reference_tails: tuple[np.ndarray, ...] = ()  # the same for the converged states (gas_reference)
    densities: PacketDensities | None = None
    design: tuple[DesignCheck, ...] = ()  # the recommendations of inputs/README.md, applied
    time_budget: TimeBudget | None = None  # how long the slow parts of the packet need (section 8)
    absorber_responses: tuple[AbsorberResponse, ...] = ()  # R and T of each absorber that is on
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
    settings_errors = check_settings(p, system, settings)
    if settings_errors:
        return Evaluation(p, settings, system, z_asymptotic, findings=tuple(settings_errors), stopped_at=SETTINGS_STAGE)

    derived = build_derived(p, system)
    findings += absorber_table_warnings(derived.absorbers)
    grids, grid_findings = make_spatial_grids(p, derived)
    findings += grid_findings
    if grids is None:
        return Evaluation(
            p, settings, system, z_asymptotic, derived, findings=tuple(findings), stopped_at="make_spatial_grids"
        )
    absorbers = absorber_summaries(p, derived, grids)

    # A propagating run refuses the flux analysis right after the grids. A write2DPES run does not propagate: it
    # writes its cut of the PES after the grids and stops (run.cpp), so it builds no packet either.
    propagates = derived.run_mode in (0, 1) and not p.write_2d_pes
    if p.do_flux_analysis and propagates:
        message = "doFluxAnalysis = true: the flux analysis is not implemented yet"
        findings.append(Finding(Severity.ERROR, "doFluxAnalysis", message, "run_program"))
        summary = packet_summary(p, derived, None)
        return Evaluation(
            p,
            settings,
            system,
            z_asymptotic,
            derived,
            grids,
            absorbers,
            None,
            summary,
            findings=tuple(findings),
            stopped_at="run_program",
        )

    # Then the program prints the memory estimates and the gas phase, and a propagating run builds the packet.
    resources = estimate_resources(p, derived, grids, settings)
    gas_phase, gas_reference, gas_phase_indicative = gas_phases(p, derived, grids, resources.angular)
    z_turn, e_n_min = time_settings(settings)
    design_input = DesignInput(p, derived, grids, resources.angular, gas_phase, gas_reference, z_turn, e_n_min)
    packet, stopped_at = None, None
    if propagates:
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
        resources=resources,
        gas_phase=gas_phase,
        gas_reference=gas_reference,
        gas_phase_indicative=gas_phase_indicative,
        bound_state_tails=bound_state_tails(gas_phase, resources),
        reference_tails=bound_state_tails(gas_reference, resources),
        densities=packet_densities(packet, grids),
        design=design_checks(design_input) + memory_checks(resources, settings),
        time_budget=time_budget(design_input),
        absorber_responses=absorber_responses(p, derived, grids, absorbers, gas_phase),
        findings=tuple(findings),
        stopped_at=stopped_at,
    )


@dataclass(frozen=True)
class Suggestion:
    """Suggested Z keys (design/suggest.py) and what the program and the README make of the changed input."""

    changes: dict[str, float | int]  # InputParameters field -> suggested value; empty if blocked
    evaluation: Evaluation  # of the input with the changes (of the input itself if blocked)
    blockers: tuple[Finding, ...]  # the program's errors on keys the suggestion does not set


def suggest(p: InputParameters, settings: Settings, dz: float | None = None) -> Suggestion:
    """The suggestion for the Z grids and absorbers of `p` (with dZ = `dz` if given). It needs an input that
    the program's consistency checks accept apart from the keys it sets: otherwise it is blocked, and
    `blockers` lists what to fix first (e.g. EMinWP <= 0, an unknown SystemID)."""
    original = evaluate(p, settings)
    keys_it_sets = {key_of_field(field) for field in SUGGESTED_FIELDS}
    blockers = tuple(
        f
        for f in original.findings
        if f.severity is Severity.ERROR and f.origin == "validate_input" and f.key not in keys_it_sets
    )
    if original.system is None or blockers:
        return Suggestion({}, original, blockers)
    system = original.system
    z_turn, e_n_min = time_settings(settings)
    changes = suggested_keys(p, build_derived(p, system), effective_z_asymptotic(p, system), dz, e_n_min, z_turn)
    return Suggestion(changes, evaluate(replace(p, **changes), settings), ())


def fitted_parameters(
    p: InputParameters, settings: Settings, scattering_length: float, specular_length: float
) -> InputParameters:
    """`p` with the Z keys that follow from its dZ and ZSpOptStart and the two absorber lengths (bohr), as the
    recipe places them (design/suggest.py fit_z_grids); `p` itself without system parameters."""
    system = system_parameters(p.system_id, settings.z_asymptotic, settings.point_group)
    if system is None:
        return p
    z_asymptotic = effective_z_asymptotic(p, system)
    return replace(p, **fit_z_grids(p, z_asymptotic, scattering_length, specular_length))


def fitted_r_grid(p: InputParameters) -> InputParameters:
    """`p` with the r absorber ending on the last point of its r grid (design/suggest.py fit_r_grid)."""
    return replace(p, **fit_r_grid(p))


def fitted_time(p: InputParameters, t_end: float) -> InputParameters:
    """`p` with tEnd >= `t_end` in whole averaging intervals, nStep to match (design/suggest.py fit_time)."""
    return replace(p, **fit_time(p, t_end))


def absorber_responses(
    p: InputParameters,
    d: DerivedParameters,
    grids: SpatialGrids,
    absorbers: tuple[AbsorberSummary, ...],
    gas: GasPhase | None,
) -> tuple[AbsorberResponse, ...]:
    """R and T of each absorber that is on (strength > 0; the r absorber only with nR > 1), over its energies
    (design/transmission.py response_window; without a gas phase the available energy is EMaxWP + EParallel)."""
    available = hartree_to_ev(gas.available_energy) if gas is not None else p.e_max_wave_packet + p.e_parallel
    responses = []
    for a in absorbers:
        if a.strength <= 0.0 or (a.name == "ROpt" and p.nr == 1):
            continue
        mass = absorber_mass(a.name, d.total_mass, d.reduced_mass)
        window = response_window(a.name, p, available, hartree_to_ev(grids.kr.max_kinetic_energy))
        responses.append(absorber_response(a.name, a.end - a.start, a.strength, mass, window))
    return tuple(responses)


def packet_densities(packet: InitialWavePacket | None, grids: SpatialGrids) -> PacketDensities | None:
    """|psi(Z)|^2 and |psi(k)|^2 of the packet, and its norm at k > 0 (from the sampled spectrum, so it holds
    the packet's own outward tail and whatever wraps around beyond -pi/dZ)."""
    if packet is None:
        return None
    z = grids.z_specular
    k = 2.0 * np.pi * np.fft.fftfreq(z.size, z.spacing)
    order = np.argsort(k)
    k_density = np.abs(packet.psi_k) ** 2
    outward = float(k_density[k > 0.0].sum() / k_density.sum())
    return PacketDensities(np.array(z.points()), np.abs(packet.psi) ** 2, k[order], k_density[order], outward)


def bound_state_tails(gas: GasPhase | None, resources: Resources) -> tuple[np.ndarray, ...]:
    """Per analysed v: at each r the largest norm beyond r over its analysed bound j (design.rules
    analysed_bound_states): what the r absorber would remove if it started there."""
    if gas is None:
        return ()
    dr = gas.r[1] - gas.r[0]
    states = analysed_bound_states(gas, resources.angular)
    tails = []
    for v in sorted({v for v, _ in states}):
        tails.append(np.max([norm_beyond(gas.wave_functions[v, j], dr) for vv, j in states if vv == v], axis=0))
    return tuple(tails)


def memory_checks(resources: Resources, settings: Settings) -> tuple[DesignCheck, ...]:
    """Whether the run fits the memory given in the settings (the RAM, the free GPU memory). With UseSymmetry but
    no representatives given, the estimate is the worst case (no symmetry reduction): it cannot show that a run
    does not fit, so the checks are then information only."""
    checks = []
    if resources.memory is not None and settings.ram_bytes is not None:
        advice = "surface symmetry (UseSymmetry), a smaller grid, or a machine with more memory"
        checks.append(memory_fits("CPU", resources.memory.peak_bytes, settings.ram_bytes, advice))
    if resources.gpu_memory is not None and settings.gpu_bytes is not None:
        advice = "spodvr refuses --device cuda: surface symmetry, a smaller grid, or --device cpu"
        checks.append(memory_fits("GPU", resources.gpu_memory.peak_bytes, settings.gpu_bytes, advice))
    if resources.symmetry_source == SYMMETRY_NOT_GIVEN:
        worst_case = (
            "the worst case, without the symmetry reduction; --representatives (the program's SURFACE SYMMETRY "
            "block) gives the real value"
        )
        checks = [
            replace(
                c, name=c.name + " (worst case: no symmetry)", passed=None, consequence=INFORMATION, advice=worst_case
            )
            for c in checks
        ]
    return tuple(checks)


def check_settings(p: InputParameters, system: SystemParameters, settings: Settings) -> list[Finding]:
    """The designer's own checks of the settings that the input fixes: the symmetry counts must be possible."""
    checks = FindingList("designer settings")
    if p.use_symmetry and settings.representatives is not None:
        n_xy = p.nx * p.ny
        checks.error_if(
            settings.representatives > n_xy,
            "",
            f"{settings.representatives} representatives, but the grid has only {n_xy} (x, y) points",
        )
        order = GROUP_ORDER[system.point_group]
        checks.error_if(
            settings.operations is not None and settings.operations > order,
            "",
            f"{settings.operations} operations, but the point group {system.point_group} has only {order}",
        )
    return checks.findings


def time_settings(settings: Settings) -> tuple[float, float]:
    """The turning point (bohr) and E_n,min (eV) of the time budget: the settings', else design/rules.py's."""
    z_turn = settings.z_turn if settings.z_turn is not None else TURNING_POINT
    energy = settings.lowest_outgoing_energy
    return z_turn, energy if energy is not None else LOWEST_OUTGOING_ENERGY


def gas_phases(
    p: InputParameters, d: DerivedParameters, grids: SpatialGrids, angular: AngularGrid
) -> tuple[GasPhase | None, GasPhase | None, bool]:
    """The gas phase as the program computes it, on the input's r grid, and the converged reference on a
    REFERENCE_REFINEMENT times finer grid over the same range: where the potential's gas phase is copied
    (program/gas_phase_potentials.py). For another potential the DZ input's curve on its own grid stands in,
    without a reference, and the states are only indicative (the third value)."""
    potential = gas_phase_potential(p, grids.r)
    if potential is None:
        curve = dz_gas_phase_curve()
        stand_in = compute_gas_phase(p, d, angular, curve) if p.v_in < len(curve.potential) else None
        return stand_in, None, True
    source = f"the {p.system_id} potential at Z_asymptotic"
    curve = GasPhaseCurve(p.r_start, p.dr, tuple(potential), p.system_id, source + ", on the input's r grid")
    n_fine = nearest_fft_size(REFERENCE_REFINEMENT * p.nr)
    fine_grid = UniformGrid(p.r_start, p.nr * p.dr / n_fine, n_fine)
    fine_potential = tuple(gas_phase_potential(p, fine_grid))
    fine = GasPhaseCurve(p.r_start, fine_grid.spacing, fine_potential, p.system_id, source + ", on a finer r grid")
    return compute_gas_phase(p, d, angular, curve), compute_gas_phase(p, d, angular, fine), False


def estimate_resources(p: InputParameters, d: DerivedParameters, g: SpatialGrids, settings: Settings) -> Resources:
    """The angular sizes and the memory estimates of the run modes that print them: the CPU estimate except
    for an analysis-only run (2), the GPU estimate for a propagating run (0, 1)."""
    angular = make_angular_grid(p, d)
    n_xy = p.nx * p.ny
    if not p.use_symmetry:
        symmetry, source = SymmetryCounts(n_xy, 1), SYMMETRY_OFF
    elif settings.representatives is None:
        symmetry, source = SymmetryCounts(n_xy, 1), SYMMETRY_NOT_GIVEN
    else:
        operations = settings.operations if settings.operations is not None else GROUP_ORDER[d.point_group]
        symmetry, source = SymmetryCounts(settings.representatives, operations), SYMMETRY_GIVEN
    memory = estimate_memory(p, d, g, angular, symmetry, settings.threads) if d.run_mode != 2 else None
    propagates = d.run_mode in (0, 1)
    gpu_memory = estimate_gpu_memory(p, g, angular, symmetry.n_representatives) if propagates else None
    return Resources(angular, rotational_basis_size(angular), symmetry, source, settings.threads, memory, gpu_memory)


def missing_system_parameters(p: InputParameters, settings: Settings) -> Evaluation:
    """Without Z_asymptotic and the point group nothing else can be computed: an unknown SystemID (spodvr
    stops in validate_input), or a CRP data set whose values are not in the settings."""
    if not is_data_set_system(p.system_id):
        finding = Finding(Severity.ERROR, "SystemID", "Unknown system ID", "validate_input")
        return Evaluation(p, settings, findings=(finding,), stopped_at="validate_input")
    message = (
        f"SystemID {p.system_id} takes Z_asymptotic and the point group from its CRP data set: give both (the "
        "program's set-up output prints them)"
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
