"""The input parameters and the derived parameters, as far as the designer uses them.

Shadows: src/input/config.hpp (InputParameters, DerivedParameters)

InputParameters is flat like the C++ struct: one field per input key, the JSON key and unit in the
comment, the C++ default as the default. Keys the designer only type-checks (Title, the output switches
such as writeSMatrix, the SCM and flux keys, …) have no field here; json_input.py lists every key. To add a
key: a field here and a row in json_input.py (tests/test_input_reading.py checks that every field has its
row).
Lengths are in bohr, energies in eV, masses in amu, times in atomic units.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class InputParameters:
    """The input keys the designer uses (config.hpp InputParameters)."""

    system_id: str  # SystemID
    new_calculation: bool  # NewCalc
    z_dof: bool  # ZDoF (the six DoF keys must all be true: 6D only)
    r_dof: bool  # RDoF
    x_dof: bool  # XDoF
    y_dof: bool  # YDoF
    theta_dof: bool  # ThetaDoF
    phi_dof: bool  # PhiDoF

    mass1: float  # Mass1 (amu)
    mass2: float  # Mass2 (amu)
    homonuclear: bool  # HomoNuc

    # Z grids: both start at z_start with spacing dz; the specular grid has more points.
    z_start: float  # ZStart (bohr)
    dz: float  # dZ (bohr)
    nz: int  # nZ: points of the scattering grid
    nz_specular: int  # nZSpecular: points of the specular grid (the initial channel only)
    r_start: float  # RStart (bohr)
    dr: float  # dR (bohr)
    nr: int  # nR
    nx: int  # nX (x and y enter only the memory estimate)
    ny: int  # nY
    j_max: int  # JMax
    projection_operator: int  # ProjOp (only 2 is supported)

    t_start: float  # tStart (a.u.)
    t_end: float  # tEnd (a.u.)
    time_step: float  # dtStep (a.u.)
    n_steps: int  # nStep

    dt_analysis: float  # dtAna (a.u.)
    dynamical_zsp_analysis: bool  # DynamZSpAna: the specular analysis line follows the packet
    dt_averages: float  # dtAve (a.u.)
    z_analysis: float  # Zana (bohr): the analysis line on the scattering grid
    n_energies: int  # NrOfEs (>= 2)

    e_min_wave_packet: float  # EMinWP (eV)
    e_max_wave_packet: float  # EMaxWP (eV)
    z0: float  # Z0 (bohr): centre of the initial packet

    # Absorbers: start, end (bohr) and either a target energy (eV) or a strength factor (hartree).
    z_opt_start: float  # ZOptStart
    z_opt_end: float  # ZOptEnd
    r_opt_start: float  # ROptStart
    r_opt_end: float  # ROptEnd
    zsp_opt_start: float  # ZSpOptStart
    zsp_opt_end: float  # ZSpOptEnd

    e_parallel: float  # EParallel (eV)
    phi_parallel: float  # PhiParallel (degrees)
    tol_specular: float  # TolSp: allowed norm of the initial packet outside its region

    # Optional keys (C++ defaults).
    mj_max: int = -1  # MJMax: absent or negative is replaced by JMax on reading (the basis also treats 0 so)
    only_analysis: bool = False  # onlyAnalysis (run mode 2)
    also_analysis: bool = False  # alsoAnalysis
    j_analysis: int | None = None  # JAna
    v_analysis: int | None = None  # VAna
    diffraction_order_analysis: int | None = None  # DifOrdAna
    z_analysis_specular: float | None = None  # ZAnaSp (bohr), only with DynamZSpAna = false
    do_flux_analysis: bool = False  # doFluxAnalysis (not ported: a propagating run stops)
    v_in: int = 0  # VIn
    j_in: int = 0  # JIn
    mj_in: int = 0  # MJIn
    only_potential: bool = False  # Potential (run mode 6)
    impose_v_min: bool = False  # ImposeVMin
    v_min: float = -5.5  # VMin (eV)
    impose_v_max: bool = False  # ImposeVMax
    v_max: float = 10.0  # VMax (eV)
    use_switching: bool = False  # UseSwitching
    z_min_switch: float = 6.0  # ZminSwitch (bohr)
    z_max_switch: float = 8.0  # ZmaxSwitch (bohr): replaces Z_asymptotic with UseSwitching
    write_2d_pes: bool = False  # write2DPES: the run writes V(r, Z) at one (x, y, theta, phi) and stops
    x_2d_pes: float = 0.0  # x2DPES (bohr, along the first lattice vector)
    y_2d_pes: float = 0.0  # y2DPES (bohr, along the second lattice vector)
    theta_2d_pes: float = 90.0  # theta2DPES (degrees)
    phi_2d_pes: float = 0.0  # phi2DPES (degrees)
    use_symmetry: bool = False  # UseSymmetry
    phis_optimal: bool = False  # nPhisOptimal
    symmetry_tolerance: float = 1e-5  # SymmetryTolerance (eV)
    crp_dataset: str = ""  # CRPlibDataSet

    z_opt_energy: float | None = None  # ZOptEnergy (eV)
    z_opt_factor: float | None = None  # ZOptFactor (hartree)
    r_opt_energy: float | None = None  # ROptEnergy
    r_opt_factor: float | None = None  # ROptFactor
    zsp_opt_energy: float | None = None  # ZSpOptEnergy
    zsp_opt_factor: float | None = None  # ZSpOptFactor
    # Low-Z absorbers (strongest at their start); all four keys absent = off.
    z_opt_low_start: float | None = None  # ZOptLowStart
    z_opt_low_end: float | None = None  # ZOptLowEnd
    z_opt_low_energy: float | None = None  # ZOptLowEnergy
    z_opt_low_factor: float | None = None  # ZOptLowFactor
    zsp_opt_low_start: float | None = None  # ZSpOptLowStart
    zsp_opt_low_end: float | None = None  # ZSpOptLowEnd
    zsp_opt_low_energy: float | None = None  # ZSpOptLowEnergy
    zsp_opt_low_factor: float | None = None  # ZSpOptLowFactor

    restart_dir: str = ""  # RestartDir
    restart_format: str = "native"  # RestartFormat
    checkpoint_format: str = "native"  # CheckpointFormat


def low_z_absorber_enabled(start: float | None, end: float | None) -> bool:
    """A low-Z absorber is on when its start and end are given (LowZAbsorberInput::enabled)."""
    return start is not None and end is not None


@dataclass(frozen=True)
class AbsorberStrength:
    """How the strength of one absorber was obtained (the factors and AbsorberTableWarning of DerivedParameters).

    The C++ keeps only the strength and, for an L/lambda outside the table, a warning; the designer keeps
    L/lambda for every energy-tuned absorber so the reports need not recompute it.
    """

    name: str  # the key prefix: ZOpt, ZOptLow, ROpt, ZSpOpt or ZSpOptLow
    target_energy: float | None  # eV, when tuned by its …Energy key (only an energy > 0 counts)
    reduced_length: float | None  # L/lambda at the target energy
    table_range: str | None  # where L/lambda lies in the table: below, inside, above or at_first_point
    strength: float  # hartree; 0 means the absorber is off


@dataclass(frozen=True)
class DerivedParameters:
    """The derived parameters the designer uses (config.hpp DerivedParameters, atomic units)."""

    z_asymptotic: float  # bohr: the PES's, or ZmaxSwitch with UseSwitching
    point_group: str
    run_mode: int  # 0 new, 1 restart, 2 analysis only, 6 potential only
    total_mass: float  # electron masses
    reduced_mass: float  # electron masses
    mj_max_basis: int  # MJMax of the rotational basis: MJMax < 1 means JMax
    phi_points_multiple: int  # the number of phi points is a multiple of this (1 unless UseSymmetry and nPhisOptimal)
    wave_packet_alpha: float  # bohr^-2: psi(Z) ~ exp(-alpha (Z - Z0)^2)
    kz: float  # bohr^-1: the packet moves with momentum -kz
    confidence_e_max: float  # eV: kinetic energy at kz + 2 sqrt(alpha), the top of the packet's 95% interval
    # The absorbers in the order of the Fortran's ProcessInput (the order of its warnings): ZOpt, ZOptLow
    # (if on), ROpt, ZSpOpt, ZSpOptLow (if on).
    absorbers: tuple[AbsorberStrength, ...]
