"""Reads the input keys like the program: known keys, value types, and the checks made while reading.

Shadows: src/input/json_input.cpp (parse_input, read_required_keys, read_optional_keys, reject_unsupported_keys,
check_while_reading, require_energy_or_factor, require_low_z_absorber, JsonReader)
Shadows: src/input/validate_outputs.cpp (plot_site_fields)

Not copied: check_scm_keys (the SCM keys are type-checked only) and the CRP data set itself.
"""

from dataclasses import dataclass, replace

from program.config import InputParameters, low_z_absorber_enabled
from program.findings import Finding, FindingList, has_errors
from program.real_compare import rlt, seq

# Value types (JsonReader::as_real, as_integer, as_logical, as_string, and the PlotSites list).
REAL, INTEGER, LOGICAL, STRING, SITE_LIST = "real", "integer", "logical", "string", "site list"
SITE_FIELDS = ("name", "x", "y")  # the fields of a PlotSites entry (validate_outputs.cpp plot_site_fields)
PLOT_DATA_KEYS = ("writeDensities", "writePotential", "PlotSites")  # not modelled by the designer

INT_MIN, INT_MAX = -(2**31), 2**31 - 1  # the range of a C++ int


@dataclass(frozen=True)
class InputKey:
    """One input key: its JSON name, its value type, and the InputParameters field that holds it."""

    name: str
    kind: str  # REAL, INTEGER, LOGICAL, STRING or SITE_LIST
    field: str | None  # None: the designer only checks the type


# In the order of read_required_keys.
REQUIRED_KEYS = (
    InputKey("Title", STRING, None),
    InputKey("SystemID", STRING, "system_id"),
    InputKey("NewCalc", LOGICAL, "new_calculation"),
    InputKey("ZDoF", LOGICAL, "z_dof"),
    InputKey("RDoF", LOGICAL, "r_dof"),
    InputKey("XDoF", LOGICAL, "x_dof"),
    InputKey("YDoF", LOGICAL, "y_dof"),
    InputKey("ThetaDoF", LOGICAL, "theta_dof"),
    InputKey("PhiDoF", LOGICAL, "phi_dof"),
    InputKey("Mass1", REAL, "mass1"),
    InputKey("Mass2", REAL, "mass2"),
    InputKey("HomoNuc", LOGICAL, "homonuclear"),
    InputKey("ZStart", REAL, "z_start"),
    InputKey("dZ", REAL, "dz"),
    InputKey("nZ", INTEGER, "nz"),
    InputKey("nZSpecular", INTEGER, "nz_specular"),
    InputKey("RStart", REAL, "r_start"),
    InputKey("dR", REAL, "dr"),
    InputKey("nR", INTEGER, "nr"),
    InputKey("nX", INTEGER, "nx"),
    InputKey("nY", INTEGER, "ny"),
    InputKey("JMax", INTEGER, "j_max"),
    InputKey("ProjOp", INTEGER, "projection_operator"),
    InputKey("tStart", REAL, "t_start"),
    InputKey("tEnd", REAL, "t_end"),
    InputKey("dtStep", REAL, "time_step"),
    InputKey("nStep", INTEGER, "n_steps"),
    InputKey("dtAna", REAL, "dt_analysis"),
    InputKey("DynamZSpAna", LOGICAL, "dynamical_zsp_analysis"),
    InputKey("dtAve", REAL, "dt_averages"),
    InputKey("Zana", REAL, "z_analysis"),
    InputKey("NrOfEs", INTEGER, "n_energies"),
    InputKey("EMinWP", REAL, "e_min_wave_packet"),
    InputKey("EMaxWP", REAL, "e_max_wave_packet"),
    InputKey("Z0", REAL, "z0"),
    InputKey("ZOptStart", REAL, "z_opt_start"),
    InputKey("ZOptEnd", REAL, "z_opt_end"),
    InputKey("ROptStart", REAL, "r_opt_start"),
    InputKey("ROptEnd", REAL, "r_opt_end"),
    InputKey("ZSpOptStart", REAL, "zsp_opt_start"),
    InputKey("ZSpOptEnd", REAL, "zsp_opt_end"),
    InputKey("EParallel", REAL, "e_parallel"),
    InputKey("PhiParallel", REAL, "phi_parallel"),
    InputKey("writeSMatrix", LOGICAL, None),
    InputKey("writeProbs", LOGICAL, None),
    InputKey("writeRDProbs", LOGICAL, None),
    InputKey("TolSp", REAL, "tol_specular"),
)

# In the order of read_optional_keys; absent keys keep the InputParameters default.
OPTIONAL_KEYS = (
    InputKey("MJMax", INTEGER, "mj_max"),
    InputKey("onlyAnalysis", LOGICAL, "only_analysis"),
    InputKey("alsoAnalysis", LOGICAL, "also_analysis"),
    InputKey("JAna", INTEGER, "j_analysis"),
    InputKey("VAna", INTEGER, "v_analysis"),
    InputKey("DifOrdAna", INTEGER, "diffraction_order_analysis"),
    InputKey("ZAnaSp", REAL, "z_analysis_specular"),
    InputKey("doFluxAnalysis", LOGICAL, "do_flux_analysis"),
    InputKey("rFluxAnalysis", REAL, None),
    InputKey("EMinFluxAna", REAL, None),
    InputKey("EMaxFluxAna", REAL, None),
    InputKey("DeltaEFluxAna", REAL, None),
    InputKey("VIn", INTEGER, "v_in"),
    InputKey("JIn", INTEGER, "j_in"),
    InputKey("MJIn", INTEGER, "mj_in"),
    InputKey("Potential", LOGICAL, "only_potential"),
    InputKey("ImposeVMin", LOGICAL, "impose_v_min"),
    InputKey("VMin", REAL, "v_min"),
    InputKey("ImposeVMax", LOGICAL, "impose_v_max"),
    InputKey("VMax", REAL, "v_max"),
    InputKey("UseSwitching", LOGICAL, "use_switching"),
    InputKey("ZminSwitch", REAL, "z_min_switch"),
    InputKey("ZmaxSwitch", REAL, "z_max_switch"),
    InputKey("UseSymmetry", LOGICAL, "use_symmetry"),
    InputKey("nPhisOptimal", LOGICAL, "phis_optimal"),
    InputKey("SymmetryTolerance", REAL, "symmetry_tolerance"),
    InputKey("CRPlibDataSet", STRING, "crp_dataset"),
    InputKey("DistortedPESMode", INTEGER, None),
    InputKey("SCMParFilepath", STRING, None),
    InputKey("PosFilepath", STRING, None),
    InputKey("PosFileSlabNum", INTEGER, None),
    InputKey("EnableSurfExpan", LOGICAL, None),
    InputKey("SCMCutoff", REAL, None),
    InputKey("ZOptEnergy", REAL, "z_opt_energy"),
    InputKey("ZOptFactor", REAL, "z_opt_factor"),
    InputKey("ROptEnergy", REAL, "r_opt_energy"),
    InputKey("ROptFactor", REAL, "r_opt_factor"),
    InputKey("ZSpOptEnergy", REAL, "zsp_opt_energy"),
    InputKey("ZSpOptFactor", REAL, "zsp_opt_factor"),
    InputKey("ZOptLowStart", REAL, "z_opt_low_start"),
    InputKey("ZOptLowEnd", REAL, "z_opt_low_end"),
    InputKey("ZOptLowEnergy", REAL, "z_opt_low_energy"),
    InputKey("ZOptLowFactor", REAL, "z_opt_low_factor"),
    InputKey("ZSpOptLowStart", REAL, "zsp_opt_low_start"),
    InputKey("ZSpOptLowEnd", REAL, "zsp_opt_low_end"),
    InputKey("ZSpOptLowEnergy", REAL, "zsp_opt_low_energy"),
    InputKey("ZSpOptLowFactor", REAL, "zsp_opt_low_factor"),
    InputKey("RestartDir", STRING, "restart_dir"),
    InputKey("RestartFormat", STRING, "restart_format"),
    InputKey("CheckpointFormat", STRING, "checkpoint_format"),
    InputKey("writeDensities", LOGICAL, None),
    InputKey("writePotential", LOGICAL, None),
    InputKey("PlotSites", SITE_LIST, None),
    InputKey("DoDump", LOGICAL, None),
    InputKey("FullSMatrix", LOGICAL, None),
)

# Fortran features this program does not have (reject_unsupported_keys).
UNSUPPORTED_LOGICAL_KEYS = ("SystemMix", "UseRCutOff", "write2DPES")  # accepted when false
UNSUPPORTED_KEYS = (
    "SystemID2",
    "SystemMixCoeff",
    "RCutOff",
    "x2DPES",
    "y2DPES",
    "theta2DPES",
    "phi2DPES",
    "PESMixingAlfa",
)  # rejected when present

# The absorbers that always need an energy or a factor.
ABSORBER_KEY_PAIRS = (("ZOptEnergy", "ZOptFactor"), ("ROptEnergy", "ROptFactor"), ("ZSpOptEnergy", "ZSpOptFactor"))


@dataclass(frozen=True)
class ParsedInput:
    """The result of reading: the parameters (None if the program would stop while reading) and the findings."""

    parameters: InputParameters | None
    findings: list[Finding]


def all_known_keys() -> tuple[InputKey, ...]:
    """Every key the program accepts, incl. the unsupported logical ones (accepted when false)."""
    unsupported = tuple(InputKey(name, LOGICAL, None) for name in UNSUPPORTED_LOGICAL_KEYS)
    return unsupported + REQUIRED_KEYS + OPTIONAL_KEYS


def parse_input(values: dict) -> ParsedInput:
    """Reads a JSON object (as from json.loads) like parse_input in json_input.cpp."""
    reading = FindingList("parse_input")
    reject_unsupported_keys(values, reading)
    fields = read_keys(values, reading)
    report_unknown_keys(values, reading)
    warn_about_plot_data(values, reading)
    if has_errors(reading.findings):
        return ParsedInput(None, reading.findings)

    parameters = InputParameters(**fields)
    check_while_reading(values, parameters, reading)
    if parameters.mj_max < 0:  # a negative MJMax means JMax
        parameters = replace(parameters, mj_max=parameters.j_max)
    return ParsedInput(None if has_errors(reading.findings) else parameters, reading.findings)


def reject_unsupported_keys(values: dict, reading: FindingList) -> None:
    """Fortran features that are not ported: a logical one set to true, or any of the others present."""
    for name in UNSUPPORTED_LOGICAL_KEYS:
        reading.error_if(values.get(name) is True, name, f"{name} = true is not supported by this program")
    for name in UNSUPPORTED_KEYS:
        reading.error_if(name in values, name, f"{name} is not supported by this program")


def read_keys(values: dict, reading: FindingList) -> dict:
    """Type-checks every known key that is present; returns the InputParameters fields that were read."""
    fields = {}
    for key in REQUIRED_KEYS + OPTIONAL_KEYS:
        if key.name not in values:
            reading.error_if(key in REQUIRED_KEYS, key.name, "required input key is missing")
            continue
        value, problem = typed_value(key.kind, values[key.name])
        if problem:
            reading.error(key.name, problem)
        elif key.field is not None:
            fields[key.field] = value
    for name in UNSUPPORTED_LOGICAL_KEYS:
        if name in values:
            _, problem = typed_value(LOGICAL, values[name])
            reading.error_if(problem is not None, name, problem or "")
    return fields


def typed_value(kind: str, value: object) -> tuple[object, str | None]:
    """The value as the program reads it, or a problem message (JsonReader::as_real, as_integer, …).

    Note: Python's bool is a kind of int, so booleans are excluded explicitly from the numbers.
    """
    is_bool = isinstance(value, bool)
    if kind == REAL:
        if is_bool or not isinstance(value, (int, float)):
            return None, "expected a real number"
        return float(value), None
    if kind == INTEGER:
        if is_bool or not isinstance(value, int):
            return None, "expected an integer"
        if not INT_MIN <= value <= INT_MAX:
            return None, "integer out of range"
        return value, None
    if kind == LOGICAL:
        return (value, None) if is_bool else (None, "expected true or false")
    if kind == SITE_LIST:
        return value, site_list_problem(value)
    return (value, None) if isinstance(value, str) else (None, "expected a string")


def site_list_problem(value: object) -> str | None:
    """What is wrong with a PlotSites value, as JsonReader::read for sites says it, or None."""
    if not isinstance(value, list):
        return "expected a list of sites"
    for site in value:
        if not isinstance(site, dict):
            return "each site must be an object with a name, x and y"
        unknown = [field for field in site if field not in SITE_FIELDS]
        if unknown:
            return f"unknown site field '{unknown[0]}'"
        missing = [field for field in SITE_FIELDS if field not in site]
        if missing:
            return f"a site without '{missing[0]}'"
        if not isinstance(site["name"], str):
            return "expected a string"
        if any(isinstance(site[c], bool) or not isinstance(site[c], (int, float)) for c in ("x", "y")):
            return "expected a real number"
    return None


def warn_about_plot_data(values: dict, reading: FindingList) -> None:
    """The plot-data keys are accepted but not modelled: their memory and disk use are not in the estimates."""
    for name in PLOT_DATA_KEYS:
        if values.get(name) not in (None, False, []):
            reading.warning(name, "plot data are not modelled by the designer: its estimates leave them out")


def report_unknown_keys(values: dict, reading: FindingList) -> None:
    """Every key that is not an input key is an error (check_no_unknown_keys)."""
    known = {key.name for key in all_known_keys()} | set(UNSUPPORTED_KEYS)
    for name in values:
        reading.error_if(name not in known, name, "unknown input key")


def check_while_reading(values: dict, p: InputParameters, reading: FindingList) -> None:
    """The checks the Fortran reader makes (check_while_reading), in its order; SCM keys excluded."""
    reading.error_if(p.projection_operator == 0, "ProjOp", "ProjOp = 0 is not supported anymore.")
    reading.error_if(
        p.projection_operator == 1, "ProjOp", "ProjOp = 1 is a legacy option and not supported; use ProjOp = 2."
    )

    crp = seq(p.system_id, "CRPDistrtSurface") or seq(p.system_id, "CRPlibrary")
    reading.error_if(
        crp and p.crp_dataset == "",
        "CRPlibDataSet",
        f"SystemID {p.system_id} needs CRPlibDataSet (the directory of the CRP data set)",
    )
    reading.error_if(
        not crp and "CRPlibDataSet" in values,
        "CRPlibDataSet",
        "CRPlibDataSet is only used with SystemID CRPlibrary or CRPDistrtSurface",
    )

    if p.do_flux_analysis:
        for name in ("rFluxAnalysis", "EMinFluxAna", "EMaxFluxAna", "DeltaEFluxAna"):
            reading.error_if(name not in values, name, f"{name} must be present")

    for energy_key, factor_key in ABSORBER_KEY_PAIRS:
        require_energy_or_factor(values, energy_key, factor_key, reading)
    for energy_key, energy in (
        ("ZOptEnergy", p.z_opt_energy),
        ("ROptEnergy", p.r_opt_energy),
        ("ZSpOptEnergy", p.zsp_opt_energy),
    ):
        reading.error_if(energy is not None and rlt(energy, 0.0), energy_key, f"{energy_key} < 0.")
    require_low_z_absorber(
        values, p, "ZOptLow", (p.z_opt_low_start, p.z_opt_low_end, p.z_opt_low_energy, p.z_opt_low_factor), reading
    )
    require_low_z_absorber(
        values,
        p,
        "ZSpOptLow",
        (p.zsp_opt_low_start, p.zsp_opt_low_end, p.zsp_opt_low_energy, p.zsp_opt_low_factor),
        reading,
    )

    if p.z_analysis_specular is not None:
        z_ana_sp = p.z_analysis_specular
        reading.error_if(p.dynamical_zsp_analysis, "ZAnaSp", "ZAnaSp can only be used with DynamZSpAna set to FALSE")
        reading.error_if(z_ana_sp > p.zsp_opt_start, "ZAnaSp", "ZAnaSp was set > ZSpOptStart")
        reading.error_if(z_ana_sp < p.z0, "ZAnaSp", "ZAnaSp was set < Z0")


def require_energy_or_factor(values: dict, energy_key: str, factor_key: str, reading: FindingList) -> None:
    """An absorber needs its energy or its factor, not both."""
    has_energy, has_factor = energy_key in values, factor_key in values
    reading.error_if(not has_energy and not has_factor, energy_key, f"{factor_key} or {energy_key} should be present.")
    reading.error_if(has_energy and has_factor, energy_key, f"give either {energy_key} or {factor_key}, not both")


def require_low_z_absorber(
    values: dict,
    p: InputParameters,
    prefix: str,
    absorber: tuple[float | None, float | None, float | None, float | None],
    reading: FindingList,
) -> None:
    """A low-Z absorber (keys prefix + Start, End, Energy, Factor): all absent (off), or Start < End inside
    [ZStart, Zana] plus Energy or Factor. The program stops at the first problem.

    The comparisons are plain < and >, as in the Fortran (not the tolerant ones).
    """
    start, end, energy, factor = absorber
    if start is None and end is None and energy is None and factor is None:
        return  # off
    start_key, end_key = prefix + "Start", prefix + "End"
    if not low_z_absorber_enabled(start, end):
        missing = end_key if start is not None else start_key
        reading.error(missing, f"{missing} is missing: a low-Z absorber needs {start_key} and {end_key}")
    elif start > end:
        reading.error(start_key, f"{start_key} > {end_key} ({start} > {end} bohr): the absorber would be empty")
    elif start == end:
        reading.error(start_key, f"{start_key} = {end_key} ({start} bohr): the absorber needs a positive length")
    elif start < p.z_start:
        reading.error(start_key, f"{start_key} was set smaller than ZStart ({start} < {p.z_start} bohr)")
    elif end > p.z_analysis:
        reading.error(end_key, f"{end_key} was set higher than Zana ({end} > {p.z_analysis} bohr)")
    else:
        require_energy_or_factor(values, prefix + "Energy", prefix + "Factor", reading)
        reading.error_if(
            energy is not None and rlt(energy, 0.0), prefix + "Energy", f"{prefix}Energy < 0 ({energy} eV)"
        )


def key_of_field(field_name: str) -> str:
    """The JSON key of an InputParameters field.

    For the front ends: they change a field with dataclasses.replace and save the new value under this key
    (input_file.save_changed_copy).
    """
    for key in REQUIRED_KEYS + OPTIONAL_KEYS:
        if key.field == field_name:
            return key.name
    raise KeyError(field_name)
