"""Tests (ctest WavePacketDesigner.InputReading): the key table, the value types and the read-time checks
(program/json_input.py)."""

import dataclasses
import re

from program.config import InputParameters
from program.findings import Severity
from program.json_input import OPTIONAL_KEYS, REQUIRED_KEYS, UNSUPPORTED_KEYS, UNSUPPORTED_LOGICAL_KEYS, parse_input
from tests.helpers import (
    DZ_INPUT,
    ER_INPUT,
    HARDWALL_INPUT,
    REPOSITORY,
    check,
    check_equal,
    check_has_error,
    load_values,
)

CPP_READER = REPOSITORY / "src/input/json_input.cpp"


def read_with(changes: dict, removed: tuple[str, ...] = ()) -> object:
    """parse_input of the DZ input with keys changed or removed."""
    values = load_values(DZ_INPUT) | changes
    for key in removed:
        values.pop(key)
    return parse_input(values)


def test_reference_inputs_are_readable() -> None:
    """The DZ, ER and hard-wall inputs read without findings."""
    for path in (DZ_INPUT, ER_INPUT, HARDWALL_INPUT):
        parsed = parse_input(load_values(path))
        check(parsed.parameters is not None and parsed.findings == [], f"{path}: {parsed.findings}")


def test_missing_required_key() -> None:
    """A missing required key stops the reading; every missing key is named."""
    parsed = read_with({}, removed=("TolSp", "ZStart"))
    check(parsed.parameters is None, "the reading should fail")
    check_has_error(parsed.findings, "TolSp", "required input key is missing", "TolSp")
    check_has_error(parsed.findings, "ZStart", "required input key is missing", "ZStart")


def test_value_types() -> None:
    """Integers must be JSON integers (160.0 is not), booleans are not numbers, a real accepts an integer."""
    check_has_error(read_with({"nZ": 160.0}).findings, "nZ", "expected an integer", "nZ = 160.0")
    check_has_error(read_with({"nZ": True}).findings, "nZ", "expected an integer", "nZ = true")
    check_has_error(read_with({"dZ": True}).findings, "dZ", "expected a real number", "dZ = true")
    check_has_error(read_with({"UseSymmetry": 1}).findings, "UseSymmetry", "expected true or false", "logical")
    check_has_error(read_with({"SystemID": 3}).findings, "SystemID", "expected a string", "string")
    check_has_error(read_with({"nZ": 2**31}).findings, "nZ", "integer out of range", "nZ = 2^31")
    parsed = read_with({"Z0": 18})
    check(parsed.parameters is not None and isinstance(parsed.parameters.z0, float), "Z0 = 18 is read as a real")


def test_unknown_and_unsupported_keys() -> None:
    """Unknown keys are errors; unsupported Fortran features are refused (logical ones only when true)."""
    check_has_error(read_with({"dz": 0.12}).findings, "dz", "unknown input key", "lower-case dz")
    check_has_error(read_with({"SystemMix": True}).findings, "SystemMix", "not supported", "SystemMix true")
    check(read_with({"SystemMix": False}).parameters is not None, "SystemMix = false is accepted")
    check_has_error(read_with({"RCutOff": 3.0}).findings, "RCutOff", "not supported", "RCutOff present")


def test_absorber_needs_energy_or_factor() -> None:
    """Each absorber needs its …Energy or its …Factor, not both; a negative energy is an error."""
    check_has_error(read_with({"ZOptFactor": 0.01}).findings, "ZOptEnergy", "not both", "both")
    check_has_error(read_with({}, removed=("ROptEnergy",)).findings, "ROptEnergy", "should be present", "neither")
    check(read_with({"ZOptFactor": 0.01}, removed=("ZOptEnergy",)).parameters is not None, "factor only")
    check_has_error(read_with({"ZSpOptEnergy": -0.1}).findings, "ZSpOptEnergy", "ZSpOptEnergy < 0.", "negative")


def test_low_z_absorber_rules() -> None:
    """A low-Z absorber: Start and End both, Start < End, inside [ZStart, Zana], Energy or Factor; plain < and >."""
    cases = (
        ({"ZOptLowStart": 0.0}, "ZOptLowEnd", "is missing"),
        ({"ZOptLowStart": 2.0, "ZOptLowEnd": 1.0, "ZOptLowEnergy": 0.1}, "ZOptLowStart", "would be empty"),
        ({"ZOptLowStart": 1.0, "ZOptLowEnd": 1.0, "ZOptLowEnergy": 0.1}, "ZOptLowStart", "positive length"),
        ({"ZOptLowStart": -1.0, "ZOptLowEnd": 1.0, "ZOptLowEnergy": 0.1}, "ZOptLowStart", "smaller than ZStart"),
        ({"ZSpOptLowStart": 0.0, "ZSpOptLowEnd": 16.0, "ZSpOptLowEnergy": 0.1}, "ZSpOptLowEnd", "higher than Zana"),
        ({"ZOptLowStart": 0.0, "ZOptLowEnd": 2.0}, "ZOptLowEnergy", "should be present"),
        ({"ZOptLowStart": 0.0, "ZOptLowEnd": 2.0, "ZOptLowEnergy": -0.2}, "ZOptLowEnergy", "< 0"),
    )
    for changes, key, message in cases:
        check_has_error(read_with(changes).findings, key, message, str(changes))
    good = read_with({"ZOptLowStart": 0.0, "ZOptLowEnd": 2.0, "ZOptLowEnergy": 0.1})
    check(good.parameters is not None and good.findings == [], f"a valid low-Z absorber: {good.findings}")


def test_specular_analysis_line_rules() -> None:
    """ZAnaSp only with DynamZSpAna = false, and between Z0 and ZSpOptStart."""
    check_has_error(read_with({"ZAnaSp": 20.0}).findings, "ZAnaSp", "DynamZSpAna", "dynamic line")
    fixed = {"DynamZSpAna": False}
    check_has_error(read_with(fixed | {"ZAnaSp": 23.0}).findings, "ZAnaSp", "> ZSpOptStart", "beyond the absorber")
    check_has_error(read_with(fixed | {"ZAnaSp": 17.0}).findings, "ZAnaSp", "< Z0", "below Z0")
    check(read_with(fixed | {"ZAnaSp": 21.0}).parameters is not None, "ZAnaSp = 21 is valid")


def test_negative_mj_max_means_j_max() -> None:
    """MJMax absent or negative becomes JMax while reading (0 stays 0 here; the basis treats it as JMax)."""
    check_equal(read_with({}, removed=("MJMax",)).parameters.mj_max, 28, "MJMax absent")
    check_equal(read_with({"MJMax": -1}).parameters.mj_max, 28, "MJMax = -1")
    check_equal(read_with({"MJMax": 0}).parameters.mj_max, 0, "MJMax = 0")


def test_key_table_matches_the_cpp_reader() -> None:
    """The key table lists exactly the keys json_input.cpp reads, the required ones in the same order."""
    source = CPP_READER.read_text()
    required_body = source[source.index("void read_required_keys") : source.index("struct LowZAbsorberKeys")]
    cpp_required = re.findall(r'in\.read\("([^"]+)"', required_body)
    check_equal([k.name for k in REQUIRED_KEYS], cpp_required, "required keys (in order)")

    cpp_read = set(re.findall(r'in\.read\("([^"]+)"', source))
    low_z = {k.name for k in OPTIONAL_KEYS if k.name.startswith(("ZOptLow", "ZSpOptLow"))}
    python_optional = {k.name for k in OPTIONAL_KEYS} - low_z
    check_equal(sorted(python_optional), sorted(cpp_read - set(cpp_required)), "optional keys")
    for name in sorted(low_z) + list(UNSUPPORTED_LOGICAL_KEYS) + list(UNSUPPORTED_KEYS):
        check(f'"{name}"' in source, f"{name} is not a key of json_input.cpp")


def test_every_parameter_field_has_a_key() -> None:
    """Each InputParameters field is filled from exactly one key of the table."""
    table_fields = [k.field for k in REQUIRED_KEYS + OPTIONAL_KEYS if k.field is not None]
    check_equal(len(table_fields), len(set(table_fields)), "fields used twice")
    check_equal(sorted(table_fields), sorted(f.name for f in dataclasses.fields(InputParameters)), "fields")


def test_plot_sites_are_read_like_the_program() -> None:
    """PlotSites: a list of objects with a name, x and y (json_input.cpp's site reader, the same messages); a
    plot-data key (each of PLOT_DATA_KEYS) gives a warning (the designer does not model plot data) but no error."""
    sites = [{"name": "top", "x": 0.0, "y": 0.0}]
    read = read_with({"writeDensities": True, "writePotential": True, "PlotSites": sites, "dtFrames": 100.0})
    check(read.parameters is not None, "a valid site list is read")
    warned = {f.key for f in read.findings if f.severity == Severity.WARNING}
    check_equal(sorted(warned), ["PlotSites", "dtFrames", "writeDensities", "writePotential"], "plot-data warnings")
    wrong = {
        "expected a list of sites": 3,
        "each site must be an object": [1],
        "unknown site field 'z'": [{"name": "a", "x": 0.0, "y": 0.0, "z": 1.0}],
        "a site without 'y'": [{"name": "a", "x": 0.0}],
        "expected a string": [{"name": 1, "x": 0.0, "y": 0.0}],
        "expected a real number": [{"name": "a", "x": True, "y": 0.0}],
    }
    for message, value in wrong.items():
        check_has_error(read_with({"writePotential": True, "PlotSites": value}).findings, "PlotSites", message, message)
