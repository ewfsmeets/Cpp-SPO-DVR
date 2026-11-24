"""Tests (ctest WavePacketDesigner.Validate): the checks of validate_input against the C++, and each with a violating
input (program/validate.py).

`test_checks_match_the_cpp` reads the (key, message) of every check from validate.py and from
src/input/validate.cpp and compares the two lists, so a check added to or changed in the C++ fails here.
VIOLATIONS holds one input change per check, starting from the DZ input (valid);
`test_every_check_has_a_violating_case` fails for a check of validate.py without a row.
"""

import ast
import re

from program.potential_systems import ANALYTIC_SYSTEMS, SystemParameters
from program.unit_conversion import angstrom_to_bohr
from program.validate import validate_input
from tests.helpers import DZ_INPUT, REPOSITORY, TOOL_DIR, check, check_equal, check_has_error, load_parameters

DZ_SYSTEM = ANALYTIC_SYSTEMS["H2Cu111DaiLight"]  # Z_asymptotic 15.0, C6v
# The keys of the flux-analysis block of validate.cpp, which the designer does not copy (not ported).
FLUX_KEYS = {"rFluxAnalysis", "EMinFluxAna", "EMaxFluxAna", "DeltaEFluxAna"}

# (InputParameters changes, expected key, part of the expected message), in the order of validate.py. Where a
# row changes more than one key, the other changes keep the input valid in every other respect that matters
# for the targeted check; the comment says which.
VIOLATIONS = [
    ({"system_id": "H2Cu100"}, "SystemID", "Unknown system ID"),
    ({"z_dof": False}, "ZDoF", "6D"),
    ({"r_dof": False}, "RDoF", "6D"),
    ({"x_dof": False}, "XDoF", "6D"),
    ({"y_dof": False}, "YDoF", "6D"),
    ({"theta_dof": False}, "ThetaDoF", "6D"),
    ({"phi_dof": False}, "PhiDoF", "6D"),
    ({"mass1": 0.0}, "Mass1", "Mass of atom 1 <= 0."),
    ({"mass2": -1.0}, "Mass2", "Mass of atom 2 <= 0."),
    ({"mass2": 2.0141}, "HomoNuc", "masses are different"),
    ({"homonuclear": False}, "HomoNuc", "masses are the same"),
    ({"nz": 0}, "nZ", "nZ < 1."),
    ({"nz": 196}, "nZ", "nZ is not a valid FFT size."),
    ({"dz": 0.0}, "dZ", "dZ <= 0."),
    ({"nz_specular": 0}, "nZSpecular", "nZSpecular < 1."),
    ({"nz_specular": 200}, "nZSpecular", "nZSpecular is not a valid FFT size."),
    ({"r_start": 0.0}, "RStart", "RStart <= 0."),
    ({"nr": 0}, "nR", "nR < 1."),
    ({"nr": 46}, "nR", "nR is not a valid FFT size."),
    ({"dr": 0.0}, "dR", "dR <= 0."),
    ({"nx": 0}, "nX", "nX < 1."),
    ({"ny": 0}, "nY", "nY < 1."),
    ({"nx": 13}, "nX", "nX is not a valid FFT size."),
    ({"ny": 13}, "nY", "nY is not a valid FFT size."),
    ({"nx": 14}, "nY", "nX /= nY."),
    ({"j_in": -1}, "JIn", "JIn < 0."),
    ({"j_max": -1}, "JMax", "JMax < 0."),
    ({"mj_max": 30}, "MJMax", "|MJMax| > JMax."),
    ({"mj_in": 1}, "MJIn", "|MJIn| > JIn."),
    ({"j_in": 30}, "JIn", "JIn > JMax."),
    ({"j_in": 28}, "JIn", "JIn >= JMax."),
    ({"j_in": 4, "mj_in": 4, "mj_max": 2}, "MJIn", "|MJIn| > |MJMax|."),  # JIn 4 keeps |MJIn| <= JIn
    ({"j_max": 27, "mj_max": 27}, "JMax", "homo-nuclear but JMax is even"),  # MJMax 27 keeps |MJMax| <= JMax
    ({"mj_max": 27}, "MJMax", "|MJMax| is odd"),
    ({"projection_operator": 3}, "ProjOp", "ProjOp must be 0, 1 or 2."),
    ({"t_start": -1.0}, "tStart", "tStart < 0."),
    # tEnd moves with tStart so that tEnd = tStart + nStep*dtStep still holds; a restart needs its RestartDir.
    ({"t_start": 10.0, "t_end": 17510.0}, "tStart", "new calculation but tStart /= 0."),
    ({"new_calculation": False, "restart_dir": "previous"}, "tStart", "restart but tStart = 0."),
    ({"new_calculation": False, "t_start": 10.0, "t_end": 17510.0}, "RestartDir", "RestartDir is not given"),
    ({"only_analysis": True}, "RestartDir", "analysis-only run"),
    ({"restart_format": "binary"}, "RestartFormat", "RestartFormat must be"),
    ({"checkpoint_format": "binary"}, "CheckpointFormat", "CheckpointFormat must be"),
    ({"t_end": -1.0}, "tEnd", "tEnd < 0."),
    ({"t_start": 20000.0}, "tStart", "tStart > tEnd."),
    ({"t_end": 17000.0}, "tEnd", "tEnd /= tStart + nStep*dtStep."),
    ({"time_step": 0.0}, "dtStep", "dtStep <= 0."),
    ({"n_steps": -1}, "nStep", "nStep < 0."),
    ({"dt_analysis": 0.0}, "dtAna", "dtAna <= 0."),
    ({"dt_analysis": 1.0}, "dtAna", "dtAna < dtStep."),
    ({"dt_averages": 0.0}, "dtAve", "dtAve <= 0."),
    ({"dt_averages": 1.0}, "dtAve", "dtAve < dtStep."),
    ({"j_analysis": -1}, "JAna", "JAna < 0."),
    ({"v_analysis": -1}, "VAna", "VAna < 0."),
    ({"diffraction_order_analysis": -1}, "DifOrdAna", "DifOrdAna < 0."),
    ({"n_energies": 1}, "NrOfEs", "NrOfEs = 1"),
    ({"n_energies": -95}, "NrOfEs", "NrOfEs = -95"),
    ({"dt_analysis": 30.0}, "dtAna", "Propagation time must be an integer times dtAna."),
    ({"dt_averages": 10.0}, "dtAve", "dtAve < dtAna."),
    ({"dt_averages": 510.0}, "dtAve", "dtAve must be an integer times dtAna."),
    ({"z_analysis": 14.88}, "Zana", "ZAna < ZAsymptotic."),
    # ZOptStart moves with Zana, which must not lie beyond it (that would be three more errors).
    ({"z_analysis": 15.05, "z_opt_start": 15.05}, "Zana", "ZAna not on a grid point."),
    ({"z_start": 15.12}, "Zana", "ZAna < ZStart."),
    # The same for a Zana beyond the grid; ZOptEnd then lies beyond it too, a second error.
    ({"z_analysis": 20.04, "z_opt_start": 20.04, "z_opt_end": 20.04}, "Zana", "too big for the scattering grid"),
    ({"z_opt_start": 14.88}, "Zana", "ZAnalysis > ZOptStart."),
    ({"j_analysis": 30}, "JAna", "JAna > JMax."),
    ({"v_analysis": 40}, "VAna", "VAna > nR-1."),
    ({"diffraction_order_analysis": 23}, "DifOrdAna", "DifOrdAna > nX+nY-2"),
    ({"dt_analysis": 100.0}, "dtAna", "dtAna is too big for fine analysis."),
    ({"e_min_wave_packet": 0.0}, "EMinWP", "EMinWP <= 0."),
    ({"e_max_wave_packet": 0.0}, "EMaxWP", "EMaxWP <= 0."),
    ({"e_min_wave_packet": 0.95}, "EMinWP", "EMinWP >= EMaxWP."),
    ({"z0": 14.0}, "Z0", "Z0 < ZAsymptotic."),
    ({"z0": 23.0}, "Z0", "Z0 > ZSpOptStart."),
    ({"v_in": -1}, "VIn", "VIn < 0."),
    ({"v_in": 40}, "VIn", "VIn > nR-1."),
    ({"j_in": 4, "mj_in": 4, "mj_max": 2}, "MJIn", "MJIn > MJMax."),  # as for |MJIn| > |MJMax| above
    ({"v_in": 1, "v_analysis": 1}, "VIn", "VIn >= VAna."),
    ({"j_in": 2, "j_analysis": 2}, "JIn", "JIn >= JAna."),
    ({"v_in": 39}, "VIn", "VIn = nR-1 = 39"),  # analysis D3 (nR 40)
    ({"impose_v_min": True, "impose_v_max": True, "v_min": 11.0}, "VMin", "VMin > VMax."),
    ({"symmetry_tolerance": 0.0}, "SymmetryTolerance", "SymmetryTolerance <= 0."),
    ({"z_opt_start": 19.5}, "ZOptStart", "ZOptStart > ZOptEnd."),
    ({"z_opt_end": 19.2}, "ZOptEnd", "ZOptEnd is too big for the scattering grid in Z."),
    ({"z_opt_start": 14.88}, "ZOptStart", "ZOptStart < ZAna."),
    ({"z_opt_start": 14.88}, "Zana", "ZAna > ZOptStart."),
    ({"r_opt_start": 9.0}, "ROptStart", "ROptStart > ROptEnd."),
    ({"r_opt_end": 8.5}, "ROptEnd", "ROptEnd is too big for the grid in r."),
    ({"zsp_opt_start": 31.0}, "ZSpOptStart", "ZSpOptStart > ZSpOptEnd."),
    ({"zsp_opt_end": 30.72}, "ZSpOptEnd", "ZSpOptEnd is too big for the specular grid in Z."),
    ({"zsp_opt_start": 14.0}, "ZSpOptStart", "ZSpOptStart < ZAsymptotic."),
    ({"e_parallel": -0.1}, "EParallel", "EParallel < 0."),
    ({"phi_parallel": -1.0}, "PhiParallel", "PhiParallel < 0."),
    ({"phi_parallel": 361.0}, "PhiParallel", "PhiParallel > 360 degrees."),
    ({"use_switching": True, "z_min_switch": 8.0, "z_max_switch": 8.0}, "ZminSwitch", "ZminSwitch >= ZmaxSwitch."),
    ({"write_2d_pes": True, "x_2d_pes": -0.1}, "x2DPES", "x2DPES < 0."),
    ({"write_2d_pes": True, "x_2d_pes": 4.85}, "x2DPES", "x2DPES > alat."),
    ({"write_2d_pes": True, "y_2d_pes": -0.1}, "y2DPES", "y2DPES < 0."),
    ({"write_2d_pes": True, "y_2d_pes": 4.85}, "y2DPES", "y2DPES > blat."),
    ({"write_2d_pes": True, "theta_2d_pes": -1.0}, "theta2DPES", "theta2DPES < 0."),
    ({"write_2d_pes": True, "theta_2d_pes": 361.0}, "theta2DPES", "theta2DPES > 360."),
    ({"write_2d_pes": True, "phi_2d_pes": -1.0}, "phi2DPES", "phi2DPES < 0."),
    ({"write_2d_pes": True, "phi_2d_pes": 361.0}, "phi2DPES", "phi2DPES > 360."),
]


def findings_of(changes: dict) -> list:
    return validate_input(load_parameters(DZ_INPUT, **changes), DZ_SYSTEM, 15.0)


def test_valid_input_has_no_findings() -> None:
    """The DZ input passes every check."""
    check_equal(findings_of({}), [], "findings of the DZ input")


def test_limits_just_inside_are_valid() -> None:
    """The analysis limits (r1726: JIn < JMax, VIn < VAna, JIn < JAna; analysis D3: VIn < nR-1) just inside."""
    for changes in ({"j_in": 26}, {"v_in": 0, "v_analysis": 1}, {"j_in": 0, "j_analysis": 2}, {"v_in": 38}):
        check_equal(findings_of(changes), [], f"findings of {changes}")


def test_each_violation_is_found() -> None:
    """Each row of VIOLATIONS gives its error."""
    for changes, key, message in VIOLATIONS:
        check_has_error(findings_of(changes), key, message, str(changes))


def python_checks() -> list[tuple[str, str]]:
    """(key, message) of every `checks.error_if(violated, key, message)` call of validate.py, in source order.

    An f-string message is returned as a format string with {} for each value, as std::format writes it.
    """
    tree = ast.parse((TOOL_DIR / "program/validate.py").read_text())
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "attr", "") == "error_if"]
    calls.sort(key=lambda call: call.lineno)
    found = []
    for call in calls:
        key, message = call.args[1], call.args[2]
        if isinstance(message, ast.JoinedStr):
            parts = [v.value if isinstance(v, ast.Constant) else "{}" for v in message.values]
            found.append((key.value, "".join(parts)))
        else:
            found.append((key.value, message.value))
    return found


# One `check(violated, "Key", message)` of validate.cpp: the message is one or more string literals,
# possibly inside std::format(...).
CPP_CHECK = re.compile(r'check\(.*?,\s*"(\w+)",\s*(?:std::format\(\s*)?((?:"(?:[^"\\]|\\.)*"\s*)+)', re.DOTALL)
CPP_STRING = re.compile(r'"((?:[^"\\]|\\.)*)"')


def cpp_checks() -> list[tuple[str, str]]:
    """(key, message) of every check of validate_input in src/input/validate.cpp, in source order."""
    source = (REPOSITORY / "src/input/validate.cpp").read_text()
    body = source.split("void validate_input(", 1)[1]  # skip the definition of check() itself
    found = []
    for key, literals in CPP_CHECK.findall(body):
        message = "".join(CPP_STRING.findall(literals)).replace('\\"', '"')
        found.append((key, message))
    return found


def test_checks_match_the_cpp() -> None:
    """validate.py has the checks of validate.cpp with the same keys and messages, in the same order; only the
    flux-analysis block is left out (FLUX_KEYS)."""
    cpp = [(key, message) for key, message in cpp_checks() if key not in FLUX_KEYS]
    python = python_checks()
    check(len(cpp) > 80, f"only {len(cpp)} checks found in validate.cpp: the pattern no longer matches")
    for i, (expected, actual) in enumerate(zip(cpp, python, strict=False)):  # lengths are compared below
        check_equal(actual, expected, f"check {i} of validate.py")
    check_equal(len(python), len(cpp), "number of checks")


def test_every_check_has_a_violating_case() -> None:
    """Every check of validate.py is triggered by some row of VIOLATIONS."""
    found = {(f.key, f.message) for changes, _, _ in VIOLATIONS for f in findings_of(changes)}
    for key, message in python_checks():
        fixed_part = message.split("{}")[0]  # the start of an f-string message, before its first value
        check(any(k == key and m.startswith(fixed_part) for k, m in found), f"no row triggers {key}: '{message}'")


def test_pes_cut_bounds_without_a_lattice() -> None:
    """A system whose lattice the designer does not know (a CRP data set) skips the write2DPES upper bounds of x and
    y; their lower bounds still apply."""
    unknown_lattice = SystemParameters(15.0, "C6v")
    beyond = load_parameters(DZ_INPUT, write_2d_pes=True, x_2d_pes=100.0, y_2d_pes=100.0)
    check_equal(validate_input(beyond, unknown_lattice, 15.0), [], "findings without a lattice")
    negative = load_parameters(DZ_INPUT, write_2d_pes=True, x_2d_pes=-1.0)
    check_has_error(validate_input(negative, unknown_lattice, 15.0), "x2DPES", "x2DPES < 0.", "lower bound")


def test_lattice_constants_match_the_cpp() -> None:
    """The lattice constants of the analytic potentials are the literals of their C++ sources (Persson's in
    angstrom, converted with the legacy factor as persson_cu111_potential.cpp converts it)."""
    cases = [  # (SystemID, C++ source, the text that holds the constants there, (a, b) in bohr)
        ("HardWallPES", "hard_wall_potential.cpp", "parameters_{6.239660, 5.239660,", (6.239660, 5.239660)),
        ("H2Cu111DaiLight", "dai_light_cu111_potential.cpp", "lattice_constant = 4.84;", (4.84, 4.84)),
        ("H2Cu111Persson", "persson_cu111_potential.cpp", "cc = 2.5173;", (angstrom_to_bohr(2.5173),) * 2),
    ]
    for name, file, text, (a, b) in cases:
        check(text in (REPOSITORY / "src/potential" / file).read_text(), f"{text} not in {file}: update this test")
        system = ANALYTIC_SYSTEMS[name]
        check_equal((system.lattice_a, system.lattice_b), (a, b), f"lattice of {name}")
