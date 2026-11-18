"""Case gas-phase: the rovibrational states and the open channels against the program's (program/gas_phase.py,
dz_gas_phase.py)."""

import re
from pathlib import Path

from dz_gas_phase import POTENTIAL_EV, dz_gas_phase_curve
from evaluation import evaluate
from program.gas_phase import GasPhase
from program.unit_conversion import hartree_to_ev
from settings import Settings
from tests.helpers import DZ_INPUT, DZ_LOG, HARDWALL_INPUT, check, check_equal, load_parameters

# E(v = 0, j = 0 .. 19) and E(v = 1, j = 0) in eV, from the block GAS-PHASE ROVIBRATIONAL STATES of spodvr for
# the DZ input (provenance: helpers.DZ_INPUT).
DZ_ENERGIES_V0 = (
    -4.4770373027650265e00,
    -4.4622216966564530e00,
    -4.4327316865222528e00,
    -4.3888444520128846e00,
    -4.3309632202508830e00,
    -4.2596035638724565e00,
    -4.1753771017120860e00,
    -4.0789738292968822e00,
    -3.9711442116327254e00,
    -3.8526819990945214e00,
    -3.7244085371063074e00,
    -3.5871591435872943e00,
    -3.4417718918911562e00,
    -3.2890788345507338e00,
    -3.1298993818567746e00,
    -2.9650353428355398e00,
    -2.7952671751101659e00,
    -2.6213512617429724e00,
    -2.4440183311510082e00,
    -2.2639732002376021e00,
)
DZ_ENERGY_V1_J0 = -3.9633059644804050e00
# numpy's LAPACK against the program's oneMKL: 6.7e-15 eV at most over all 58 printed states (measured
# 2026-11-18); the tolerance leaves room for another LAPACK.
ENERGY_TOLERANCE = 1e-13  # eV


# psi(r_first) and psi(r_last) of four states, printed with 9 digits in the same block; (1, 2) has both
# negative, which pins the sign rule. numpy's eigenvectors differ from oneMKL's by up to 9e-17 more than half a
# printed unit there (measured 2026-11-18), hence the floor.
DZ_EDGE_VALUES = {
    (0, 0): ("9.50283027e-06", "4.85812633e-05"),
    (0, 3): ("-2.35260083e-06", "5.34211609e-05"),
    (1, 0): ("1.20029145e-04", "4.45749959e-05"),
    (1, 2): ("-8.85940983e-05", "-6.88085376e-05"),
}
EDGE_FLOOR = 1e-15


def gas_phase_of(**changes: object) -> GasPhase:
    """The gas phase of the DZ input with `changes` to InputParameters fields."""
    gas = evaluate(load_parameters(DZ_INPUT, **changes), Settings()).gas_phase
    check(gas is not None, "no gas phase")
    return gas


def test_dz_states_match_the_program(work: Path) -> None:
    """The FGH on the DZ curve gives the program's energies for v = 0, j = 0 .. 19 and v = 1, j = 0."""
    gas = gas_phase_of()
    for j, expected in enumerate(DZ_ENERGIES_V0):
        energy = float(hartree_to_ev(gas.energies[0, j]))
        check(abs(energy - expected) <= ENERGY_TOLERANCE, f"E(0, {j}) = {energy!r}, program {expected!r}")
    energy = float(hartree_to_ev(gas.energies[1, 0]))
    check(abs(energy - DZ_ENERGY_V1_J0) <= ENERGY_TOLERANCE, f"E(1, 0) = {energy!r}, program {DZ_ENERGY_V1_J0!r}")


def test_hard_coded_curve_is_the_logged_one(work: Path) -> None:
    """dz_gas_phase.py holds the 40 values of the block GAS-PHASE POTENTIAL of the DZ run's output, exactly."""
    log = DZ_LOG.read_text()
    block = log[log.index("GAS-PHASE POTENTIAL") : log.index("GAS-PHASE ROVIBRATIONAL STATES")]
    logged = [float(v) for v in re.findall(r"^\s+\d+\s+\S+\s+(\S+)\s*$", block, re.MULTILINE)]
    check_equal(list(POTENTIAL_EV), logged, "V(r) in eV")


def test_dz_edge_values_of_the_wave_functions(work: Path) -> None:
    """psi at the first and last r point of four states, signs included, to half a printed unit plus EDGE_FLOOR."""
    gas = gas_phase_of()
    for (v, j), printed in DZ_EDGE_VALUES.items():
        for value, text in zip((gas.wave_functions[v, j][0], gas.wave_functions[v, j][-1]), printed, strict=True):
            half_unit = 0.5e-8 * 10 ** int(text.split("e")[1])
            check(abs(value - float(text)) <= half_unit + EDGE_FLOOR, f"psi of ({v}, {j}): {value!r}, printed {text}")


def test_dz_open_channels(work: Path) -> None:
    """As the program prints for the DZ input: 2 open vibrational levels, highest open j 10, both analysed."""
    gas = gas_phase_of()
    check_equal((gas.v_max_open, gas.j_max_open, gas.n_v_analysis, gas.j_max_analysis), (2, 10, 2, 10), "channels")
    check_equal(gas.warnings, (), "warnings")


def test_available_energy_includes_the_initial_state(work: Path) -> None:
    """The open channels count from EMaxWP plus the initial state's internal and parallel energy. By hand from
    the program's printed E(0, j): VIn 1 adds E(1, 0) - E(0, 0) = 0.514 eV, and 1.347 <= 0.95 + 0.514 < 1.682
    (j = 14, 16), so j 14 is the highest open; EParallel 0.3 gives 1.035 <= 1.25 < 1.347, so j 12."""
    check_equal(gas_phase_of(v_in=1).j_max_open, 14, "highest open j with VIn 1")
    check_equal(gas_phase_of(e_parallel=0.3).j_max_open, 12, "highest open j with EParallel 0.3")


def test_too_small_vana_and_jana_are_reported(work: Path) -> None:
    """Channels within the 95% interval (up to 0.95 eV for DZ) but not analysed give the program's warnings."""
    check_equal(gas_phase_of(j_analysis=4).warnings, ("Not all open j channels analysed!",), "JAna 4")
    check_equal(gas_phase_of(v_analysis=1).warnings, ("Not all open v channels analysed!",), "VAna 1")
    d = evaluate(load_parameters(DZ_INPUT), Settings()).derived
    check_equal(f"{d.confidence_e_max:.6f}", "0.950000", "95% interval (INPUT block: 0.400000 .. 0.950000 eV)")


def test_states_are_indicative_off_the_curve(work: Path) -> None:
    """The states are the program's only for the curve's own potential and r grid: DZ yes; DZ with nR 48 or the
    hard wall no."""
    check_equal(evaluate(load_parameters(DZ_INPUT), Settings()).gas_phase_indicative, False, "DZ")
    check_equal(evaluate(load_parameters(DZ_INPUT, nr=48), Settings()).gas_phase_indicative, True, "nR 48")
    check_equal(evaluate(load_parameters(HARDWALL_INPUT), Settings()).gas_phase_indicative, True, "hard wall")


def test_analysed_levels_follow_vana_and_jana(work: Path) -> None:
    """VAna >= 1 and JAna >= 0 replace the open levels as the analysed ones; VAna 0 means the open levels."""
    gas = gas_phase_of(v_analysis=1, j_analysis=4)
    check_equal((gas.n_v_analysis, gas.j_max_analysis), (1, 4), "VAna 1, JAna 4")
    check_equal(gas_phase_of(v_analysis=0).n_v_analysis, 2, "VAna 0")


def test_small_rotational_basis_is_reported(work: Path) -> None:
    """With JMax 6 (and MJMax 6) the open j reach beyond the basis: the program's warning."""
    gas = gas_phase_of(j_max=6, mj_max=6)
    check("Rotational basis too small to describe all open channels!" in gas.warnings, f"warnings {gas.warnings}")


def test_wave_functions_are_normalised(work: Path) -> None:
    """Each state is normalised to sum psi^2 dr = 1 on the curve's grid."""
    gas = gas_phase_of()
    dr = dz_gas_phase_curve().dr
    for v, j in ((0, 0), (1, 0), (0, 10)):
        norm = float((gas.wave_functions[v, j] ** 2).sum() * dr)
        check(abs(norm - 1.0) <= 1e-12, f"norm of ({v}, {j}): {norm!r}")


TESTS = [
    test_dz_states_match_the_program,
    test_hard_coded_curve_is_the_logged_one,
    test_dz_edge_values_of_the_wave_functions,
    test_dz_open_channels,
    test_available_energy_includes_the_initial_state,
    test_too_small_vana_and_jana_are_reported,
    test_states_are_indicative_off_the_curve,
    test_analysed_levels_follow_vana_and_jana,
    test_small_rotational_basis_is_reported,
    test_wave_functions_are_normalised,
]
