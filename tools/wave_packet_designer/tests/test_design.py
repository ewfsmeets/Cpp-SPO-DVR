"""Case design: the recommendations of inputs/README.md (design/rules.py), the absorbers' R and T
(design/transmission.py) and the suggested Z grids (design/suggest.py).

FAILING holds one input change per rule that the rule must flag, starting from the DZ input (which passes
every rule); `test_every_rule_has_a_failing_case` fails for a rule without a row.
"""

import math
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

from design.rules import (
    DesignCheck,
    DesignInput,
    analysed_bound_states,
    dz_holds_the_momenta,
    energy_range,
    energy_resolution,
    gaussian_tail,
    lowest_converged_energy,
)
from design.suggest import fit_time
from design.transmission import reflection_transmission
from evaluation import Evaluation, evaluate, suggest
from program.potential_systems import ANALYTIC_SYSTEMS
from program.process import build_derived
from program.unit_conversion import ev_to_hartree
from report import verdict_lines
from settings import Settings
from tests.helpers import DZ_INPUT, ER_INPUT, HARDWALL_INPUT, check, check_equal, load_parameters

DAI_LIGHT_SYSTEM = ANALYTIC_SYSTEMS["H2Cu111DaiLight"]
# (InputParameters changes, the name of the check that must fail).
FAILING = [
    ({"e_min_wave_packet": 0.05, "e_max_wave_packet": 1.0}, "energy range: k0 / sigma_k"),
    ({"dz": 0.15, "e_min_wave_packet": 0.35, "e_max_wave_packet": 1.05}, "dZ: (k_grid - k0) / sigma_k"),
    ({"z_analysis": 15.12, "z_opt_start": 15.12}, "Zana above the first possible point"),
    ({"z0": 17.5}, "Z0 - Zana in sigma_Z"),
    ({"zsp_opt_start": 21.0}, "ZSpOptStart - outer tail of the packet"),
    ({"z_opt_end": 18.0}, "scattering absorber length"),
    ({"z_opt_end": 18.0}, "ZOptEnd on the last grid point"),
    ({"zsp_opt_end": 30.0}, "ZSpOptEnd on the last grid point"),
    ({"r_opt_end": 8.0}, "ROptEnd on the last grid point"),
    ({"r_opt_start": 2.0}, "norm of analysed bound states beyond ROptStart"),
    # The DZ input itself: its tEnd 17 500 is below the 25 000 of the README's recipe (section 2, step 11).
    ({}, "lowest outgoing normal energy converged by tEnd"),
]


def evaluation_of(**changes: object) -> Evaluation:
    """The evaluation of the DZ input with `changes` (the design checks need the grids)."""
    e = evaluate(load_parameters(DZ_INPUT, **changes), Settings())
    check(e.design != (), f"no design checks: {e.stopped_at} {e.findings}")
    return e


def failed_checks(e: Evaluation) -> list[str]:
    return [c.name for c in e.design if c.passed is False]


def test_dz_passes_every_rule_but_its_time(work: Path) -> None:
    """The DZ input, built by the recipe, passes every rule except the time: its tEnd 17 500 converges outgoing
    energies down to 0.06 eV, as inputs/README.md notes (step 11: 25 000; "the DZ input uses 17 500"). With
    tEnd 25 000 it passes that one too."""
    check_equal(failed_checks(evaluation_of()), ["lowest outgoing normal energy converged by tEnd"], "failed checks")
    check_equal(failed_checks(evaluation_of(t_end=25000.0, n_steps=10000)), [], "failed checks with tEnd 25 000")


def test_time_budget_of_the_readme_example(work: Path) -> None:
    """Section 8's DZ example by hand: the slowest analysed part (0.40 eV, M/k_min = 354 a.u. per bohr) runs
    16 bohr in and, reflected elastically, 13 bohr out: 29 x 354 = 10 266 a.u. (the exact M/k_min gives
    10 252); outgoing E_n = 0.03 eV takes 16 x 354 + 13 x 1290 = 22 434 (22 438). With Z_turn = ZStart = 0:
    33 bohr, 11 666."""
    b = evaluation_of().time_budget
    check(abs(b.t_elastic - 29 * 354) < 20 and abs(b.t_needed - (16 * 354 + 13 * 1290)) < 20, f"{b}")
    b = evaluate(load_parameters(DZ_INPUT), Settings(z_turn=0.0)).time_budget
    check_equal(round(b.t_elastic), 11666, "elastic with Z_turn 0")


def test_time_is_fitted_to_whole_averaging_intervals(work: Path) -> None:
    """A requested tEnd becomes the next whole number of dtAve intervals (DZ: 500 a.u., 200 steps of 2.5) so the
    program's time checks hold: 25 000 stays, 25 100 becomes 25 500 with nStep 10 200."""
    p = load_parameters(DZ_INPUT)
    for requested, expected in ((25000.0, (25000.0, 10000)), (25100.0, (25500.0, 10200))):
        fitted = fit_time(p, requested)
        check_equal((fitted["t_end"], fitted["n_steps"]), expected, f"tEnd {requested}")
        check(evaluate(replace(p, **fitted), Settings()).accepted, f"tEnd {requested} accepted")


def test_memory_must_fit(work: Path) -> None:
    """With the memory of a 12 GiB GPU given, the ER input (14.8 GiB on the GPU, no symmetry) fails the GPU
    check, which the program would refuse; DZ with its symmetry counts (8.8 GiB) fits 9 GiB free."""
    gpu = 12 * 2**30
    er = evaluate(load_parameters(ER_INPUT), Settings(gpu_bytes=gpu))
    check("GPU memory: peak" in failed_checks(er), f"ER failed checks {failed_checks(er)}")
    dz = evaluate(load_parameters(DZ_INPUT), Settings(representatives=19, gpu_bytes=9 * 2**30, ram_bytes=16 * 2**30))
    check_equal([c.passed for c in dz.design if c.name.endswith("memory: peak")], [True, True], "DZ fits")


def test_worst_case_memory_is_information_only(work: Path) -> None:
    """DZ uses UseSymmetry: without --representatives the estimate is the worst case (all 144 (x, y) points,
    17.2 GiB on the GPU). It cannot prove that the run does not fit, so with 9 GiB free GPU memory the check is
    information, named as the worst case, and the verdict says so."""
    e = evaluate(load_parameters(DZ_INPUT), Settings(gpu_bytes=9 * 2**30))
    gpu = next(c for c in e.design if c.name.startswith("GPU memory"))
    check(gpu.passed is None and "worst case" in gpu.name, f"{gpu}")
    check(any("WORST CASE" in line for line in verdict_lines(e)), "verdict")


def test_analysed_states_follow_the_angular_grid(work: Path) -> None:
    """The analysed j are those of the program's analysis basis: for a homonuclear molecule with odd JIn only odd
    j (hardwall_small with JIn 1, JMax 5)."""
    e = evaluate(load_parameters(HARDWALL_INPUT, j_max=5, mj_max=5, j_in=1, mj_in=0), Settings())
    js = {j for _, j in analysed_bound_states(e.gas_phase, e.resources.angular)}
    check(js and all(j % 2 == 1 for j in js), f"analysed j {sorted(js)}")


def test_each_rule_flags_its_case(work: Path) -> None:
    """Each row of FAILING fails its rule."""
    for changes, name in FAILING:
        check(name in failed_checks(evaluation_of(**changes)), f"{changes} does not fail '{name}'")


def test_every_rule_has_a_failing_case(work: Path) -> None:
    """Every check with a limit (not information only) has a row in FAILING."""
    names = {c.name for c in evaluation_of().design if c.passed is not None}
    check_equal(names - {name for _, name in FAILING}, set(), "rules without a failing case")


def rule_value(rule: Callable[[DesignInput], DesignCheck], **changes: object) -> float:
    """The value of one rule for the DZ input with `changes`, without the program's checks (some README rows
    put Zana off the grid)."""
    p = load_parameters(DZ_INPUT, **changes)
    base = evaluation_of()
    return rule(
        DesignInput(p, build_derived(p, DAI_LIGHT_SYSTEM), base.grids, base.resources.angular, base.gas_phase)
    ).value


def test_momentum_margins_of_the_readme_table(work: Path) -> None:
    """inputs/README.md section 4: (k_grid - k0)/sigma_k and the wrapped fraction for its four H2 rows, to the
    digits printed there; and section 3: k0/sigma_k = 3.2 for 0.05-1.0 eV."""
    rows = [
        (0.15, 0.35, 1.05, "4.3", "8e-06"),
        (0.15, 0.40, 1.00, "5.0", "3e-07"),
        (0.14, 0.35, 1.05, "5.1", "1.3e-07"),
        (0.12, 0.40, 0.95, "9.2", "1e-20"),
    ]
    for dz, e_min, e_max, n_sigma, wrapped in rows:
        value = rule_value(dz_holds_the_momenta, dz=dz, e_min_wave_packet=e_min, e_max_wave_packet=e_max)
        digits = len(wrapped.split("e")[0].replace(".", ""))
        check_equal((f"{value:.1f}", f"{gaussian_tail(value):.{digits - 1}e}"), (n_sigma, wrapped), f"dZ {dz}")
    value = rule_value(energy_range, e_min_wave_packet=0.05, e_max_wave_packet=1.0)
    check_equal(f"{value:.1f}", "3.2", "k0/sigma_k for 0.05-1.0 eV")


def test_lowest_converged_energy_of_the_readme_example(work: Path) -> None:
    """Section 8: T = (Z0 - Z_turn) M/k_min + (Zana - Z_turn) M/k_n; DZ: 16 * 354 + 13 * 1290 a.u. for
    E_n = 0.03 eV. With exactly that tEnd the rule gives back 0.03 eV."""
    d = evaluation_of().derived
    k_min = d.kz - 2.0 * math.sqrt(d.wave_packet_alpha)
    k_n = math.sqrt(2.0 * d.total_mass * ev_to_hartree(0.03))
    check_equal((round(d.total_mass / k_min), round(d.total_mass / k_n, -1)), (354, 1290), "a.u. per bohr")
    t_end = 16.0 * d.total_mass / k_min + 13.0 * d.total_mass / k_n
    e = evaluation_of()  # the rule alone: another tEnd would also need new time keys
    energy = lowest_converged_energy(
        DesignInput(replace(e.parameters, t_end=t_end), d, e.grids, e.resources.angular, e.gas_phase)
    ).value
    check(abs(energy - 0.03) <= 1e-12, f"lowest converged energy {energy!r}")


def test_free_wave_passes_unchanged(work: Path) -> None:
    """With strength 0 the transfer matrix gives R = 0 and T = 1 (to the rounding of its 2000 slices)."""
    reflection, transmission = reflection_transmission(ev_to_hartree(0.3), 4.0, 0.0, 3674.2)
    check(reflection < 1e-24 and abs(transmission - 1.0) < 1e-12, f"R {reflection!r}, T {transmission!r}")


def test_dz_absorbers_as_the_readme_says(work: Path) -> None:
    """Section 7, for the DZ absorbers: scattering R + T <= 4e-5 for 0.05-1.0 eV and 8e-5 at 1.5 eV; specular
    R + T <= 4e-9 for 0.2-1.0 eV (one printed digit)."""
    e = evaluation_of()
    mass = e.derived.total_mass
    z_opt, zsp_opt = (next(a for a in e.absorbers if a.name == name) for name in ("ZOpt", "ZSpOpt"))

    def r_plus_t(absorber: object, energy_ev: float) -> float:
        r, t = reflection_transmission(ev_to_hartree(energy_ev), absorber.end - absorber.start, absorber.strength, mass)
        return r + t

    scattering = max(r_plus_t(z_opt, energy) for energy in (0.05, 0.1, 0.2, 0.4, 0.7, 1.0))
    specular = max(r_plus_t(zsp_opt, energy) for energy in (0.2, 0.4, 0.7, 1.0))
    printed = (f"{scattering:.0e}", f"{r_plus_t(z_opt, 1.5):.0e}", f"{specular:.0e}")
    check_equal(printed, ("4e-05", "8e-05", "4e-09"), "R + T")


def test_suggestion_for_the_dz_input(work: Path) -> None:
    """The DZ input's dZ 0.12 holds the momenta, so the recipe keeps it and gives the DZ grid sizes (nZ 160,
    nZSpecular 256); with --dz 0.14 (the largest step that holds them) it gives a smaller grid. Both inputs are
    accepted by the program and pass every rule."""
    for dz, expected_dz in ((None, 0.12), (0.14, 0.14)):
        s = suggest(load_parameters(DZ_INPUT), Settings(), dz)
        check_equal(s.changes["dz"], expected_dz, "dZ")
        check(s.evaluation.accepted and s.evaluation.findings == (), f"{dz}: {s.evaluation.findings}")
        check_equal(failed_checks(s.evaluation), [], f"failed checks with dZ {dz}")
    s = suggest(load_parameters(DZ_INPUT), Settings())
    check_equal((s.changes["nz"], s.changes["nz_specular"]), (160, 256), "grid sizes")


def test_suggestion_with_another_grid_start(work: Path) -> None:
    """With ZStart -2 (as hardwall_small) the grid points are ZStart + i dZ: the suggested Zana, absorber ends
    and grid sizes fit that grid, and the program accepts the input."""
    s = suggest(load_parameters(DZ_INPUT, z_start=-2.04), Settings())
    check(s.evaluation.accepted and failed_checks(s.evaluation) == [], f"{s.evaluation.findings}")
    check_equal(s.changes["z_analysis"], 15.0, "Zana (-2.04 + 142 * 0.12)")


def test_suggestion_is_blocked_by_other_errors(work: Path) -> None:
    """An input the program rejects on a key the suggestion keeps (EMinWP 0) gets no suggestion: the error is
    listed instead."""
    s = suggest(load_parameters(DZ_INPUT, e_min_wave_packet=0.0), Settings())
    check_equal((s.changes, [f.key for f in s.blockers]), ({}, ["EMinWP"]), "blocked suggestion")


def test_limits_hold_at_their_boundary(work: Path) -> None:
    """The README's limits are inclusive: EMaxWP = 4 EMinWP (k0/sigma_k = 6 up to rounding, 5.999999999999999 for
    0.25-1.0 eV) passes, and so does the largest dZ that holds the momenta; a step beyond each fails."""
    check(energy_range_passes(0.25, 1.0) and not energy_range_passes(0.25, 1.001), "EMaxWP = 4 EMinWP")
    p = load_parameters(DZ_INPUT)
    d = build_derived(p, DAI_LIGHT_SYSTEM)
    largest_dz = math.pi / (d.kz + 6.0 * math.sqrt(d.wave_packet_alpha))
    base = evaluation_of()
    for dz, passes in ((largest_dz, True), (largest_dz * 1.001, False)):
        check_equal(
            dz_holds_the_momenta(DesignInput(replace(p, dz=dz), d, base.grids, None, None)).passed, passes, f"{dz}"
        )


def energy_range_passes(e_min: float, e_max: float) -> bool:
    """Whether the energy-range rule passes for the DZ input with this energy range."""
    p = load_parameters(DZ_INPUT, e_min_wave_packet=e_min, e_max_wave_packet=e_max)
    base = evaluation_of()
    return energy_range(DesignInput(p, build_derived(p, DAI_LIGHT_SYSTEM), base.grids, None, None)).passed


def test_energy_resolution_of_the_readme(work: Path) -> None:
    """Section 8: the energy resolution 2 pi / T is 6.8 meV at 25 000 a.u."""
    e = evaluation_of()
    c = energy_resolution(DesignInput(replace(e.parameters, t_end=25000.0), e.derived, e.grids, None, None))
    check_equal(f"{c.value * 1000:.1f}", "6.8", "meV")


def test_absorber_responses_use_the_right_mass_and_energies(work: Path) -> None:
    """ROpt is computed with the reduced mass, the Z absorbers with the total mass, each over its energies."""
    e = evaluation_of()
    d = e.derived
    for response in e.absorber_responses:
        a = next(a for a in e.absorbers if a.name == response.name)
        mass = d.reduced_mass if a.name == "ROpt" else d.total_mass
        r, t = reflection_transmission(ev_to_hartree(response.energies[3]), a.end - a.start, a.strength, mass)
        check_equal(response.total[3], r + t, f"R + T of {a.name}")
        # ZSpOpt: the packet's energies; ZOpt: up to EMaxWP (the DZ ground state has no internal energy to give);
        # ROpt: up to the r kinetic energy at pi/dR, 3.654734 eV in the program's GRIDS block.
        expected = {"ZSpOpt": (0.4, 0.95), "ZOpt": (0.05, 0.95), "ROpt": (0.05, 3.654734)}[a.name]
        window = (round(response.energies[0], 6), round(response.energies[-1], 6))
        check_equal(window, expected, f"energies of {a.name}")


def test_suggestion_for_the_er_input(work: Path) -> None:
    """The recipe on the ER input (Z_asymptotic 8): accepted, every rule passed."""
    s = suggest(load_parameters(ER_INPUT), Settings())
    check(s.evaluation.accepted and failed_checks(s.evaluation) == [], f"{s.evaluation.findings}")


TESTS = [
    test_dz_passes_every_rule_but_its_time,
    test_time_budget_of_the_readme_example,
    test_time_is_fitted_to_whole_averaging_intervals,
    test_memory_must_fit,
    test_worst_case_memory_is_information_only,
    test_analysed_states_follow_the_angular_grid,
    test_each_rule_flags_its_case,
    test_every_rule_has_a_failing_case,
    test_momentum_margins_of_the_readme_table,
    test_lowest_converged_energy_of_the_readme_example,
    test_free_wave_passes_unchanged,
    test_dz_absorbers_as_the_readme_says,
    test_suggestion_for_the_dz_input,
    test_suggestion_with_another_grid_start,
    test_suggestion_is_blocked_by_other_errors,
    test_suggestion_for_the_er_input,
    test_limits_hold_at_their_boundary,
    test_energy_resolution_of_the_readme,
    test_absorber_responses_use_the_right_mass_and_energies,
]
