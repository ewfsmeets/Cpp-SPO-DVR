"""The verdict's findings on synthetic values: each fires just above its limit and not at it, NaN always
fires, the sections exist in inputs/README.md; and E_n,min against the example of inputs/README.md section 8."""

import re
from pathlib import Path

import numpy as np

from compute import diagnostics as d

H2_MASS = (1.0078 + 1.0078) * 1822.89  # electron masses, the program's legacy factor
EV = 27.2114  # the program's legacy hartree -> eV
NAN = float("nan")
README = Path(__file__).resolve().parents[3] / "inputs" / "README.md"


def test_lowest_converged_energy_of_the_readme_example() -> None:
    """DZ: 16·354 + 13·1290 ~ 22 000 a.u. for E_n >= 0.03 eV (Z0 18, Zana 15, Z_turn 2): M/k_min = 354, and
    M/k_n,min = 1290 at 0.03 eV to the README's last digit; inf when the slowest part has not reached Z_turn."""
    k_min = H2_MASS / 354.0
    slowest_ev = k_min * k_min / (2 * H2_MASS) * EV
    m_over_k_n = H2_MASS / np.sqrt(2 * H2_MASS * 0.03 / EV)
    assert abs(m_over_k_n - 1290) < 1
    geometry = {"z0_bohr": 18.0, "z_analysis_bohr": 15.0, "total_mass_me": H2_MASS, "hartree_to_ev": EV}
    e_n_min = d.lowest_converged_normal_energy_ev(
        16 * 354.0 + 13 * m_over_k_n, slowest_normal_energy_ev=slowest_ev, **geometry
    )
    assert abs(e_n_min - 0.03) < 1e-12
    assert d.lowest_converged_normal_energy_ev(1000.0, slowest_normal_energy_ev=slowest_ev, **geometry) == float("inf")


def test_findings_sorted_worst_first() -> None:
    """Errors before warnings before information; findings that do not apply are dropped."""
    findings = d.sorted_findings([d.norm_left(0.1), None, d.unitarity(1e-2, 0.3), d.diverged(100.0)])
    assert [f.severity for f in findings] == ["error", "warning", "information"]


def test_a_clean_run_has_no_error_or_warning() -> None:
    """A complete, converged, unitary run with packet weight everywhere, its absorbers returning little: only
    information (the absorber reflection, the norm left)."""
    energies = np.linspace(0.2, 0.6, 5)
    p = np.full(5, 0.9)
    findings = d.sorted_findings(
        [
            d.stopped_early([], 10000.0, 10000.0),
            d.diverged(None),
            d.unanalysed_flux(1e-6),
            d.unitarity(0.0, 0.3),
            d.top_rotational_level(np.zeros(5), energies, 2),
            d.outer_diffraction_channels(np.zeros(5), energies),
            d.r_grid_ends(1e-9, (0, 0), False),
            d.still_moving(p, p, energies),
            reflection_finding([1e-5] * 5, list(p), energy_index=list(range(5))),
            d.norm_left(1e-5),
        ]
    )
    assert [f.name for f in findings] == ["absorber reflection", "norm left"]
    assert d.slowest_converged(0.03, np.zeros(5), energies).severity == "information"


def test_each_finding_fires_just_above_its_limit_and_for_nan() -> None:
    """Each threshold: just above it the finding applies, at it not; NaN always applies."""
    assert d.unanalysed_flux(1.01 * d.UNANALYSED_LIMIT) and not d.unanalysed_flux(d.UNANALYSED_LIMIT)
    assert d.unanalysed_flux(NAN)
    assert d.unitarity(1.01 * d.UNITARITY_LIMIT, 0.3) and not d.unitarity(d.UNITARITY_LIMIT, 0.3)
    assert d.unitarity(NAN, 0.3)
    assert d.r_grid_ends(1.01 * d.R_EDGE_LIMIT, (0, 0), False) and not d.r_grid_ends(d.R_EDGE_LIMIT, (0, 0), False)
    assert "a smaller RStart or a longer r grid" in d.r_grid_ends(NAN, (0, 0), None).meaning
    inner = d.r_grid_ends(1e-2, (0, 0), True)
    assert "at the first r point" in inner.value and "a smaller RStart" in inner.meaning
    assert "a longer r grid" in d.r_grid_ends(1e-2, (0, 0), False).meaning
    energies = np.array([0.3, 0.4])
    above, at, nan = (np.array([0.0, x]) for x in (1.01 * d.EDGE_LIMIT, d.EDGE_LIMIT, NAN))
    top = d.top_rotational_level
    assert top(above, energies, 2) and not top(at, energies, 2) and top(nan, energies, 2)
    outer = d.outer_diffraction_channels
    assert outer(above, energies) and not outer(at, energies) and outer(nan, energies)
    energies = np.array([0.3, 0.4])
    slow = d.slowest_converged
    assert slow(0.03, np.array([0.0, 1.01 * d.SLOW_CHANNELS_LIMIT]), energies).severity == "warning"
    assert slow(0.03, np.array([0.0, d.SLOW_CHANNELS_LIMIT]), energies).severity == "information"
    assert slow(0.03, np.array([NAN, 0.0]), energies).severity == "warning"
    assert d.diverged(0.0) and not d.diverged(None)


def reflection_finding(
    amplitude: list[float],
    probability: list[float],
    *,
    reflection: list[float] | None = None,
    energy_index: list[int] | None = None,
) -> d.Finding:
    """absorber_reflection for rows given by their returned amplitude a and P (R default a^2: below the
    first-order limit), all at one collision energy unless energy_index says otherwise; the channel of row k is
    (0, k, 0, 0, 0)."""
    n_rows = len(amplitude)
    index = np.zeros(n_rows, dtype=int) if energy_index is None else np.array(energy_index)
    return d.absorber_reflection(
        reflection=np.array(amplitude) ** 2 if reflection is None else np.array(reflection),
        amplitude=np.array(amplitude),
        probability=np.array(probability),
        collision_energy_ev=0.2 + 0.1 * np.arange(index.max() + 1),
        channels=np.array([[0, k, 0, 0, 0] for k in range(n_rows)]),
        energy_index=index,
    )


def test_absorber_reflection_fires_just_above_its_limit_and_for_nan() -> None:
    """One row with P = 1: delta P = 2 a P + 2 sqrt(P a^2 P) = 4 a (exact for a = limit / 4); just above the
    limit a warning, at it information; a NaN in a or P, and an absorber that is off (a = 1), warn."""
    assert reflection_finding([1.01 * d.REFLECTION_LIMIT / 4], [1.0]).severity == "warning"
    at_limit = reflection_finding([d.REFLECTION_LIMIT / 4], [1.0])
    assert at_limit.severity == "information" and at_limit.value.startswith("up to 1.00e-03 in P")
    assert reflection_finding([NAN, 1e-6], [0.5, 0.5], reflection=[0.0, 0.0]).severity == "warning"
    assert reflection_finding([1e-6, 1e-6], [NAN, 0.5]).severity == "warning"
    off = reflection_finding([1.0], [0.5], reflection=[0.0])
    assert off.severity == "warning" and "direct 1.00e+00" in off.value
    assert off.section == d.SECTION_ABSORBERS


def test_absorber_reflection_terms_by_hand() -> None:
    """Three rows at one energy and one at another: Q sums min(a, 1)^2 P per energy (the near-threshold row,
    a = 2 with R = 0.5, capped at 1 but included), the re-scattered term is 2 sqrt(P Q), the direct one 2 a P;
    the near-threshold row is not estimated."""
    amplitude, probability = np.array([0.01, 0.02, 2.0, 0.03]), np.array([0.5, 0.3, 0.1, 0.2])
    error = d.reflection_error(
        reflection=np.array([1e-4, 4e-4, 0.5, 9e-4]),
        amplitude=amplitude,
        probability=probability,
        energy_index=np.array([0, 0, 0, 1]),
        n_energies=2,
    )
    q = np.array([1e-4 * 0.5 + 4e-4 * 0.3 + 1.0 * 0.1, 9e-4 * 0.2])
    assert np.allclose(error.rescattered, 2.0 * np.sqrt(probability * q[[0, 0, 0, 1]]), rtol=1e-15, atol=0.0)
    assert np.allclose(error.direct, 2.0 * amplitude * probability, rtol=1e-15, atol=0.0)
    assert list(error.near_threshold) == [False, False, True, False]
    finding = reflection_finding(list(amplitude), list(probability), reflection=[1e-4, 4e-4, 0.5, 9e-4])
    assert "1 rows with R > 0.1 (near a threshold) not estimated, the largest P among them 1.00e-01" in finding.value


def test_absorber_reflection_near_threshold_fires_just_above_its_limit() -> None:
    """A row above the first-order limit (R 0.5, a = sqrt(R) + T 0.75) is not estimated but can carry an error of
    order its P: with 2 P just above REFLECTION_LIMIT a warning, at it information; the value quotes its P."""
    just_above = reflection_finding([0.75], [1.01 * d.REFLECTION_LIMIT / 2], reflection=[0.5])
    assert just_above.severity == "warning"
    at_limit = reflection_finding([0.75], [d.REFLECTION_LIMIT / 2], reflection=[0.5])
    assert at_limit.severity == "information"
    assert at_limit.value == (
        "no channel estimated; no estimated channel with P >= 0.01; "
        "1 rows with R > 0.1 (near a threshold) not estimated, the largest P among them 5.00e-04"
    )


def test_absorber_reflection_relative_part_only_for_large_p() -> None:
    """A row with P below REFLECTION_RELATIVE_P and a large a (0.05) is left out of the relative error the value
    quotes (2a of the other row), but still counts in the absolute check (it warns when its error is large)."""
    small_p = 0.1 * d.REFLECTION_RELATIVE_P
    finding = reflection_finding([1e-4, 0.05], [0.5, small_p], energy_index=[0, 1])
    assert finding.severity == "information"
    expected = (
        "direct relative error 2a up to 2.00e-04 where P >= 0.01 (E_coll 0.2 eV, (v, j, mJ, n, m) = (0, 0, 0, 0, 0))"
    )
    assert expected in finding.value
    assert "1.00e-01" not in finding.value
    transmitting = reflection_finding([1e-4, 0.5], [0.5, small_p], reflection=[1e-8, 0.01], energy_index=[0, 1])
    assert transmitting.severity == "warning"  # 2 a P + 2 sqrt(P a^2 P) = 2e-3 for a = 0.5 (T large, R 0.01)
    none = reflection_finding([0.05], [small_p])
    assert "no estimated channel with P >= 0.01" in none.value


def moving(p_end: float, change: float) -> d.Finding | None:
    """still_moving for one energy."""
    return d.still_moving(np.array([p_end]), np.array([p_end - change]), np.array([0.3]))


def test_findings_name_the_worst_energy() -> None:
    """The still-moving and the unitarity findings say at which collision energy they are worst."""
    moved = d.still_moving(np.array([0.5, 0.5, 0.5]), np.array([0.5, 0.4, 0.49]), np.array([0.2, 0.3, 0.4]))
    assert "up to 1.00e-01 at E_coll 0.3 eV" in moved.value and moved.value.startswith("2 of 3 energies")
    assert "E_coll 0.32 eV" in d.unitarity(1e-2, 0.32).value


def test_still_moving_at_its_relative_and_its_absolute_limit() -> None:
    """P_s = 0.5: the relative limit binds (5e-4); P_s = 0.01: the floor binds (1e-4); NaN moves."""
    relative = d.MOVING_RELATIVE * 0.5
    assert moving(0.5, 1.01 * relative) and not moving(0.5, 0.99 * relative)
    assert moving(0.01, 1.01 * d.MOVING_FLOOR) and not moving(0.01, 0.99 * d.MOVING_FLOOR)
    assert moving(NAN, 0.0)


def test_stopped_early() -> None:
    """A missing end-of-run file is an error that says how far the BK file got."""
    finding = d.stopped_early(["averages.dat"], 500.0, 10000.0)
    assert finding.severity == "error" and "BK records up to t = 500 a.u." in finding.value
    assert "no BK record" in d.stopped_early(["averages.dat"], None, 10000.0).value


def test_slow_probability_sums_the_channels_below_e_n_min() -> None:
    """Rows of two energies: only the channels leaving with a normal energy below E_n,min = 0.05 eV count."""
    energy_index = np.array([0, 0, 1, 1, 1])
    outgoing_ev = np.array([0.01, 0.30, 0.049, 0.05, 0.2])
    probability = np.array([0.1, 0.8, 0.2, 0.3, 0.4])
    assert list(d.slow_probability(energy_index, outgoing_ev, probability, 3, 0.05)) == [0.1, 0.2, 0.0]


def test_sections_are_headings_of_the_inputs_readme() -> None:
    """Every section a finding cites (every SECTION_* of compute/diagnostics.py) is a '## ' heading of
    inputs/README.md."""
    headings = set(re.findall(r"^## (.+)$", README.read_text(encoding="utf-8"), re.MULTILINE))
    sections = [value for name, value in vars(d).items() if name.startswith("SECTION_")]
    assert len(sections) >= 3
    for section in sections:
        assert section in headings, section


def test_norm_by_region_names_its_time_and_regions() -> None:
    """Information with the record's time, every region's norm and the specular wave inside its analysis line;
    its meaning names what is not yet analysed."""
    finding = d.norm_by_region(9900.0, {"r >= r_c": 1.5e-6, "r < r_c, Z < zasymp (7 bohr)": 2.0e-3}, 3.0e-5)
    assert finding.severity == "information" and finding.name == "norm by region at the end"
    assert finding.value == (
        "t = 9900 a.u.: r >= r_c: 1.5e-06; r < r_c, Z < zasymp (7 bohr): 2.0e-03; of the bands beyond Zana, the "
        "specular wave inside its analysis line: 3.0e-05"
    )
    assert "rest norm" in finding.meaning and "not yet analysed" in finding.meaning


def test_frame_edge() -> None:
    """A basis edge in the frames: none at or below FRAME_EDGE_LIMIT; a finding of the given severity above it
    naming its worst time; one for NaN."""
    times = np.array([0.0, 100.0, 200.0])
    assert d.frame_edge("warning", "top J", np.array([0.0, 1e-7, d.FRAME_EDGE_LIMIT]), times, "m") is None
    finding = d.frame_edge("warning", "top J", np.array([0.0, 3e-6, 1e-6]), times, "m")
    assert finding.severity == "warning" and finding.value == "up to 3.00e-06 of the initial norm (t = 100 a.u.)"
    assert d.frame_edge("warning", "top J", np.array([0.0, np.nan, 0.0]), times, "m") is not None


def test_rotational_and_grid_edges() -> None:
    """A rotational edge is information (a passing population is normal; the S matrix's j' check decides); a grid
    edge is a warning on its rise above the first frame (compute.frames.grid_edge_rises): no rise raises nothing, a
    rise of 3e-6 warns with that value and time, NaN warns."""
    times = np.array([0.0, 100.0, 200.0])
    rotational = d.rotational_edge("top J", np.array([0.0, 3e-6, 0.0]), times)
    assert rotational.severity == "information" and "S matrix" in rotational.meaning
    assert d.grid_edge("first r points", np.array([0.0, 0.0, -1e-6]), times) is None
    grid = d.grid_edge("first r points", np.array([0.0, 3e-6, 0.0]), times)
    assert grid.severity == "warning" and "above the first frame's" in grid.name
    assert grid.value == "up to 3.00e-06 of the initial norm (t = 100 a.u.)"
    assert d.grid_edge("first r points", np.array([0.0, np.nan, 0.0]), times) is not None


def test_momentum_edge() -> None:
    """A momentum grid's edge is a warning on its rise above the first frame (compute.frame_extras.momentum_edges),
    citing section 4 (the grid must hold the momenta): a rise of 3e-6 warns with its time; 1e-6 does not (the limit
    is inclusive)."""
    times = np.array([0.0, 100.0])
    finding = d.momentum_edge("k_r", np.array([0.0, 3e-6]), times)
    assert finding.severity == "warning" and finding.section == d.SECTION_STEPS
    assert finding.name == "packet at the edge of the basis: k_r (momentum grid, above the first frame's)"
    assert finding.value == "up to 3.00e-06 of the initial norm (t = 100 a.u.)" and "aliasing" in finding.meaning
    assert d.momentum_edge("k_r", np.array([0.0, d.FRAME_EDGE_LIMIT]), times) is None
