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
    """A complete, converged, unitary run with packet weight everywhere: only the information on the norm left."""
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
            d.r_grid_ends(1e-9, (0, 0)),
            d.still_moving(p, p, energies),
            d.norm_left(1e-5),
        ]
    )
    assert [f.name for f in findings] == ["norm left"]
    assert d.slowest_converged(0.03, np.zeros(5), energies).severity == "information"


def test_each_finding_fires_just_above_its_limit_and_for_nan() -> None:
    """Each threshold: just above it the finding applies, at it not; NaN always applies."""
    assert d.unanalysed_flux(1.01 * d.UNANALYSED_LIMIT) and not d.unanalysed_flux(d.UNANALYSED_LIMIT)
    assert d.unanalysed_flux(NAN)
    assert d.unitarity(1.01 * d.UNITARITY_LIMIT, 0.3) and not d.unitarity(d.UNITARITY_LIMIT, 0.3)
    assert d.unitarity(NAN, 0.3)
    assert d.r_grid_ends(1.01 * d.R_EDGE_LIMIT, (0, 0)) and not d.r_grid_ends(d.R_EDGE_LIMIT, (0, 0))
    assert d.r_grid_ends(NAN, (0, 0))
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
    """Every section a finding cites is a '## ' heading of inputs/README.md."""
    headings = set(re.findall(r"^## (.+)$", README.read_text(encoding="utf-8"), re.MULTILINE))
    for section in (d.SECTION_TIME, d.SECTION_STEPS):
        assert section in headings, section
