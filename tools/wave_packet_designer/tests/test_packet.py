"""Case packet: the initial wave packet against the program's, and the TolSp check
(program/initial_wave_packet.py)."""

from pathlib import Path

from evaluation import evaluate
from program.potential_systems import ANALYTIC_SYSTEMS
from program.process import build_derived
from settings import Settings
from tests.helpers import DZ_INPUT, check, check_equal, check_has_error, load_parameters


def test_dz_packet_parameters_match_the_program(work: Path) -> None:
    """The block INPUT of spodvr for the DZ input (provenance: helpers.DZ_INPUT): the masses, alpha and kz,
    printed with 17 digits, bitwise (same operations in the same order)."""
    d = build_derived(load_parameters(DZ_INPUT), ANALYTIC_SYSTEMS["H2Cu111DaiLight"])
    check_equal(d.total_mass, 3.6742170840000003e03, "total mass")
    check_equal(d.reduced_mass, 9.1855427100000009e02, "reduced mass")
    check_equal(d.wave_packet_alpha, 1.9767175470965832e00, "alpha")
    check_equal(d.kz, 1.3205179688916404e01, "kz")


def test_dz_packet_matches_the_program(work: Path) -> None:
    """The block ROTATIONAL BASIS AND INITIAL WAVE PACKET of spodvr for the DZ input (provenance:
    helpers.DZ_INPUT): the norms bitwise (same exp and summation order); the |psi(k)| values to FFT rounding
    (numpy here, FFTW or oneMKL there), absolute 1e-15 for the tiny ones at k = 0 and the edge."""
    info = evaluate(load_parameters(DZ_INPUT), Settings()).packet.info
    check_equal((info.i_low, info.i_high), (125, 183), "allowed region")
    check_equal(info.norm, 9.9999999999999978e-01, "norm")
    check_equal(info.norm_inside, 9.9999999999999978e-01, "norm inside")
    check(abs(info.k_space_maximum - 5.3214583713275465e-01) <= 1e-12 * 5.3214583713275465e-01, "|psi(k)| maximum")
    check(abs(info.k_space_at_zero - 1.4080154529207057e-10) <= 1e-15, "|psi(k)| at k = 0")
    check(abs(info.k_space_at_edge - 3.0196068748602284e-10) <= 1e-15, "|psi(k)| at the edge")


def test_tolsp_error_when_the_packet_reaches_the_analysis_line(work: Path) -> None:
    """Z0 at Zana + 5 sigma_Z (16.78): about 3e-7 of the norm lies below Zana > TolSp = 1e-7, so spodvr stops."""
    e = evaluate(load_parameters(DZ_INPUT, z0=16.78), Settings())
    check_equal(e.stopped_at, "make_initial_wave_packet", "stage")
    check_has_error(e.findings, "TolSp", "too much norm", "TolSp error")


def test_fixed_specular_analysis_line_limits_the_region(work: Path) -> None:
    """With DynamZSpAna false and ZAnaSp 21.0 the region ends one point before Z = 21.0 (index 175)."""
    e = evaluate(load_parameters(DZ_INPUT, dynamical_zsp_analysis=False, z_analysis_specular=21.0), Settings())
    check_equal(e.packet.info.i_high, 174, "last point of the region")


def test_no_packet_in_a_potential_only_run(work: Path) -> None:
    """Run mode 6 builds no wave packet; the set-up passes, and the summary has the widths but no norms."""
    e = evaluate(load_parameters(DZ_INPUT, only_potential=True), Settings())
    check(e.accepted and e.packet is None, f"accepted {e.accepted}, packet {e.packet is not None}")
    s = e.packet_summary
    check(s.sigma_k > 0.0 and s.norm_outside is None and s.tolsp_measure is None, f"summary {s}")


TESTS = [
    test_dz_packet_parameters_match_the_program,
    test_dz_packet_matches_the_program,
    test_tolsp_error_when_the_packet_reaches_the_analysis_line,
    test_fixed_specular_analysis_line_limits_the_region,
    test_no_packet_in_a_potential_only_run,
]
