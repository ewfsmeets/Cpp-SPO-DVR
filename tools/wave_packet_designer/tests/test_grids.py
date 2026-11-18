"""Case grids: grid indices against the program's, the Z_asymptotic boundary, the grid-point trap
(program/spatial_grids.py)."""

from fractions import Fraction
from pathlib import Path

from evaluation import evaluate
from program.spatial_grids import (
    IndexRange,
    SplitGridPoint,
    UniformGrid,
    add_split_points,
    fortran_grid_point,
    round_to_quad,
)
from program.unit_conversion import hartree_to_ev
from settings import Settings
from tests.helpers import DZ_INPUT, HARDWALL_INPUT, check, check_equal, load_parameters


def test_dz_grids_match_the_program(work: Path) -> None:
    """The block GRIDS of spodvr for the DZ input (provenance: helpers.DZ_INPUT): analysis at index 125, 125
    points in the interaction region, absorbers at 125..159, 184..255 and 18..39, and the largest kinetic energy
    in r 3.654734 eV."""
    g = evaluate(load_parameters(DZ_INPUT), Settings()).grids
    check_equal((g.iz_analysis, g.z_analysis, g.nz_potential), (125, 15.0, 125), "analysis index, Zana, nzpot")
    check_equal(g.z_optical, IndexRange(125, 159), "ZOpt indices")
    check_equal(g.z_specular_optical, IndexRange(184, 255), "ZSpOpt indices")
    check_equal(g.r_optical, IndexRange(18, 39), "ROpt indices")
    check_equal(f"{hartree_to_ev(g.kr.max_kinetic_energy):.6f}", "3.654734", "kinetic energy at pi/dR")


def test_hardwall_interaction_region(work: Path) -> None:
    """113 of the 176 Z points of hardwall_small lie below Z_asymptotic = 7 (TEST(MemoryEstimate,
    ReferenceCaseByHand))."""
    check_equal(evaluate(load_parameters(HARDWALL_INPUT), Settings()).grids.nz_potential, 113, "nzpot")


def test_point_within_tolerance_of_z_asymptotic_is_gas_phase(work: Path) -> None:
    """With ZStart = -1e-13 point 125 lies 1e-13 below Z_asymptotic = 15: rge counts it as gas phase, so
    nzpot stays 125 (a plain < would give 126)."""
    grids = evaluate(load_parameters(DZ_INPUT, z_start=-1e-13), Settings()).grids
    check(grids.z_scattering.point(125) < 15.0, "point 125 lies below 15.0")
    check_equal(grids.nz_potential, 125, "nzpot")


def test_quad_rounding(work: Path) -> None:
    """round_to_quad keeps 113 significant bits and rounds ties to even."""
    check_equal(round_to_quad(Fraction(2**113 + 1)), Fraction(2**113), "tie rounds to even")
    check_equal(round_to_quad(Fraction(2**113 + 3)), Fraction(2**113 + 4), "tie rounds to even (up)")
    check_equal(round_to_quad(Fraction(3, 4)), Fraction(3, 4), "exact value")


def test_grid_point_trap_of_the_input_guide(work: Path) -> None:
    """inputs/README.md: with ZStart -2 and dZ 0.08, point 30 is 0.3999999999999999 here and 0.4 in the
    Fortran, so a low-Z absorber from 0.4 contains it only in the Fortran."""
    grid = UniformGrid(-2.0, 0.08, 64)
    check_equal(grid.point(30), 0.3999999999999999, "this program's point 30")
    check_equal(fortran_grid_point(grid, 30), 0.4, "the Fortran's point 30")
    split: list[SplitGridPoint] = []
    add_split_points(split, grid, "ZOptLowStart", 0.4, high_side=False)
    check_equal([(s.index, s.inside_here) for s in split], [(30, False)], "split points")


TESTS = [
    test_dz_grids_match_the_program,
    test_hardwall_interaction_region,
    test_point_within_tolerance_of_z_asymptotic_is_gas_phase,
    test_quad_rounding,
    test_grid_point_trap_of_the_input_guide,
]
