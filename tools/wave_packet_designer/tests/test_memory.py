"""Tests (ctest WavePacketDesigner.Memory): the angular sizes and the CPU and GPU memory estimates against the program's
(program/angular_grid.py, memory_estimate.py; report.format_bytes)."""

from pathlib import Path

from evaluation import SYMMETRY_OFF, Resources, evaluate
from program.point_group import PHI_POINTS_MULTIPLE
from report import format_bytes
from settings import Settings
from tests.helpers import DZ_INPUT, HARDWALL_INPUT, check, check_equal, load_parameters

COMPLEX_BYTES, REAL_BYTES = 16, 8


def resources_of(path: Path, settings: Settings, **changes: object) -> Resources:
    """The resources of an input with `changes` to InputParameters fields (the set-up must reach the grids)."""
    e = evaluate(load_parameters(path, **changes), settings)
    check(e.resources is not None, f"no estimate: {e.stopped_at} {e.findings}")
    return e.resources


def test_hardwall_by_hand() -> None:
    """TEST(MemoryEstimate, ReferenceCaseByHand), 4 threads: potential 3 * 10 * 16 * 113 * 56 (theta, phi, x y,
    nzpot, r), freed after the set-up; wave function 15 * 16 * 176 * 56; specular 2 copies of 56 * 280; work
    arrays the (r, Z) unit 4 * 56 * 176 per thread."""
    memory = resources_of(HARDWALL_INPUT, Settings(threads=4)).memory
    check_equal(len(memory.items), 5, "items")
    potential = 3 * 10 * 16 * 113 * 56
    check_equal((memory.items[0].values, memory.items[0].bytes), (potential, potential * REAL_BYTES), "potential")
    check_equal((memory.items[0].during_setup, memory.items[0].during_propagation), (True, False), "potential flags")
    check_equal(memory.items[1].bytes, potential * COMPLEX_BYTES, "exp table")
    scattering = 15 * 16 * 176 * 56
    check_equal(memory.items[2].values, scattering, "scattering wave function")
    check_equal(memory.items[3].bytes, 2 * 56 * 280 * COMPLEX_BYTES, "specular wave function")
    work_values = 4 * 56 * 176
    check_equal((memory.items[4].values, memory.items[4].bytes), (work_values, 4 * work_values * COMPLEX_BYTES), "work")
    wave_functions = (scattering + 2 * 56 * 280) * COMPLEX_BYTES
    check_equal(memory.setup_bytes, potential * (REAL_BYTES + COMPLEX_BYTES) + wave_functions, "set-up total")
    expected = potential * COMPLEX_BYTES + wave_functions + 4 * work_values * COMPLEX_BYTES
    check_equal(memory.propagation_bytes, expected, "propagation total")


def test_dz_matches_the_program() -> None:
    """The blocks GRIDS (theta and phi points), ROTATIONAL BASIS, MEMORY ESTIMATE and GPU MEMORY ESTIMATE of
    spodvr for the DZ input (provenance: helpers.DZ_INPUT, except the work arrays: the committed log prints 194400
    values, the current program 25600 (expected_cases.json, case dz); 8 threads, 19 of 144 (x, y) points with
    C6v's 12 operations): every row's values and printed size, and the totals."""
    r = resources_of(DZ_INPUT, Settings(threads=8, representatives=19))
    check_equal((r.angular.n_theta, r.angular.n_phi, r.rotational_basis_size), (15, 60, 435), "angular sizes")
    cpu = [
        ("Potential V", 85500000, "652.3 MiB"),
        ("Surface symmetry maps", 1008, "3.9 KiB"),
        ("exp(-i V dt) table", 85500000, "1.3 GiB"),
        ("Wave function, scattering grid", 400896000, "6.0 GiB"),
        ("Wave function, specular grid", 10240, "320.0 KiB"),
        ("Work arrays", 25600, "3.1 MiB"),
    ]
    check_equal([(i.name, i.values, format_bytes(i.bytes)) for i in r.memory.items], cpu, "CPU rows")
    totals = (format_bytes(r.memory.setup_bytes), format_bytes(r.memory.propagation_bytes))
    check_equal(totals, ("7.9 GiB", "7.3 GiB"), "CPU totals")
    gpu = [
        ("Wave function, scattering grid", 400896000, "6.0 GiB"),
        ("Wave function, specular grid", 10240, "160.0 KiB"),
        ("exp(-i V dt) table", 85500000, "1.3 GiB"),
        ("DVR scratch (one slab)", 64800000, "988.8 MiB"),
        ("cuFFT work area", 33554432, "512.0 MiB"),
        ("BK plane", 2505600, "38.2 MiB"),
        ("Copy buffer", 4134240, "63.1 MiB"),
    ]
    check_equal([(i.name, i.values, format_bytes(i.bytes)) for i in r.gpu_memory.items], gpu, "GPU rows")
    check_equal(format_bytes(r.gpu_memory.peak_bytes), "8.8 GiB", "GPU total")


def test_large_potential_is_applied_on_the_fly() -> None:
    """TEST(MemoryEstimate, LargePotentialComputedOnTheFly): hardwall_small with nX = nY = 64 needs more than
    4 GiB for exp(-i V dt), so the table is not stored and the potential is kept."""
    memory = resources_of(HARDWALL_INPUT, Settings(threads=1), nx=64, ny=64).memory
    check(memory.items[0].values * COMPLEX_BYTES > 4 << 30, "the table would exceed 4 GiB")
    check(memory.items[0].during_propagation, "the potential is kept")
    check_equal(memory.items[1].bytes, 0, "no exp table")


def test_on_the_fly_with_symmetry_adds_the_exp_scratch() -> None:
    """The same with UseSymmetry and 2000 of 4096 representatives (still above 4 GiB): each thread also holds
    exp(i A) of the representatives of one (Z, r) point. By hand: one (Z, r) point per T_xy group (15 * 4096
    values, 960 KiB, a multiple of 64 bytes) and no momentum shift, so T_xy runs in place without a group array;
    DVR blocks of one angular-transform pass 11 * 3 * 64 = 2112 (64 columns of 528 bytes fit 58 KiB), scratch
    3 * 10 * 2000 = 60000, in all 62112 (more than the (r, Z) unit 39424)."""
    settings = Settings(threads=1, representatives=2000)
    memory = resources_of(HARDWALL_INPUT, settings, nx=64, ny=64, use_symmetry=True).memory
    check_equal(memory.items[2].bytes, 0, "no exp table")
    check_equal(memory.items[-1].values, 62112, "work arrays")


def test_potential_only_run() -> None:
    """TEST(MemoryEstimate, PotentialOnlyRun): run mode 6 holds only the potential; no GPU estimate."""
    r = resources_of(HARDWALL_INPUT, Settings(threads=4), only_potential=True)
    check_equal(len(r.memory.items), 1, "items")
    check_equal((r.memory.peak_bytes, r.memory.propagation_bytes), (r.memory.items[0].bytes, 0), "totals")
    check(r.gpu_memory is None, "no GPU estimate")


def test_mj_max_zero_means_j_max() -> None:
    """MJMax = 0 gives the basis of MJMax = JMax (process_input), so the DZ grid keeps its 60 phi points."""
    r = resources_of(DZ_INPUT, Settings(), mj_max=0)
    check_equal((r.angular.mj_max, r.angular.n_phi, r.rotational_basis_size), (28, 60, 435), "MJMax 0")


def test_phi_points_multiple_of_each_point_group() -> None:
    """TEST(PointGroup, PhiPointsMultiple): 2, 4, 6, 6, 3, 1 for C2v, C4v, C6v, C3v, C3, none."""
    expected = {"C2v": 2, "C4v": 4, "C6v": 6, "C3v": 6, "C3": 3, "none": 1}
    check_equal(PHI_POINTS_MULTIPLE, expected, "multiples")


def test_phi_points_multiple_needs_use_symmetry_and_nphisoptimal() -> None:
    """TEST(Input, SurfaceSymmetryKeys): hardwall_small (C2v) gets the multiple 2 only with UseSymmetry and
    nPhisOptimal together."""
    for changes, multiple in (
        ({}, 1),
        ({"use_symmetry": True}, 1),
        ({"use_symmetry": True, "phis_optimal": True}, 2),
        ({"phis_optimal": True}, 1),
    ):
        e = evaluate(load_parameters(HARDWALL_INPUT, **changes), Settings())
        check_equal(e.derived.phi_points_multiple, multiple, f"multiple with {changes}")


def test_basis_of_the_parity_branches() -> None:
    """The J ranges of make_angular_grid. Heteronuclear, JMax 3, MJMax 1 (TEST(RotationalBasis,
    MjMaxBelowJMax)): 4 theta points, mJ = -1: J = 1..3, mJ = 0: J = 0..3, mJ = 1: J = 1..3, so 10 functions.
    Homonuclear with odd JIn 1, JMax 5, MJMax 5: J = 1, 3, 5, so 3 theta points and 3 + 2 * (3 + 2 + 2 + 1 + 1)
    = 21 functions (mJ = 0; |mJ| = 1; 2 and 3; 4 and 5)."""
    hetero = resources_of(
        HARDWALL_INPUT, Settings(), homonuclear=False, mass2=2.0141, j_max=3, mj_max=1, j_in=1, mj_in=0
    )
    check_equal((hetero.angular.n_theta, hetero.rotational_basis_size), (4, 10), "heteronuclear")
    odd = resources_of(HARDWALL_INPUT, Settings(), j_max=5, mj_max=5, j_in=1, mj_in=0)
    check_equal((odd.angular.j_min, odd.angular.n_theta, odd.rotational_basis_size), (1, 3, 21), "odd JIn")


def test_symmetry_counts_need_use_symmetry() -> None:
    """Counts given for an input without UseSymmetry change nothing: the program stores every (x, y) point."""
    r = resources_of(HARDWALL_INPUT, Settings(representatives=3))
    check_equal((r.symmetry_source, r.symmetry.n_representatives), (SYMMETRY_OFF, 16), "symmetry")


def test_format_bytes_units() -> None:
    """TEST(FormatBytes, Units)."""
    cases = [(0, "0 B"), (1023, "1023 B"), (1024, "1.0 KiB"), (3 << 19, "1.5 MiB"), (4 << 30, "4.0 GiB")]
    cases.append((2048 << 30, "2.0 TiB"))
    for n_bytes, text in cases:
        check_equal(format_bytes(n_bytes), text, f"{n_bytes} bytes")
