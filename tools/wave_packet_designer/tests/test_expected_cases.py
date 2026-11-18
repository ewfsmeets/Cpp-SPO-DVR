"""Case expected-cases, the drift guard: the designer's grid indices, angular sizes and CPU and GPU memory
estimates against expected_cases.json. The program's tests check the same file against the C++:
TEST(MemoryEstimate, MatchesTheWavePacketDesignerCases) in tests/test_memory_estimate.cpp and
TEST(CudaBackend, MemoryEstimateMatchesTheWavePacketDesignerCases) in tests/cuda/test_cuda_backend.cpp. So when
the C++ changes one of these numbers, its test fails: update the copy in program/, then regenerate the file
with `test_designer.py expected-cases <scratch directory> --write` (README, "Keeping up with the C++").

Each case of the file is a definition written by hand (a base input, key changes, the threads, and the
symmetry counts the program computes for it) and the results that --write fills in.
"""

import json
from pathlib import Path

import program.spatial_grids
from evaluation import Evaluation, evaluate
from program.json_input import parse_input
from program.memory_estimate import MemoryEstimate
from program.spatial_grids import IndexRange
from settings import Settings
from tests.helpers import REPOSITORY, TOOL_DIR, check, check_equal, load_values

EXPECTED_CASES = TOOL_DIR / "expected_cases.json"
DEFINITION_KEYS = ("name", "about", "input", "changes", "threads", "symmetry")  # by hand; --write keeps them
PLANTED_CASE = "dz_point_below_z_asymptotic"  # the case the negative test must fail


def load_cases() -> list[dict]:
    """The cases of expected_cases.json."""
    return json.loads(EXPECTED_CASES.read_text())["cases"]


def case_input(case: dict) -> dict:
    """The case's input: its base input file (relative to the repository root) with its key changes; a change
    to null removes the key."""
    values = load_values(REPOSITORY / case["input"])
    for key, value in case["changes"].items():
        if value is None:
            values.pop(key, None)
        else:
            values[key] = value
    return values


def evaluate_case(case: dict) -> Evaluation:
    """evaluate() on the case's input, with its threads and, for UseSymmetry, its symmetry counts."""
    parsed = parse_input(case_input(case))
    check(parsed.parameters is not None, f"{case['name']}: spodvr cannot read it: {parsed.findings}")
    settings = Settings(threads=case["threads"])
    if parsed.parameters.use_symmetry:
        symmetry = case["symmetry"]
        settings = Settings(
            threads=case["threads"], representatives=symmetry["representatives"], operations=symmetry["operations"]
        )
    e = evaluate(parsed.parameters, settings)
    check(e.grids is not None and e.resources is not None, f"{case['name']}: stopped at {e.stopped_at}")
    return e


def designer_results(case: dict) -> dict:
    """What the designer computes for the case, in the file's form (the C++ tests read the same names)."""
    e = evaluate_case(case)
    g, r = e.grids, e.resources
    return {
        "grids": {
            "iz_analysis": g.iz_analysis,
            "nz_potential": g.nz_potential,
            "z_optical": index_range(g.z_optical),
            "z_specular_optical": index_range(g.z_specular_optical),
            "r_optical": index_range(g.r_optical),
            "z_optical_low": index_range(g.z_optical_low),
            "z_specular_optical_low": index_range(g.z_specular_optical_low),
        },
        "angular": {
            "n_theta": r.angular.n_theta,
            "n_phi": r.angular.n_phi,
            "rotational_basis_size": r.rotational_basis_size,
        },
        "cpu": estimate(r.memory),
        "gpu": estimate(r.gpu_memory),
    }


def index_range(points: IndexRange | None) -> list[int] | None:
    """[first, last] of an index range, None for an absorber that is off."""
    return None if points is None else [points.first, points.last]


def estimate(memory: MemoryEstimate | None) -> dict | None:
    """A memory estimate with the field names of the C++ MemoryEstimate (None: not printed for this run mode)."""
    if memory is None:
        return None
    items = [
        {
            "name": item.name,
            "values": item.values,
            "bytes": item.bytes,
            "during_setup": item.during_setup,
            "during_propagation": item.during_propagation,
        }
        for item in memory.items
    ]
    return {
        "items": items,
        "setup_bytes": memory.setup_bytes,
        "propagation_bytes": memory.propagation_bytes,
        "peak_bytes": memory.peak_bytes,
    }


def differences(expected: object, actual: object, where: str) -> list[str]:
    """Every place where the file's value and the designer's differ (exact, type included: True is not 1)."""
    if isinstance(expected, dict) and isinstance(actual, dict) and expected.keys() == actual.keys():
        return [found for key in expected for found in differences(expected[key], actual[key], f"{where}.{key}")]
    if isinstance(expected, list) and isinstance(actual, list) and len(expected) == len(actual):
        pairs = enumerate(zip(expected, actual, strict=True))
        return [found for i, (x, y) in pairs for found in differences(x, y, f"{where}[{i}]")]
    if type(expected) is type(actual) and expected == actual:
        return []
    return [f"{where}: file {expected!r}, designer {actual!r}"]


def all_differences(cases: list[dict]) -> list[str]:
    """The differences of every case, each starting with the case's name."""
    found = []
    for case in cases:
        results = designer_results(case)
        found += differences({key: case.get(key) for key in results}, results, case["name"])
    return found


def test_designer_matches_the_expected_cases(work: Path) -> None:
    """Every case of expected_cases.json: the designer's grid indices, angular sizes and memory estimates equal
    the file's, which the program's tests check against the C++."""
    found = all_differences(load_cases())
    check(
        not found, f"{len(found)} difference(s); regenerate with --write once the copy is right:\n" + "\n".join(found)
    )


def test_the_cases_cover_the_estimate_paths(work: Path) -> None:
    """The cases keep covering what the estimates distinguish: UseSymmetry with a reduction, a potential-only
    run (no GPU estimate), a potential above the 4 GiB table limit (kept, no exp table), the low-Z absorbers,
    and both molecule kinds (half or all theta points)."""
    cases = load_cases()
    inputs = [case_input(case) for case in cases]
    check(
        any(c["symmetry"]["representatives"] < v["nX"] * v["nY"] for c, v in zip(cases, inputs, strict=True)),
        "symmetry",
    )
    check(any(c["gpu"] is None and c["cpu"] is not None for c in cases), "a potential-only run")
    check(any(c["gpu"] is not None and c["cpu"]["items"][1]["bytes"] == 0 for c in cases), "exp table not stored")
    check(any(c["grids"]["z_optical_low"] is not None for c in cases), "low-Z absorbers")
    check(any(v["Mass1"] != v["Mass2"] for v in inputs) and any(v["Mass1"] == v["Mass2"] for v in inputs), "masses")


def test_a_planted_tolerance_slip_fails_one_case(work: Path) -> None:
    """Negative test: the program places the analysis line and nzpot with the tolerant comparison rge (1e-12).
    A copy with plain >= instead (a realistic porting slip) must fail exactly the case whose Z grid has a point
    1e-13 below Z_asymptotic and Zana (on the other DZ grids 125 * 0.12 is exactly 15, so they cannot see it), on
    16 fields: the two indices, the potential and the exp table (values, bytes), the CPU totals (3), the GPU exp
    table and DVR slab (values, bytes) and the GPU totals (3)."""
    original = program.spatial_grids.rge
    program.spatial_grids.rge = lambda a, b: a >= b
    try:
        found = all_differences(load_cases())
    finally:
        program.spatial_grids.rge = original
    check_equal(sorted({line.split(".")[0] for line in found}), [PLANTED_CASE], "cases that fail")
    check_equal(len(found), 16, "fields that fail:\n" + "\n".join(found))


def write() -> None:
    """Regenerates expected_cases.json: keeps each case's definition and recomputes the designer's results."""
    data = json.loads(EXPECTED_CASES.read_text())
    data["cases"] = [
        {**{key: case[key] for key in DEFINITION_KEYS}, **designer_results(case)} for case in data["cases"]
    ]
    EXPECTED_CASES.write_text(to_json(data) + "\n")


def to_json(value: object, indent: str = "") -> str:
    """JSON that reads as a table: an object or list of plain values (a grid range, a memory item) on one line,
    everything else one entry per line."""
    if is_flat(value):
        return json.dumps(value)
    inner = indent + "  "
    if isinstance(value, dict):
        rows = [f"{inner}{json.dumps(key)}: {to_json(entry, inner)}" for key, entry in value.items()]
        return "{\n" + ",\n".join(rows) + "\n" + indent + "}"
    rows = [inner + to_json(entry, inner) for entry in value]
    return "[\n" + ",\n".join(rows) + "\n" + indent + "]"


def is_flat(value: object) -> bool:
    """A plain value, a list of plain values, or an object whose entries are plain values or such lists."""
    if isinstance(value, list):
        return all(is_plain(entry) for entry in value)
    if isinstance(value, dict):
        return all(is_plain(entry) or is_flat(entry) and isinstance(entry, list) for entry in value.values())
    return True


def is_plain(value: object) -> bool:
    """A number, string, boolean or null."""
    return not isinstance(value, (dict, list))


TESTS = [
    test_designer_matches_the_expected_cases,
    test_the_cases_cover_the_estimate_paths,
    test_a_planted_tolerance_slip_fails_one_case,
]
