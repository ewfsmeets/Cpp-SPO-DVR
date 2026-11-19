#!/usr/bin/env python3
"""Tests of the wave-packet designer, one ctest per case (WavePacketDesigner.<Case>, tests/CMakeLists.txt).

Each case is a module in tests/ with an explicit list TESTS of test functions `test_x(work: Path)` (work: a
scratch directory of the case); a test raises TestFailure (tests/helpers.py) when it fails. To add a test:
write the function in the case's module and add it to TESTS.

Cases: see CASES below.
Usage: python3 -I -B tools/wave_packet_designer/test_designer.py <case> <work directory> [--write]
--write regenerates the expected data of a case that has a `write` function (expected-cases: the drift
guard's expected_cases.json) instead of testing it.
Exit status: 0 every test of the case passes, 1 a test fails, 2 usage error.
Environment: none read. Run time: about 20 s (plots), 10 s (design), under 5 s for the other cases.
Needs: numpy and matplotlib.
Docs: tools/wave_packet_designer/README.md ("Tests")
"""

import argparse
import importlib
import os
import shutil
import sys
import traceback
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))

from tests.helpers import TestFailure  # noqa: E402

# Case name -> module in tests/.
CASES = {
    "fft-sizes": "tests.test_fft_sizes",
    "real-compare": "tests.test_real_compare",
    "input-reading": "tests.test_input_reading",
    "validate": "tests.test_validate",
    "absorbers": "tests.test_absorbers",
    "grids": "tests.test_grids",
    "packet": "tests.test_packet",
    "memory": "tests.test_memory",
    "expected-cases": "tests.test_expected_cases",
    "gas-phase": "tests.test_gas_phase",
    "design": "tests.test_design",
    "plots": "tests.test_plots",
    "input-file": "tests.test_input_file",
    "cli": "tests.test_cli",
    "structure": "tests.test_structure",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", choices=sorted(CASES))
    parser.add_argument("work", type=Path, help="scratch directory (deleted and recreated)")
    parser.add_argument("--write", action="store_true", help="regenerate the case's expected data")
    args = parser.parse_args()

    shutil.rmtree(args.work, ignore_errors=True)
    args.work.mkdir(parents=True)
    # Figures without a window, and matplotlib's cache in the scratch directory (set before matplotlib loads).
    os.environ["MPLBACKEND"] = "Agg"
    os.environ["MPLCONFIGDIR"] = str(args.work / "matplotlib")
    module = importlib.import_module(CASES[args.case])
    if args.write:
        return write_expected_data(module, args.case)
    failures = 0
    for test in module.TESTS:
        try:
            test(args.work)
        except TestFailure as failure:
            failures += 1
            print(f"FAIL: {args.case}/{test.__name__}: {failure}")
        except Exception:  # an unexpected exception is a failure too, with its traceback
            failures += 1
            print(f"FAIL: {args.case}/{test.__name__}: unexpected exception\n{traceback.format_exc()}")
    if failures:
        print(f"FAIL: {args.case}: {failures} of {len(module.TESTS)} tests failed")
        return 1
    print(f"PASS: {args.case} ({len(module.TESTS)} tests)")
    return 0


def write_expected_data(module: object, case: str) -> int:
    """Runs the case's `write` function, which regenerates its expected data from the designer."""
    if not hasattr(module, "write"):
        print(f"error: case {case} has no expected data to write", file=sys.stderr)
        return 2
    try:
        module.write()
    except TestFailure as failure:
        print(f"FAIL: {case}: {failure}")
        return 1
    print(f"WROTE: {case}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
