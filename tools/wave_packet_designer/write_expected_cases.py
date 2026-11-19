"""Regenerates expected_cases.json, the drift guard between the designer and the C++ (README, "Keeping up
with the C++"): keeps each case's definition and recomputes the designer's grid indices, angular sizes
and memory estimates (tests/test_expected_cases.py, function write). The test that checks the file is
the ctest WavePacketDesigner.ExpectedCases.

Usage: uv run python -I -B tools/wave_packet_designer/write_expected_cases.py
Exit status: 0 written, 1 a check failed while computing the results.
Environment: none read. Run time: a few seconds. Needs: numpy.
Docs: tools/wave_packet_designer/README.md ("Keeping up with the C++")
"""

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's modules and its tests package

from tests.helpers import TestFailure  # noqa: E402
from tests.test_expected_cases import write  # noqa: E402


def main() -> int:
    """Writes expected_cases.json; returns the exit status."""
    try:
        write()
    except TestFailure as failure:
        print(f"FAIL: {failure}")
        return 1
    print("WROTE: expected_cases.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
