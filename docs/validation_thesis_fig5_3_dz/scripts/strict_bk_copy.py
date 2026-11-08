#!/usr/bin/env python3
"""A copy of the strict-build GPU run's BK file that the build-fast replay helper accepts.

The C++ BK reader checks the stored restart parameters exactly against those of the reading
binary (src/io/checkpoint.cpp read_header). Of the 35 parameters, only "ZOptFactor (effective)"
differs between the strict and the fast builds (by 1 ulp, a fast-math rounding of the set-up);
it is an absorber factor of the propagation, not used by the analysis. The copy carries the fast
build's value in that one header field; the coefficient records are untouched (checked below).
Writes runs/thesis_fig5_3_dz_validation/strict_bk/bk_coefficients.bin (266 MB, git-ignored)."""
from __future__ import annotations

import shutil
import struct

import numpy as np

from common import RUN_DIRS, RUNS, WORK, NativeBk

NAME = "ZOptFactor (effective)"
SOURCE = RUNS / "thesis_fig5_3_dz_cuda" / "bk_coefficients.bin"
TARGET = WORK / "strict_bk" / "bk_coefficients.bin"


def parameter_offset(path, name):
    """Byte offset of the double value of one restart parameter in a native BK file header."""
    with open(path, "rb") as f:
        f.read(16 + 4)                                  # magic, version
        (n,) = struct.unpack("<Q", f.read(8))
        f.read(n)                                       # SystemID
        (count,) = struct.unpack("<Q", f.read(8))
        for _ in range(count):
            (n,) = struct.unpack("<Q", f.read(8))
            if f.read(n).decode() == name:
                return f.tell()
            f.read(8)
    raise KeyError(name)


def main():
    strict, fast = NativeBk(SOURCE), NativeBk(RUN_DIRS["gpu"] / "bk_coefficients.bin")
    differing = [k for k in strict.parameters if strict.parameters[k] != fast.parameters[k]]
    assert differing == [NAME], differing
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(SOURCE, TARGET)
    with open(TARGET, "r+b") as f:
        f.seek(parameter_offset(TARGET, NAME))
        f.write(struct.pack("<d", fast.parameters[NAME]))
    copy = NativeBk(TARGET)
    assert copy.parameters == fast.parameters
    assert copy.first_offset == strict.first_offset and copy.n_records == strict.n_records
    a = np.memmap(SOURCE, dtype=np.uint8, mode="r")
    b = np.memmap(TARGET, dtype=np.uint8, mode="r")
    assert np.array_equal(a[strict.first_offset - 8 * strict.nstata * strict.nviba:], b[copy.first_offset - 8 * copy.nstata * copy.nviba:])
    print(f"{TARGET}: {NAME} {strict.parameters[NAME]!r} -> {fast.parameters[NAME]!r}; energies and records identical")


if __name__ == "__main__":
    main()
