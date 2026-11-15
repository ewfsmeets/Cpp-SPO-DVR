#!/usr/bin/env python3
"""Tests of extract_pes.py on synthetic .pes files (ctest FortranReference.ExtractPes*).

Each case writes small .pes files in the record layout of Potential.f90 WritePesToFile under WORK
and runs extract_pes.py as a program:
  without-scm  a (Z, r)-only, a 6D and a symmetry-reduced 6D file: each output has the sha256 that
               extract_pes.py gave before it learnt the SCM records (measured 2026-11-14), so the
               files of every case without SCM are extracted as before
  scm          the 6D file with the six SCM records of SCM_Wrapper.f90 SetPESVariables_SCM: the
               same output plus their '#' lines, which the C++ reader skips
  refusals     five or seven extra records, an SCM record of the wrong size, a truncated file,
               a trailer that differs from its record's length, a subrecord marker (-1), a
               two-byte fragment, an empty file, a missing file: exit 1 with "error: <what>; <fix>"
               and no output file; bad arguments: exit 2

Usage: python3 -I -B tools/fortran_reference/test_extract_pes.py CASE WORK
  CASE  one of the cases above;  WORK  scratch directory (recreated)
Exit status: 0 the case behaved as expected, 1 it did not, 2 bad arguments.
Environment: none read. Run time: under a second. Docs: tools/fortran_reference/README.md.
"""

import hashlib
import shutil
import struct
import subprocess
import sys
from pathlib import Path

EXTRACT_PES = Path(__file__).resolve().parent / "extract_pes.py"
CASES = ("without-scm", "scm", "refusals")

# nzpot nr nx ny ntheta nphi: 216 values, all in the subsample (fewer than 40000).
DIMS = (3, 2, 2, 3, 2, 3)

# sha256 of the output of extract_pes.py at the commit before the SCM records (efe0034) for the
# files of without_scm_files().
DIGESTS_BEFORE_SCM = {
    "zr_only": "3e2773a5455e7a72ee99d7362e3039a8ef8d5947d503fbea19372e90fe134979",
    "full": "a6bc78cab00536d3177a904e303dffa674ebcfb48bd4a9f3c9276cf115e7e71e",
    "symmetric": "ebae34993d1b7533927f28edbc0e44eb85fa370c0c25448f33e09b11529120ea",
}

# The six SCM records as SetPESVariables_SCM writes them: ParameterVersion; EnableExpansion,
# EnabReadPosFile; IdLatticeConstant, DistLatticeConstant; CutOff, UsedLayers, NumSlabConfigs;
# SlabAtomSize(1:3); SCMSlabRefocusMethod.
SCM_RECORDS = (struct.pack("<i", 2), struct.pack("<2i", 1, 1), struct.pack("<2d", 6.8164, 6.8333),
               struct.pack("<d2i", 16.0, 2, 1), struct.pack("<3i", 4, 4, 6), struct.pack("<i", 0))
SCM_LINES = """\
# SCM records (SCM_Wrapper.f90 SetPESVariables_SCM; logicals as 0/1):
#   ParameterVersion: 2
#   EnableExpansion EnabReadPosFile: 1 1
#   IdLatticeConstant DistLatticeConstant: 6.81639999999999979e+00 6.83330000000000037e+00
#   CutOff UsedLayers NumSlabConfigs: 1.60000000000000000e+01 2 1
#   SlabAtomSize(1:3): 4 4 6
#   SCMSlabRefocusMethod: 0
"""


def record(payload):
    """A gfortran unformatted sequential record: the payload between two 4-byte lengths."""
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


def pes_records(kind):
    """The records of a .pes file: 'zr_only' (V constant over x, y, theta, phi), 'full' (6D) or
    'symmetric' (6D, three of the six (x, y) points stored, symmetryIndex expanding them)."""
    nzpot, nr, nx, ny, ntheta, nphi = DIMS
    n_xy = 3 if kind == "symmetric" else nx * ny
    n_angles = ntheta * nphi
    if kind == "symmetric":  # (x, y) point ixy is stored at ixy mod 3
        index = [(ixy % 3) * n_angles + itp + 1 for ixy in range(nx * ny) for itp in range(n_angles)]
    else:
        index = list(range(1, n_angles * nx * ny + 1))
    header = [struct.pack("<4i", nzpot, nr, nx, ny), struct.pack("<2i", ntheta, nphi),
              struct.pack("<4d", -1.0, 0.4, 2.0, 1.4), struct.pack("<4d", 0.0, 0.0, 2.4, 3.6),
              struct.pack(f"<{ntheta}d", 0.7, 2.4), struct.pack(f"<{nphi}d", 0.0, 2.1, 4.2),
              struct.pack("<i2d", 1, 11.3, 11.9), struct.pack("<2i", int(kind == "symmetric"), n_xy),
              struct.pack(f"<{len(index)}i", *index)]
    blocks = []
    for izr in range(nzpot * nr):
        if kind == "zr_only":
            values = [0.1 * izr - 0.25] * (n_angles * n_xy)
        else:
            values = [(k + 1) / 7.0 - 0.3 * izr for k in range(n_angles * n_xy)]
        blocks.append(struct.pack(f"<{len(values)}d", *values))
    return header, blocks


def write_pes(path, header, blocks, extra=()):
    """The .pes file with the records `extra` between symmetryIndex and the potential blocks."""
    path.write_bytes(b"".join(record(r) for r in [*header, *extra, *blocks]))


def extract(work, name):
    """Runs extract_pes.py on WORK/<name>.pes into WORK/<name>.txt: (exit status, stderr)."""
    result = subprocess.run([sys.executable, "-I", "-B", str(EXTRACT_PES), str(work / f"{name}.pes"),
                             str(work / f"{name}.txt")], capture_output=True, text=True)
    return result.returncode, result.stderr


def check(condition, message):
    if not condition:
        print(f"FAIL: {message}", file=sys.stderr)
        sys.exit(1)


def extracted(work, name):
    status, stderr = extract(work, name)
    check(status == 0, f"{name}: exit {status}: {stderr}")
    return (work / f"{name}.txt").read_text()


def without_scm_case(work):
    for kind, digest in DIGESTS_BEFORE_SCM.items():
        write_pes(work / f"{kind}.pes", *pes_records(kind))
        text = extracted(work, kind)
        actual = hashlib.sha256(text.encode()).hexdigest()
        check(actual == digest, f"{kind}: output sha256 {actual}, before the SCM records {digest}")
        check("SCM" not in text, f"{kind}: SCM lines in the output of a file without SCM records")


def scm_case(work):
    header, blocks = pes_records("full")
    write_pes(work / "full.pes", header, blocks)
    write_pes(work / "scm.pes", header, blocks, SCM_RECORDS)
    plain = extracted(work, "full").splitlines(keepends=True)
    expected = "".join(plain[:2]) + SCM_LINES + "".join(plain[2:])
    check(extracted(work, "scm") == expected, "the SCM output is not the 6D output plus the SCM lines")


def expect_refusal(work, name, pattern):
    status, stderr = extract(work, name)
    check(status == 1 and stderr.startswith("error: ") and pattern in stderr and "; name a " in stderr,
          f"{name}: exit {status}, expected 1 with \"{pattern}\": {stderr}")
    check(not (work / f"{name}.txt").exists(), f"{name}: a refused extraction wrote an output file")
    print(f"refused as expected: {stderr.strip()}")


def refusals_case(work):
    header, blocks = pes_records("full")
    write_pes(work / "five.pes", header, blocks, SCM_RECORDS[:5])
    expect_refusal(work, "five", "20 records, not 9 + 6 (nzpot * nr) or 6 more (SCM)")
    write_pes(work / "seven.pes", header, blocks, (*SCM_RECORDS, struct.pack("<i", 0)))
    expect_refusal(work, "seven", "22 records, not 9 + 6")
    # CutOff written as an INTEGER: the 4th SCM record has 12 bytes, not those of "<d2i".
    wrong = (*SCM_RECORDS[:3], struct.pack("<3i", 16, 2, 1), *SCM_RECORDS[4:])
    write_pes(work / "wrong_size.pes", header, blocks, wrong)
    expect_refusal(work, "wrong_size", "record size 12 does not match format <d2i")
    # The last block without its trailing length.
    write_pes(work / "truncated.pes", header, blocks)
    data = (work / "truncated.pes").read_bytes()
    (work / "truncated.pes").write_bytes(data[:-4])
    expect_refusal(work, "truncated", f"the record at byte {len(data) - len(record(blocks[-1]))} has no matching end marker")
    # The trailer of the second record one byte longer than the record.
    write_pes(work / "trailer.pes", header, blocks)
    data = bytearray((work / "trailer.pes").read_bytes())
    second = len(record(header[0]))
    data[second + 4 + len(header[1])] += 1
    (work / "trailer.pes").write_bytes(bytes(data))
    expect_refusal(work, "trailer", f"the record at byte {second} has no matching end marker")
    # A record marker of -1: gfortran's subrecords (records over 2 GiB).
    write_pes(work / "subrecords.pes", header, blocks)
    data = (work / "subrecords.pes").read_bytes()
    (work / "subrecords.pes").write_bytes(struct.pack("<i", -1) + data[4:])
    expect_refusal(work, "subrecords", "subrecords (records over 2 GiB) at byte 0 are not supported")
    # Two bytes, shorter than a record marker.
    (work / "fragment.pes").write_bytes(b"\x01\x00")
    expect_refusal(work, "fragment", "fragment.pes: a truncated record marker at byte 0")
    (work / "empty.pes").write_bytes(b"")
    expect_refusal(work, "empty", "empty.pes: the file holds no records")
    expect_refusal(work, "missing", "No such file or directory")
    result = subprocess.run([sys.executable, "-I", "-B", str(EXTRACT_PES), "only-one-argument"],
                            capture_output=True, text=True)
    check(result.returncode == 2 and result.stderr.startswith("usage: "),
          f"one argument: exit {result.returncode}, expected 2 with a usage line: {result.stderr}")


def main(argv):
    if len(argv) != 2 or argv[0] not in CASES:
        print(f"usage: test_extract_pes.py {{{','.join(CASES)}}} WORK", file=sys.stderr)
        return 2
    work = Path(argv[1])
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    {"without-scm": without_scm_case, "scm": scm_case, "refusals": refusals_case}[argv[0]](work)
    print(f"PASS: {argv[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
