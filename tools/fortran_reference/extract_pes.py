#!/usr/bin/env python3
"""Extract a compact text reference from a Fortran SPO-DVR binary .pes file.

The .pes file is written by Potential.f90 WritePesToFile as gfortran unformatted sequential
records (4-byte record markers, default INTEGER/LOGICAL = 4 bytes, REAL = 8 bytes because of
-fdefault-real-8). The stored potential is V in Hartree *before* premultiplication by -dt,
layout vdvr(itp, ixy, izr) with itp = iphi + nphi*(itheta-1), izr = ir + nr*(iz-1).

With surface symmetry (UseSymmetry) only the (x, y) points of a fundamental region are stored
and symmetryIndex maps every full-grid point to its stored equivalent; the script expands it.

An SCM run (SystemID CRPDistrtSurface) writes six more records between symmetryIndex and the
potential (DistortedSurfaceWrapper.f90 WriteDistPESSetup, SCM_Wrapper.f90 SetPESVariables_SCM).
They are detected by counting records, checked against their formats and copied into '#' lines.
The C++ has no .pes cache (SCM D11, docs/fortran_deviations.md): only this script reads them.

Output (text, full precision):
  header line with nzpot nr nx ny ntheta nphi
  for a potential that depends only on (Z, r) (the hard wall, no symmetry): the theta and phi
  grid points and V(iz, ir) for every (Z, r) point;
  otherwise (format 6d): use_symmetry and the number of computed (x, y) points, the theta and
  phi grid points, and V(itheta, iphi, ix, iy, iz, ir) for a deterministic subsample (every k-th
  flat index, k coprime to the (theta, phi, x, y) block size, about 40000 values) plus the
  first point of every (Z, r) block.

Usage: extract_pes.py <file.pes> <output.txt>
Exit status: 0 written, 1 a file of another layout or one that cannot be read or written
(error: <what>; <fix>), 2 bad arguments.
Environment: none read. Run time: seconds (the reference cases). Called by run_setup_reference.sh;
docs: tools/fortran_reference/README.md.
"""
import math
import struct
import sys
from pathlib import Path

# The gfortran record reader is refdata.py's, in this directory (python3 -I leaves it off sys.path).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from refdata import FormatError, read_records  # noqa: E402

# The records before the potential: dimensions, angular sizes, Z/r and x/y grid ends, thetas,
# phis, the switching function, the symmetry flag and size, symmetryIndex.
HEADER_RECORDS = 9

# The six records of an SCM run, in file order: the names SetPESVariables_SCM writes and their
# formats (logicals as 4-byte integers 0/1).
SCM_RECORDS = (
    ("ParameterVersion", "<i"),
    ("EnableExpansion EnabReadPosFile", "<2i"),
    ("IdLatticeConstant DistLatticeConstant", "<2d"),
    ("CutOff UsedLayers NumSlabConfigs", "<d2i"),
    ("SlabAtomSize(1:3)", "<3i"),
    ("SCMSlabRefocusMethod", "<i"),
)


def unpack(fmt, payload):
    if struct.calcsize(fmt) != len(payload):
        raise FormatError(f"record size {len(payload)} does not match format {fmt}")
    return struct.unpack(fmt, payload)


def format_value(value):
    """A value of an SCM record as text: a REAL as the potential values (.17e: 18 significant
    digits), an INTEGER as is."""
    return f"{value:.17e}" if isinstance(value, float) else str(value)


def scm_lines(records):
    """The six SCM records as '#' lines (empty for a file without them)."""
    if not records:
        return ""
    lines = ["# SCM records (SCM_Wrapper.f90 SetPESVariables_SCM; logicals as 0/1):\n"]
    for (names, fmt), payload in zip(SCM_RECORDS, records):
        values = " ".join(format_value(v) for v in unpack(fmt, payload))
        lines.append(f"#   {names}: {values}\n")
    return "".join(lines)


def expand_symmetric(block, symmetry_index):
    """Full-grid block from the reduced one: V(i) = reduced(symmetryIndex(i)) (1-based)."""
    return [block[k - 1] for k in symmetry_index]


def sample_indices(n_values, block_size):
    """Deterministic subsample of the flat (itp, ixy, izr) indices: every k-th one, with k
    coprime to the block size so that all (theta, phi, x, y) combinations are visited."""
    stride = max(1, n_values // 40000)
    while math.gcd(stride, block_size) != 1:
        stride += 1
    return set(range(0, n_values, stride))


def write_zr_only(out, header, thetas, phis, values, nr, scm):
    """Hard wall format: V(iz, ir), the potential being constant over (x, y, theta, phi)."""
    out.write("# Extracted by tools/fortran_reference/extract_pes.py from the Fortran .pes file.\n")
    out.write("# V in Hartree (before the -dt premultiplication); constant over (x, y, theta, phi).\n")
    out.write(scm)
    out.write("# nzpot nr nx ny ntheta nphi\n")
    out.write(header)
    out.write("# theta points\n")
    out.write("".join(f"{t:.17e}\n" for t in thetas))
    out.write("# phi points\n")
    out.write("".join(f"{p:.17e}\n" for p in phis))
    out.write("# iz ir V(iz,ir)  (1-based indices, ir fastest)\n")
    for izr, block in enumerate(values):
        iz, ir = divmod(izr, nr)
        out.write(f"{iz + 1} {ir + 1} {block[0]:.17e}\n")


def write_6d(out, header, thetas, phis, values, dims, use_symmetry, n_xy_computed, scm):
    """6D format: a deterministic subsample of V(itheta, iphi, ix, iy, iz, ir) plus the whole
    (Z, r) table at the first (x, y, theta, phi) point."""
    nzpot, nr, nx, ny, ntheta, nphi = dims
    block_size = ntheta * nphi * nx * ny
    samples = sample_indices(block_size * nzpot * nr, block_size)
    out.write("# Extracted by tools/fortran_reference/extract_pes.py from the Fortran .pes file.\n")
    out.write("# V in Hartree (before the -dt premultiplication); 6D subsample (format 6d).\n")
    out.write(scm)
    out.write("# nzpot nr nx ny ntheta nphi\n")
    out.write(header)
    out.write("# use_symmetry (0/1) and the number of (x, y) points for which V was computed\n")
    out.write(f"{int(bool(use_symmetry))} {n_xy_computed}\n")
    out.write("# theta points\n")
    out.write("".join(f"{t:.17e}\n" for t in thetas))
    out.write("# phi points\n")
    out.write("".join(f"{p:.17e}\n" for p in phis))
    out.write("# itheta iphi ix iy iz ir V  (1-based; full grid, expanded through symmetryIndex)\n")
    for izr, block in enumerate(values):
        iz, ir = divmod(izr, nr)
        for i, v in enumerate(block):
            if i == 0 or (i + block_size * izr) in samples:
                ixy, itp = divmod(i, ntheta * nphi)
                itheta, iphi = divmod(itp, nphi)
                ix, iy = divmod(ixy, ny)
                out.write(f"{itheta + 1} {iphi + 1} {ix + 1} {iy + 1} {iz + 1} {ir + 1} {v:.17e}\n")


def extract(pes_path, output_path):
    records = read_records(Path(pes_path).read_bytes())
    if not records:
        raise FormatError("the file holds no records")
    nzpot, nr, nx, ny = unpack("<4i", records[0])
    # Records beyond the header and the nzpot * nr potential blocks: none, or the six of SCM.
    n_extra = len(records) - HEADER_RECORDS - nzpot * nr
    if n_extra not in (0, len(SCM_RECORDS)):
        raise FormatError(f"{len(records)} records, not {HEADER_RECORDS} + {nzpot * nr} (nzpot * nr) "
                          f"or {len(SCM_RECORDS)} more (SCM)")
    ntheta, nphi = unpack("<2i", records[1])
    # records[2]: zgrid(1), rgrid(1), zgrid(nzpot), rgrid(nr); records[3]: the x, y grid ends
    thetas = unpack(f"<{ntheta}d", records[4])
    phis = unpack(f"<{nphi}d", records[5])
    # records[6]: lswitch, zmin_switch, zmax_switch
    use_symmetry, n_xy_computed = unpack("<2i", records[7])
    full_size = ntheta * nphi * nx * ny
    symmetry_index = unpack(f"<{full_size}i", records[8])
    scm = scm_lines(records[HEADER_RECORDS : HEADER_RECORDS + n_extra])

    stored_size = ntheta * nphi * n_xy_computed
    values = []
    for payload in records[HEADER_RECORDS + n_extra :]:
        block = unpack(f"<{stored_size}d", payload)
        values.append(expand_symmetric(block, symmetry_index) if use_symmetry else list(block))

    header = f"{nzpot} {nr} {nx} {ny} {ntheta} {nphi}\n"
    zr_only = not use_symmetry and all(v == block[0] for block in values for v in block)
    with open(output_path, "w") as out:
        if zr_only:
            write_zr_only(out, header, thetas, phis, values, nr, scm)
        else:
            dims = (nzpot, nr, nx, ny, ntheta, nphi)
            write_6d(out, header, thetas, phis, values, dims, use_symmetry, n_xy_computed, scm)


def main(argv):
    if len(argv) != 2:
        print("usage: extract_pes.py <file.pes> <output.txt>", file=sys.stderr)
        return 2
    try:
        extract(argv[0], argv[1])
    except (FormatError, struct.error) as error:
        print(f"error: {argv[0]}: {error}; name a .pes file written by the reference install", file=sys.stderr)
        return 1
    except OSError as error:
        print(f"error: {error}; name a .pes file that exists and an output path that can be written", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
