#!/usr/bin/env python3
"""Extract a compact text reference from a Fortran SPO-DVR binary .pes file.

The .pes file is written by Potential.f90 WritePesToFile as gfortran unformatted sequential
records (4-byte record markers, default INTEGER/LOGICAL = 4 bytes, REAL = 8 bytes because of
-fdefault-real-8). The stored potential is V in Hartree *before* premultiplication by -dt,
layout vdvr(itp, ixy, izr) with itp = iphi + nphi*(itheta-1), izr = ir + nr*(iz-1).

With surface symmetry (UseSymmetry) only the (x, y) points of a fundamental region are stored
and symmetryIndex maps every full-grid point to its stored equivalent; the script expands it.

Output (text, full precision):
  header line with nzpot nr nx ny ntheta nphi
  for a potential that depends only on (Z, r) (the hard wall, no symmetry): the theta and phi
  grid points and V(iz, ir) for every (Z, r) point;
  otherwise (format 6d): use_symmetry and the number of computed (x, y) points, the theta and
  phi grid points, and V(itheta, iphi, ix, iy, iz, ir) for a deterministic subsample (every k-th
  flat index, k coprime to the (theta, phi, x, y) block size, about 40000 values) plus the
  first point of every (Z, r) block.

Usage: extract_pes.py <file.pes> <output.txt>
"""
import math
import struct
import sys


def read_records(path):
    """Yield the payload of every record of a gfortran unformatted sequential file."""
    with open(path, "rb") as f:
        data = f.read()
    pos = 0
    while pos < len(data):
        (length,) = struct.unpack_from("<i", data, pos)
        if length < 0:
            raise ValueError("subrecords (records > 2 GiB) are not supported")
        payload = data[pos + 4 : pos + 4 + length]
        (trailer,) = struct.unpack_from("<i", data, pos + 4 + length)
        if trailer != length:
            raise ValueError(f"corrupt record at byte {pos}")
        yield payload
        pos += length + 8


def unpack(fmt, payload):
    if struct.calcsize(fmt) != len(payload):
        raise ValueError(f"record size {len(payload)} does not match format {fmt}")
    return struct.unpack(fmt, payload)


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


def write_zr_only(out, header, thetas, phis, values, nr):
    """Hard wall format: V(iz, ir), the potential being constant over (x, y, theta, phi)."""
    out.write("# Extracted by tools/fortran_reference/extract_pes.py from the Fortran .pes file.\n")
    out.write("# V in Hartree (before the -dt premultiplication); constant over (x, y, theta, phi).\n")
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


def write_6d(out, header, thetas, phis, values, dims, use_symmetry, n_xy_computed):
    """6D format: a deterministic subsample of V(itheta, iphi, ix, iy, iz, ir) plus the whole
    (Z, r) table at the first (x, y, theta, phi) point."""
    nzpot, nr, nx, ny, ntheta, nphi = dims
    block_size = ntheta * nphi * nx * ny
    samples = sample_indices(block_size * nzpot * nr, block_size)
    out.write("# Extracted by tools/fortran_reference/extract_pes.py from the Fortran .pes file.\n")
    out.write("# V in Hartree (before the -dt premultiplication); 6D subsample (format 6d).\n")
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


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    records = read_records(sys.argv[1])
    nzpot, nr, nx, ny = unpack("<4i", next(records))
    ntheta, nphi = unpack("<2i", next(records))
    next(records)  # zgrid(1), rgrid(1), zgrid(nzpot), rgrid(nr)
    next(records)  # xgrid(1), ygrid(1), xgrid(nx), ygrid(ny)
    thetas = unpack(f"<{ntheta}d", next(records))
    phis = unpack(f"<{nphi}d", next(records))
    next(records)  # lswitch, zmin_switch, zmax_switch
    use_symmetry, n_xy_computed = unpack("<2i", next(records))
    full_size = ntheta * nphi * nx * ny
    symmetry_index = unpack(f"<{full_size}i", next(records))

    stored_size = ntheta * nphi * n_xy_computed
    values = []
    for _ in range(nzpot * nr):
        block = unpack(f"<{stored_size}d", next(records))
        values.append(expand_symmetric(block, symmetry_index) if use_symmetry else list(block))
    if next(records, None) is not None:
        sys.exit("error: unexpected trailing records")

    header = f"{nzpot} {nr} {nx} {ny} {ntheta} {nphi}\n"
    zr_only = not use_symmetry and all(v == block[0] for block in values for v in block)
    with open(sys.argv[2], "w") as out:
        if zr_only:
            write_zr_only(out, header, thetas, phis, values, nr)
        else:
            dims = (nzpot, nr, nx, ny, ntheta, nphi)
            write_6d(out, header, thetas, phis, values, dims, use_symmetry, n_xy_computed)


if __name__ == "__main__":
    main()
