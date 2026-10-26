#!/usr/bin/env python3
"""Extract a compact text reference from a Fortran SPO-DVR binary .pes file.

The .pes file is written by Potential.f90 WritePesToFile as gfortran unformatted sequential
records (4-byte record markers, default INTEGER/LOGICAL = 4 bytes, REAL = 8 bytes because of
-fdefault-real-8). The stored potential is V in Hartree *before* premultiplication by -dt,
layout vdvr(itp, ixy, izr) with itp = iphi + nphi*(itheta-1), izr = ir + nr*(iz-1).

Output (text, full precision):
  header line with nzpot nr nx ny ntheta nphi
  the theta and phi grid points
  for a potential that depends only on (Z, r) (the hard wall): V(iz, ir) for every (Z, r)
  point; the script verifies that each (Z, r) block is constant over (x, y, theta, phi) and
  aborts otherwise.

Usage: extract_pes.py <file.pes> <output.txt>
"""
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
    use_symmetry, _ = unpack("<2i", next(records))
    next(records)  # symmetry indices
    if use_symmetry:
        sys.exit("error: surface symmetry was used; the (x,y) layout is not supported here")

    block_size = ntheta * nphi * nx * ny
    values = []
    for izr in range(nzpot * nr):
        block = unpack(f"<{block_size}d", next(records))
        if any(v != block[0] for v in block):
            sys.exit(f"error: potential block izr = {izr + 1} is not constant over (x,y,theta,phi)")
        values.append(block[0])
    if next(records, None) is not None:
        sys.exit("error: unexpected trailing records")

    with open(sys.argv[2], "w") as out:
        out.write("# Extracted by tools/fortran_reference/extract_pes.py from the Fortran .pes file.\n")
        out.write("# V in Hartree (before the -dt premultiplication); constant over (x, y, theta, phi).\n")
        out.write("# nzpot nr nx ny ntheta nphi\n")
        out.write(f"{nzpot} {nr} {nx} {ny} {ntheta} {nphi}\n")
        out.write("# theta points\n")
        out.write("".join(f"{t:.17e}\n" for t in thetas))
        out.write("# phi points\n")
        out.write("".join(f"{p:.17e}\n" for p in phis))
        out.write("# iz ir V(iz,ir)  (1-based indices, ir fastest)\n")
        for izr, v in enumerate(values):
            iz, ir = divmod(izr, nr)
            out.write(f"{iz + 1} {ir + 1} {v:.17e}\n")


if __name__ == "__main__":
    main()
