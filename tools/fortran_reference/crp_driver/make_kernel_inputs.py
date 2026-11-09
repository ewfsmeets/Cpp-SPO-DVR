"""Writes the inputs of the spline-kernel reference runs of the CRP driver (crp_reference_driver).

Usage: python3 -I make_kernel_inputs.py <output_dir>

Every real value is written as its 16 hexadecimal IEEE digits, so the C++ tests read exactly the
inputs the Fortran driver got. Deterministic (no randomness).
"""
import math
import os
import struct
import sys


def h(x):
    return struct.pack(">d", float(x)).hex()


def hex_line(values):
    return " ".join(h(v) for v in values) + "\n"


def spline1d(x, y, ic1, cm1, icn, cmn, t):
    return f"{len(x)} {ic1} {icn} {len(t)}\n" + hex_line([cm1, cmn]) + hex_line(x) + hex_line(y) + hex_line(t)


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)

    def write(name, text):
        with open(os.path.join(out, name), "w") as f:
            f.write(text)

    # 1D splines: non-uniform knots, every end condition, queries outside, at and between knots.
    x = [0.0, 0.3, 0.75, 1.2, 1.9, 2.5, 3.1, 4.0]
    y = [math.sin(1.3 * v) + 0.2 * v * v for v in x]
    t = [-0.5, 0.0, 0.1, 0.3, 0.6, 1.2, 1.55, 2.2, 3.1, 3.7, 4.0, 4.6]
    write("spline1d_not_a_knot.in", spline1d(x, y, 0, 0.0, 0, 0.0, t))
    write("spline1d_clamped.in", spline1d(x, y, 2, 1.1, 2, -0.7, t))
    write("spline1d_zdep.in", spline1d(x, y, 0, 0.0, 2, 0.0, t))  # as Molecule/extrapol.dat
    write("spline1d_second_derivative.in", spline1d(x, y, 1, 0.4, 1, -0.3, t))
    write("spline1d_natural.in", spline1d(x, y, 3, 9.0, 5, 9.0, t))
    write("spline1d_decreasing.in", spline1d(x[::-1], y[::-1], 0, 0.0, 0, 0.0, t))

    # Bicubic spline (r, Z cut): queries inside, on the edges and outside.
    gx = [0.5, 0.8, 1.0, 1.4, 1.7, 2.3, 3.0]
    gy = [-0.5, 0.25, 1.0, 2.0, 3.5, 5.0]
    f2 = [math.cos(a) * math.exp(-0.3 * b) + 0.1 * a * b for b in gy for a in gx]  # x fastest
    q2 = [(0.6, 0.0), (1.0, 1.0), (2.9, 4.9), (0.5, -0.5), (3.0, 5.0), (1.33, 2.71), (0.4, 1.0), (1.0, 5.2)]
    text = f"{len(gx)} {len(gy)} {len(q2)}\n" + hex_line(gx) + hex_line(gy) + hex_line(f2)
    text += "".join(hex_line(q) for q in q2)
    write("spline2d.in", text)

    # Periodic 3D spline (C4v atom surface): 5 x 5 cell of period 4.8, clamped in z.
    a = 4.8
    px = [i * a * 0.25 for i in range(5)]
    pz = [-1.0, -0.4, 0.3, 1.0, 1.8, 2.7, 3.9, 5.0]
    f3 = []
    for k, z in enumerate(pz):
        for j in range(5):
            for i in range(5):
                ii, jj = i % 4, j % 4  # periodic: the last row/column repeats the first
                u, v = ii * 0.25, jj * 0.25
                f3.append(math.exp(-0.9 * z) * (1 + 0.3 * math.cos(2 * math.pi * u) * math.cos(2 * math.pi * v))
                          + 0.05 * math.sin(2 * math.pi * u) * z)
    d1 = [-0.9 * (1 + 0.1 * ((i % 4) + (j % 4))) for j in range(5) for i in range(5)]
    d2 = [0.002 * (1 + (i % 4) * (j % 4)) for j in range(5) for i in range(5)]
    q3 = [(0.1, 0.2, 0.0), (2.4, 1.1, -0.7), (4.8, 4.8, 5.0), (0.0, 0.0, -1.0), (3.3, 0.7, 2.2), (1.9, 3.9, 4.4)]
    text = f"5 5 {len(pz)} {len(q3)}\n" + hex_line(px) + hex_line(px) + hex_line(pz) + hex_line(f3)
    text += hex_line(d1) + hex_line(d2) + "".join(hex_line(q) for q in q3)
    write("spline3d.in", text)

    # Periodic 1D spline (TRISPL): non-uniform period, y[n-1] = y[0].
    tx = [0.0, 0.4, 1.1, 1.5, 2.2, 2.9, 3.3, 4.1, 4.8]
    ty = [math.cos(2 * math.pi * v / 4.8) + 0.3 * math.sin(4 * math.pi * v / 4.8) for v in tx[:-1]]
    ty.append(ty[0])
    write("trispl.in", f"{len(tx)}\n" + hex_line(tx) + hex_line(ty))

    # UBICAR: below, at and between the knots, above.
    ux = [0.5, 1.0, 1.25, 2.0, 3.5, 5.0]
    ut = [0.4, 0.5, 0.75, 1.0, 1.1, 3.5, 4.99, 5.0, 5.01]
    write("ubicar.in", f"{len(ux)} {len(ut)}\n" + hex_line(ux) + hex_line(ut))


if __name__ == "__main__":
    main()
