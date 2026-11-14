"""Writes the inputs of the SCM reference cases (tests/scm/reference/<case>/) for the SCM driver
(scm_reference_driver.f90, run by tools/fortran_reference/run_scm_reference.sh --generate).

  example    points.in: evaluation points for the example files of the SRP48 data set
             (EXAMPLE_DIR/params.dat and examplePosfile.dat, used where they are)
  synthetic  params.in: layer 1 of the example, layer 2 changed (in the example both layers are
             equal, which would hide a mix-up of the layers); slab_v3.in and slab_v2.in: slab
             files of versions 3 and 2 with three configurations each (random displacements from a
             fixed seed), version 3 with three atoms at the cut-off of 9.5 bohr (just inside, on it,
             just outside); points.in for them

The points (u v Z r theta phi, u and v skewed lattice coordinates in bohr) exercise the correction:
the four sites at heights from below the CRP's Zmin to 6 bohr in four orientations; an r sweep
through P17 and P18 (the clamp of the parameters); H atoms 1e-3 and 1e-9 bohr inside and outside
p7 and p8 (the switch) of the layer-1 ideal atom at the origin, for r below, inside and above
[P17, P18] (at 1e-9 the switch weight rounds to exactly 0 or 1, at 1e-3 it does not); the
V_asy points (0, 0, zasymp, r, pi/2, 0); random points in the cell; and, last, 20 random points far
from the surface: 14 with Z in [8, zasymp) and 6 on the plane Z = zasymp at general (u, v), where
the correction is tiny but not always zero. Every real value of a point
list is written as its 16 hexadecimal IEEE digits, so the C++ tests read exactly the driver's input.
Only +, -, *, /, sqrt and random() are used: the output does not depend on the C library.

Usage: python3 -I -B make_scm_inputs.py example|synthetic CASE_DIR EXAMPLE_DIR
  CASE_DIR     the case directory to write (created if needed)
  EXAMPLE_DIR  potentials/H2Cu111SRP_C3v.dir (params.dat, examplePosfile.dat)
Environment: none read. Run time: under a second. Docs: tools/fortran_reference/README.md.
"""
import math
import random
import struct
import sys
from pathlib import Path

Z_ASYMPTOTIC = 12.283218925763002  # zasymp of the SRP48 data set (its GetCRPDistSurfParams value)
SYNTHETIC_V3_LATTICE = 4.95  # a_T of slab_v3.in (bohr; > a0 = 4.916 of SRP48: s > 1)
SYNTHETIC_V2_LATTICE = 4.85  # a_T of slab_v2.in
SYNTHETIC_CUTOFF = 9.5  # the cut-off of the version-3 run (run_scm_reference.sh); c*c = 90.25 exactly


def hex_line(values):
    return " ".join(struct.pack(">d", float(v)).hex() for v in values) + "\n"


def read_parameters(path):
    """The 36 values P1..P18 of layers 1 and 2 of a params.dat ("VERSION 2", then the values)."""
    values = [float(token) for line in path.read_text().splitlines()[1:] for token in line.split()]
    if len(values) != 36:
        sys.exit(f"error: {path} has {len(values)} values, not 36; name the example params.dat")
    return values[:18], values[18:]


def switch_distance(layer, k, r):
    """p_k = P(8+k)*t + P(k) at bond length r, t = r clamped to [P17, P18] (k = 7: p7, 8: p8)."""
    t = min(max(r, layer[16]), layer[17])
    return layer[7 + k] * t + layer[k - 1]


def points(lattice, layer1, seed, random_count):
    """The evaluation points (see the module docstring)."""
    pi = math.pi
    sites = [(0.0, 0.0), (lattice / 2, 0.0), (lattice / 3, lattice / 3), (2 * lattice / 3, 2 * lattice / 3)]
    orientations = [(0.0, 0.0), (pi / 2, 0.0), (pi / 2, pi / 6), (pi, 0.0)]
    result = []
    for u, v in sites:
        for z in (-1.0, 0.25, 1.0, 1.75, 2.5, 4.0, 6.0):
            result += [(u, v, z, 1.4, theta, phi) for theta, phi in orientations]
    sweep = [0.3, 0.6, layer1[16], 0.9, 1.2, 1.5, layer1[17], 2.1, 2.8, 3.5, 4.5]
    result += [(lattice / 2, 0.0, 1.5, r, pi / 2, pi / 6) for r in sweep]
    result += [(lattice / 3, lattice / 3, 2.5, r, pi / 3, 1.0) for r in sweep]
    # Upright molecule over the origin: the lower H atom at height p_k + delta above the ideal atom.
    for r in (0.6, 1.4, 2.5):
        for k in (7, 8):
            result += [(0.0, 0.0, switch_distance(layer1, k, r) + 0.5 * r + delta, r, 0.0, 0.0)
                       for delta in (-1e-3, -1e-9, 1e-9, 1e-3)]
    result += [(0.0, 0.0, Z_ASYMPTOTIC, r, pi / 2, 0.0) for r in (0.4, 0.8, 1.4, 2.0, 3.0, 4.0)]
    rng = random.Random(seed)
    for _ in range(random_count):
        result.append((lattice * rng.random(), lattice * rng.random(), -1.0 + 9.0 * rng.random(),
                       0.3 + 4.2 * rng.random(), pi * rng.random(), 2 * pi * rng.random()))
    # Far from the surface; drawn after the points above, which keep their values.
    for k in range(20):
        z = Z_ASYMPTOTIC if k >= 14 else 8.0 + (Z_ASYMPTOTIC - 8.0) * rng.random()
        result.append((lattice * rng.random(), lattice * rng.random(), z, 0.3 + 4.2 * rng.random(),
                       pi * rng.random(), 2 * pi * rng.random()))
    return f"{len(result)}\n" + "".join(hex_line(p) for p in result)


def ideal_atom(lattice, i, j, layer):
    """Ideal FCC(111) position (x, y, z) of lattice site (i, j) of a layer (1 = top, at z = 0):
    layer 2 under the hcp site (1/3, 1/3), layer 3 under the fcc site, as in the example file."""
    shift = (layer - 1) / 3.0
    u, v = (i + shift) * lattice, (j + shift) * lattice
    return (u + 0.5 * v, math.sqrt(3.0) / 2 * v, -(layer - 1) * math.sqrt(2.0 / 3.0) * lattice)


def displaced(rng, ideal, layer):
    """The ideal position plus a uniform random displacement (larger in the top layer)."""
    size = 0.3 if layer == 1 else 0.15
    return tuple(c + size * (2.0 * rng.random() - 1.0) for c in ideal)


def at_cutoff(ideal, side, c):
    """A distorted position near ideal (x, y) whose Fortran cut-off test ((x*x + y*y) + z*z < c*c)
    is just true (side -1), exactly on the boundary (0) or just false (+1): x is moved by ulps."""
    z = 0.05
    radius = math.sqrt(c * c - z * z)
    scale = radius / math.sqrt(ideal[0] * ideal[0] + ideal[1] * ideal[1])
    x0, y = ideal[0] * scale, ideal[1] * scale
    candidates = []
    x = x0
    for _ in range(64):
        x = math.nextafter(x, -math.inf)
    for _ in range(129):
        candidates.append(((x * x + y * y) + z * z, x))
        x = math.nextafter(x, math.inf)
    c2 = c * c
    if side == 0:
        on = [x for s, x in candidates if s == c2]
        if not on:
            sys.exit("error: no position exactly on the cut-off sphere; change z in at_cutoff")
        return (on[0], y, z)
    inside = [(s, x) for s, x in candidates if s < c2]
    outside = [(s, x) for s, x in candidates if s > c2]
    return (max(inside)[1] if side < 0 else min(outside)[1], y, z)


def atom_line(distorted, ideal, layer):
    """xd yd zd xi yi zi layer; xd with a Fortran D exponent (17 digits: the same double)."""
    xd = f"{distorted[0]:.16E}".replace("E", "D")
    return f"  {xd} {distorted[1]!r} {distorted[2]!r} {ideal[0]!r} {ideal[1]!r} {ideal[2]!r} {layer}\n"


def slab_v3(seed):
    """5 x 5 x 3 slab, three configurations; in each, three layer-1 atoms at the cut-off."""
    rng = random.Random(seed)
    lines = ["synthetic_v3 FCC111\n", " 5 5 3 3 3 75\n", f" {SYNTHETIC_V3_LATTICE!r}\n"]
    boundary = {(2, 0): -1, (0, 2): 0, (-2, 2): 1}  # site -> side of the cut-off
    for _ in range(3):
        for layer in (1, 2, 3):
            for j in range(-2, 3):
                for i in range(-2, 3):
                    ideal = ideal_atom(SYNTHETIC_V3_LATTICE, i, j, layer)
                    distorted = displaced(rng, ideal, layer)
                    if layer == 1 and (i, j) in boundary:
                        distorted = at_cutoff(ideal, boundary[(i, j)], SYNTHETIC_CUTOFF)
                    lines.append(atom_line(distorted, ideal, layer))
    return "".join(lines)


def slab_v2(seed):
    """Version 2: 13 atoms per configuration (layer 1: the origin and its 6 neighbours; three each
    of layers 2 and 3, which the correction skips), layers interleaved; three configurations."""
    rng = random.Random(seed)
    layer1 = [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    below = [(0, 0), (-1, 0), (0, -1)]
    order = [(1, layer1[0]), (2, below[0]), (1, layer1[1]), (3, below[0]), (1, layer1[2]), (2, below[1]),
             (1, layer1[3]), (3, below[1]), (1, layer1[4]), (2, below[2]), (1, layer1[5]), (3, below[2]),
             (1, layer1[6])]
    lines = ["synthetic_v2 FCC111\n", f" 5 5 3 3 2 {len(order)}\n", f" {SYNTHETIC_V2_LATTICE!r}\n"]
    for _ in range(3):
        for layer, (i, j) in order:
            ideal = ideal_atom(SYNTHETIC_V2_LATTICE, i, j, layer)
            lines.append(atom_line(displaced(rng, ideal, layer), ideal, layer))
    return "".join(lines)


def synthetic_parameters(layer1):
    """params.dat with layer 1 of the example (one record over two lines) and a changed layer 2."""
    layer2 = ([layer1[0] + 0.35, 0.85 * layer1[1]] + [0.6 * p for p in layer1[2:6]]
              + [layer1[6] - 0.6, layer1[7] - 0.5] + [0.5 * p for p in layer1[8:16]] + [0.8, 1.7])
    text = "VERSION 2\n"
    text += " ".join(repr(p) for p in layer1[:9]) + "\n" + " ".join(repr(p) for p in layer1[9:]) + "\n"
    return text + " ".join(repr(p) for p in layer2) + "\n", layer2


def main(argv):
    if len(argv) != 4 or argv[1] not in ("example", "synthetic"):
        print("usage: make_scm_inputs.py example|synthetic CASE_DIR EXAMPLE_DIR", file=sys.stderr)
        return 2
    case, case_dir, example_dir = argv[1], Path(argv[2]), Path(argv[3])
    layer1, _ = read_parameters(example_dir / "params.dat")
    case_dir.mkdir(parents=True, exist_ok=True)
    if case == "example":
        lattice = float((example_dir / "examplePosfile.dat").read_text().splitlines()[2].split()[0])
        (case_dir / "points.in").write_text(points(lattice, layer1, 20261114, 200))
        return 0
    params, _ = synthetic_parameters(layer1)
    (case_dir / "params.in").write_text(params)
    (case_dir / "slab_v3.in").write_text(slab_v3(3))
    (case_dir / "slab_v2.in").write_text(slab_v2(2))
    (case_dir / "points.in").write_text(points(SYNTHETIC_V3_LATTICE, layer1, 20261115, 100))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
