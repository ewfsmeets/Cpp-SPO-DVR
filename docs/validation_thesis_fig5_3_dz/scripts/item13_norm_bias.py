#!/usr/bin/env python3
"""Item 13/14: a candidate for the uniform amplitude difference between the C++ and Fortran
propagations (C++ S-matrix elements larger by ~6e-12 in every channel, item18b).

Every kinetic block T(dt') does an FFT round trip (coordinate -> momentum -> coordinate) per
grid direction, with explicit scale factors s (forward) and 1/(s N) (inverse). In exact
arithmetic N s (1/(s N)) = 1; with the double-precision factors each round trip multiplies the
wave by 1 + delta, the same delta at every step. delta differs between the codes: the Fortran
transforms Z, r, x, y separately (four factors per direction pair), computes sqrt(dx/nx) in
quad precision (scalFacX, a documented quad-precision site), while the C++ transforms (Z, r) and
(x, y) as 2D FFTs with one combined forward and one combined inverse factor. This script
computes delta of both codes exactly (fractions of the double-precision factors, IEEE
semantics, i.e. without the fast-math reassociation of the C++ build) and the predicted
amplitude ratio after the 7875 kinetic blocks of the run (9 per analysis step: T/2 V (T V)^7
T/2). Writes data/items/item13_norm_bias.json."""
from __future__ import annotations

import json
import math
from fractions import Fraction as Fr

from common import INPUT_DIR, load_item, save_item


def exact_round_trip(n, fwd, inv):
    return Fr(n) * Fr(fwd) * Fr(inv) - 1


def quad_sqrt_ratio(a, n):
    """sqrt(a / n) with a, n exact, correctly rounded to double (Fortran quad: dfloat)."""
    q = Fr(a) / n
    # integer square root on a scaled fraction, then exact rounding to double
    x = math.sqrt(float(q))
    # refine: pick the double among neighbours whose square is closest to q
    cands = [x, math.nextafter(x, 0), math.nextafter(x, 2 * x + 1)]
    return min(cands, key=lambda c: abs(Fr(c) * Fr(c) - q))


def main():
    js = json.loads((INPUT_DIR / "input.json").read_text())
    nz, nzsp, nr, nx, ny = js["nZ"], js["nZSpecular"], js["nR"], js["nX"], js["nY"]
    dz, dr = js["dZ"], js["dR"]
    alat = 4.84
    dx = alat / nx
    dy = alat / ny
    pi = 3.141592653589793
    sq2pi = math.sqrt(2.0 * pi)
    sz, sr = dz / sq2pi, dr / sq2pi
    # C++ (src/propagation/propagation_tables.cpp, kinetic_operators.cpp, specular_operator.cpp), IEEE double
    sx_c, sy_c = math.sqrt(dx / nx), math.sqrt(dy / ny)
    cpp = {
        "zr": exact_round_trip(nz * nr, sz * sr, (1.0 / (nz * sz)) * (1.0 / (nr * sr))),
        "xy": exact_round_trip(nx * ny, sx_c * sy_c, (1.0 / (nx * sx_c)) * (1.0 / (ny * sy_c))),
        "specular": exact_round_trip(nzsp * nr, sz * sr, (1.0 / (nzsp * sz)) * (1.0 / (nr * sr))),
    }
    # Fortran (KineticOp.f90 SetupScaleFactors, ffts.f90 InverseComplexFFT): 1D transforms
    sx_f, sy_f = quad_sqrt_ratio(dx, nx), quad_sqrt_ratio(dy, ny)
    f1 = {d: exact_round_trip(n, s, 1.0 / (s * float(n))) for d, (n, s) in
          {"z": (nz, sz), "r": (nr, sr), "x": (nx, sx_f), "y": (ny, sy_f), "zsp": (nzsp, sz)}.items()}
    fortran = {"zr": (1 + f1["z"]) * (1 + f1["r"]) - 1, "xy": (1 + f1["x"]) * (1 + f1["y"]) - 1,
               "specular": (1 + f1["zsp"]) * (1 + f1["r"]) - 1}
    n_blocks = 875 * 9
    out = {"n_kinetic_blocks": n_blocks, "scale_factors": {"sz": sz, "sr": sr, "sx_cpp": sx_c, "sx_fortran_quad": sx_f,
                                                           "sx_ulps_apart": abs(sx_c - sx_f) / math.ulp(sx_c)},
           "delta_per_block": {}, "predicted_amplitude_ratio_minus_1": {}}
    for g in ("zr", "xy", "specular"):
        out["delta_per_block"][g] = {"cpp": float(cpp[g]), "fortran": float(fortran[g])}
    scat_c = (1 + cpp["zr"]) * (1 + cpp["xy"])
    scat_f = (1 + fortran["zr"]) * (1 + fortran["xy"])
    out["predicted_amplitude_ratio_minus_1"]["scattering_grid"] = float((scat_c / scat_f) ** 1 - 1) * n_blocks
    out["predicted_amplitude_ratio_minus_1"]["specular_grid"] = float(((1 + cpp["specular"]) / (1 + fortran["specular"])) - 1) * n_blocks
    # observed (item18b): alpha from the replays, per energy
    out["note"] = ("IEEE double semantics for the C++ factors; the run used -ffast-math, which may evaluate "
                   "divisions as multiplications by reciprocals (different last bits)")
    save_item("item13_norm_bias", out)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
