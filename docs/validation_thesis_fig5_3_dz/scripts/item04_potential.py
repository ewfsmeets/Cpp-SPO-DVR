#!/usr/bin/env python3
"""Item 4 (and the grid part of item 2): the 6D potential of the three programs on the full grid.

Inputs (extract mode, heavy): the Fortran .pes of the run (V in hartree before the -dt factor,
symmetry-reduced; layout: tools/fortran_reference/extract_pes.py) and the outputs of the helper
`potential_check` (C++ library of the build that ran): v_stored.bin (V as the run stores it, at
the C++ representatives), symmetry_maps.txt (the C++ gather), v_full.bin (the PES at every
full-grid point, no symmetry), grids.txt.

Comparisons (scaled |dV| = |dV| / max(|V|, 1 eV), V in hartree, 1 eV = 1/27.2114 hartree):
  same_point : Fortran stored value vs C++ PES at the same physical point (the Fortran's own
               representative points: identical coordinates), all 85.5 million stored values
  used       : what each program uses at every full-grid point (Fortran: expanded through
               symmetryIndex; C++: gathered through its maps), 648 million values
  asym_cpp / asym_fortran : the PES's own asymmetry, C++ PES at a point vs at the
               representative image chosen by the C++ / Fortran maps
Writes data/items/item04.json."""
from __future__ import annotations

import hashlib
import struct
import sys

import numpy as np

from common import CHAEV, RUN_DIRS, WORK, fortran_record_offsets, read_record, save_item, stats

ONE_EV = 1.0 / CHAEV
TWO_EV = 2.0 / CHAEV


def read_pes(path):
    recs = fortran_record_offsets(path)
    h = [read_record(path, o, n) for o, n in recs[:9]]
    nzpot, nr, nx, ny = struct.unpack("<4i", h[0])
    ntheta, nphi = struct.unpack("<2i", h[1])
    z1, r1, zn, rn = struct.unpack("<4d", h[2])
    x1, y1, xn, yn = struct.unpack("<4d", h[3])
    thetas = np.frombuffer(h[4], dtype="<f8")
    phis = np.frombuffer(h[5], dtype="<f8")
    use_symmetry, n_xy = struct.unpack("<2i", h[7])
    symidx = np.frombuffer(h[8], dtype="<i4").astype(np.int64) - 1
    blocks = recs[9:]
    assert len(blocks) == nzpot * nr
    stored = ntheta * nphi * n_xy
    assert {n for _, n in blocks} == {8 * stored}
    stride = blocks[1][0] - blocks[0][0]
    # strided float64 view over the record payloads
    base = np.memmap(path, dtype="<f8", mode="r", offset=blocks[0][0],
                     shape=((len(blocks) - 1) * stride // 8 + stored,))
    values = np.lib.stride_tricks.as_strided(base, shape=(len(blocks), stored), strides=(stride, 8))
    return {"nzpot": nzpot, "nr": nr, "nx": nx, "ny": ny, "ntheta": ntheta, "nphi": nphi,
            "grid_ends": {"z1": z1, "r1": r1, "zn": zn, "rn": rn, "x1": x1, "y1": y1, "xn": xn, "yn": yn},
            "thetas": thetas, "phis": phis, "use_symmetry": use_symmetry, "n_xy_computed": n_xy,
            "symidx": symidx, "values": values}


def read_grids(path):
    out = {}
    lines = open(path).read().split("\n")
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if len(t) == 2 and not t[0].replace("_", "").isdigit():
            name, n = t[0], int(t[1]) if t[1].isdigit() else None
            if n is not None and name not in ("nz_potential", "iz_analysis"):
                out[name] = np.array([float(x) for x in lines[i + 1:i + 1 + n]])
                i += n + 1
                continue
            out[name] = int(t[1])
        i += 1
    return out


def read_symmetry(path, ntp, nphi):
    lines = open(path).read().splitlines()
    n_xy = int(lines[0].split()[1])
    nrep = int(lines[0].split()[3])
    rep_of = np.zeros(n_xy, dtype=np.int64)
    op_of = np.zeros(n_xy, dtype=np.int64)
    rep_xy = np.zeros(nrep, dtype=np.int64)
    phimap = {}
    for line in lines[1:]:
        t = line.split()
        if t[0] == "rep":
            rep_xy[int(t[1])] = int(t[2])
        elif t[0] == "phimap":
            phimap[int(t[1])] = np.array([int(x) for x in t[2:]])
        else:
            rep_of[int(t[0])], op_of[int(t[0])] = int(t[1]), int(t[2])
    ntheta = ntp // nphi
    itp = np.arange(ntp)
    iphi, itheta = itp % nphi, itp // nphi
    gather = np.empty(ntp * n_xy, dtype=np.int64)       # full (itp, ixy) -> stored (itp', irep)
    for ixy in range(n_xy):
        image = phimap[int(op_of[ixy])][iphi] + nphi * itheta
        gather[ixy * ntp:(ixy + 1) * ntp] = image + ntp * rep_of[ixy]
    stored_to_full = np.empty(ntp * nrep, dtype=np.int64)  # stored (itp', irep) -> full (itp', rep_xy)
    for irep in range(nrep):
        stored_to_full[irep * ntp:(irep + 1) * ntp] = itp + ntp * rep_xy[irep]
    return {"n_xy": n_xy, "n_rep": nrep, "gather": gather, "stored_to_full": stored_to_full,
            "rep_xy": rep_xy.tolist(), "n_operations_used": len(set(op_of.tolist()))}


def fortran_representatives(symidx, ntp, n_xy_full, n_xy_red):
    """Full-grid ixy of each Fortran stored (x, y) point: the point whose block maps onto it by the
    identity (symmetryIndex(itp, ixy) = (itp, k) for every itp)."""
    reps = {}
    itp = np.arange(ntp)
    for ixy in range(n_xy_full):
        s = symidx[ixy * ntp:(ixy + 1) * ntp]
        k = s[0] // ntp
        if np.array_equal(s, itp + ntp * k):
            reps.setdefault(int(k), ixy)
    assert sorted(reps) == list(range(n_xy_red)), reps
    rep = np.array([reps[k] for k in range(n_xy_red)])
    stored_to_full = (np.arange(ntp)[None, :] + ntp * rep[:, None]).reshape(-1)
    return rep, stored_to_full


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def describe(a, where_fn, mask=None):
    """n, median, p95, max (+ location) of a; partitions a in place (a is scratch afterwards)."""
    if mask is not None:
        idx = np.flatnonzero(mask)
        a = a[idx]
    else:
        idx = None
    if a.size == 0:
        return {"n": 0}
    i = int(np.argmax(a))
    out = {"n": int(a.size), "max": float(a[i]), "argmax": where_fn(int(idx[i]) if idx is not None else i),
           "n_zero": int(np.count_nonzero(a == 0))}
    k50, k95 = (a.size - 1) // 2, int(round(0.95 * (a.size - 1)))
    a.partition([k50, k95])
    out["median"] = float(a[k50])
    out["p95"] = float(a[k95])
    return out


def main():
    build = sys.argv[1] if len(sys.argv) > 1 else "build-fast"
    hdir = WORK / "potential" / build
    pes_path = next(RUN_DIRS["fortran"].glob("*.pes"))
    F = read_pes(pes_path)
    g = read_grids(hdir / "grids.txt")
    nzpot, nr, nx, ny, ntheta, nphi = (F[k] for k in ("nzpot", "nr", "nx", "ny", "ntheta", "nphi"))
    ntp, nxy, nzr = ntheta * nphi, nx * ny, nzpot * nr
    sym = read_symmetry(hdir / "symmetry_maps.txt", ntp, nphi)
    c_stored = np.memmap(hdir / "v_stored.bin", dtype="<f8", mode="r").reshape(nzr, ntp * sym["n_rep"])
    c_full = np.memmap(hdir / "v_full.bin", dtype="<f8", mode="r").reshape(nzr, ntp * nxy)
    f_rep, f_stored_to_full = fortran_representatives(F["symidx"], ntp, nxy, F["n_xy_computed"])

    res = {"build": build, "pes_file": pes_path.name,
           "dims": {k: F[k] for k in ("nzpot", "nr", "nx", "ny", "ntheta", "nphi", "n_xy_computed", "use_symmetry")},
           "cpp_n_representatives": sym["n_rep"], "cpp_representatives_xy": sym["rep_xy"],
           "fortran_representatives_xy": f_rep.tolist(),
           "n_common_representatives": len(set(sym["rep_xy"]) & set(f_rep.tolist()))}

    # grids at full precision: .pes header records vs C++ grids.txt (item 2)
    ge = F["grid_ends"]
    cg = {"z1": g["z_scattering"][0], "zn": g["z_scattering"][nzpot - 1], "r1": g["r"][0], "rn": g["r"][-1],
          "x1": g["x"][0], "xn": g["x"][-1], "y1": g["y"][0], "yn": g["y"][-1]}
    res["grid_ends"] = {k: {"fortran": ge[k], "cpp": cg[k], "diff": cg[k] - ge[k],
                            "ulps": abs(cg[k] - ge[k]) / np.spacing(abs(ge[k])) if ge[k] != 0 else (0.0 if cg[k] == 0 else float("inf"))}
                        for k in ge}
    dth = g["theta"] - F["thetas"]
    dph = g["phi"] - F["phis"]
    res["theta"] = {"n": len(dth), "max_abs_diff": float(np.max(np.abs(dth))), "n_identical": int(np.sum(dth == 0)),
                    "max_ulps": float(np.max(np.abs(dth) / np.spacing(np.abs(F["thetas"]))))}
    res["phi"] = {"n": len(dph), "max_abs_diff": float(np.max(np.abs(dph))), "n_identical": int(np.sum(dph == 0))}

    n_used = nzr * ntp * nxy
    n_same = nzr * ntp * F["n_xy_computed"]
    symidx = F["symidx"]
    gather = sym["gather"]
    f_image = f_stored_to_full[symidx]          # full point -> Fortran representative image (full index)
    c_image = sym["stored_to_full"][gather]     # full point -> C++ representative image (full index)
    below_2ev = np.empty(n_used, dtype=np.bool_)
    same_below_2ev = np.empty(n_same, dtype=np.bool_)

    def scan(kind):
        """One pass over the (Z, r) blocks; returns the scaled differences of one comparison
        (one large array in memory at a time)."""
        n_out = n_same if kind == "same" else n_used
        out = np.empty(n_out, dtype=np.float64 if kind == "same" else np.float32)
        max_abs = 0.0
        for izr in range(nzr):
            cf = np.asarray(c_full[izr])
            lo, hi = izr * ntp * nxy, (izr + 1) * ntp * nxy
            if kind == "same":
                fs = np.asarray(F["values"][izr])
                lo2, hi2 = izr * ntp * F["n_xy_computed"], (izr + 1) * ntp * F["n_xy_computed"]
                out[lo2:hi2] = np.abs(fs - cf[f_stored_to_full]) / np.maximum(np.abs(fs), ONE_EV)
                same_below_2ev[lo2:hi2] = fs < TWO_EV
            elif kind == "used":
                f_exp = np.asarray(F["values"][izr])[symidx]
                c_sym = np.asarray(c_stored[izr])[gather]
                d = np.abs(f_exp - c_sym)
                max_abs = max(max_abs, float(d.max()))
                out[lo:hi] = d / np.maximum(np.abs(f_exp), ONE_EV)
                below_2ev[lo:hi] = cf < TWO_EV
            elif kind == "asym_cpp":
                out[lo:hi] = np.abs(cf - cf[c_image]) / np.maximum(np.abs(cf), ONE_EV)
            else:
                out[lo:hi] = np.abs(cf - cf[f_image]) / np.maximum(np.abs(cf), ONE_EV)
        return out, max_abs

    def where_full(i):
        izr, rest = divmod(i, ntp * nxy)
        ixy, itp = divmod(rest, ntp)
        iz, ir = divmod(izr, nr)
        ix, iy = divmod(ixy, ny)
        itheta, iphi = divmod(itp, nphi)
        return {"iz": iz, "ir": ir, "ix": ix, "iy": iy, "itheta": itheta, "iphi": iphi,
                "Z": float(g["z_scattering"][iz]), "r": float(g["r"][ir]),
                "V_eV": float(c_full[izr, rest]) * CHAEV}

    def where_same(i):
        izr, rest = divmod(i, ntp * F["n_xy_computed"])
        k, itp = divmod(rest, ntp)
        return where_full(izr * ntp * nxy + int(f_rep[k]) * ntp + itp)

    # each describe() partitions its argument: masked statistics first
    same, _ = scan("same")
    res["same_point_v_below_2eV"] = describe(same, where_same, mask=same_below_2ev)
    res["same_point"] = describe(same, where_same)
    del same
    used, res["used_max_abs_hartree"] = scan("used")
    res["used_v_below_2eV"] = describe(used, where_full, mask=below_2ev)
    res["used"] = describe(used, where_full)
    del used
    asym_c, _ = scan("asym_cpp")
    res["asym_cpp_maps_v_below_2eV"] = describe(asym_c, where_full, mask=below_2ev)
    res["asym_cpp_maps"] = describe(asym_c, where_full)
    del asym_c
    asym_f, _ = scan("asym_fortran")
    res["asym_fortran_maps"] = describe(asym_f, where_full)
    del asym_f
    res["n_points_below_2eV"] = int(np.count_nonzero(below_2ev))
    vmin, vmax = np.inf, -np.inf
    for izr in range(nzr):
        cf = np.asarray(c_full[izr])
        vmin, vmax = min(vmin, float(cf.min())), max(vmax, float(cf.max()))
    res["v_range_eV"] = [vmin * CHAEV, vmax * CHAEV]
    res["criterion"] = 1e-13
    res["sha256"] = {"v_stored.bin": sha256(hdir / "v_stored.bin")}
    save_item("item04" if build == "build-fast" else f"item04_{build}", res)
    print("same point :", {k: res["same_point"][k] for k in ("n", "median", "p95", "max")}, res["same_point"]["argmax"])
    print("same <2 eV :", {k: res["same_point_v_below_2eV"][k] for k in ("n", "median", "p95", "max")})
    print("used       :", {k: res["used"][k] for k in ("n", "median", "p95", "max")}, res["used"]["argmax"])
    print("asym C++   :", {k: res["asym_cpp_maps"][k] for k in ("median", "p95", "max")}, res["asym_cpp_maps"]["argmax"])
    print("asym F     :", {k: res["asym_fortran_maps"][k] for k in ("median", "p95", "max")})
    print("reps C++", sym["rep_xy"], "\nreps F  ", f_rep.tolist(), "common", res["n_common_representatives"])
    print("theta", res["theta"], "phi", res["phi"])
    for k, v in res["grid_ends"].items():
        print(" ", k, v)


if __name__ == "__main__":
    main()
