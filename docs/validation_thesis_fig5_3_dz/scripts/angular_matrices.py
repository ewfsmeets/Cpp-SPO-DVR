#!/usr/bin/env python3
"""Angular matrices (criteria.md, 2026-11-08): does the codes' different angular round trip
explain the specular rotationally inelastic family excess?

The potential operator takes every (Z, r) point with Z < Zana through FBR -> DVR ->
exp(-iV dt) -> DVR -> FBR, after adding the specular wave to j = 0, (n, m) = (0, 0). What the round
trip leaks from j = 0 into (J, mJ = 0) stays on the scattering grid, in the family's channels.

Inputs (runs/thesis_fig5_3_dz_validation/angular/<code>/, make_all.sh --heavy):
  fortran          angular_fortran.f90 linked against the reference executable's object files
  build-cuda       angular_check.cpp linked against build-cuda (the strict GPU run)
  build-fast       ... build-fast (the fast CPU run)
  build-cuda-fast  ... build-cuda-fast (the fast GPU run)
each with thetas.bin, weights.bin, matrix.bin, basis.txt, roundtrip.bin; the four BK files; the
C++ potential on the full grid (item 4, potential/build-fast/v_full.bin) for the robustness check.
Writes data/items/angular_matrices.json."""
from __future__ import annotations

import json
import re
from fractions import Fraction

import numpy as np

from common import CHAEV, INPUT_DIR, RUN_DIRS, RUNS, WORK, FortranBk, NativeBk, save_item
from item04_potential import read_pes
from parse_logs import cpp_section

ANG = WORK / "angular"
CODES = {"fortran": ANG / "fortran", "strict": ANG / "build-cuda", "fast": ANG / "build-fast",
         "fast_gpu": ANG / "build-cuda-fast"}
STRICT = RUNS / "thesis_fig5_3_dz_cuda"
BK_FILES = {"fortran": RUN_DIRS["fortran"] / "fort.21", "strict": STRICT / "bk_coefficients.bin",
            "fast": RUN_DIRS["cpu"] / "bk_coefficients.bin", "fast_gpu": RUN_DIRS["gpu"] / "bk_coefficients.bin"}
LOGS = {"strict": WORK / "logs" / "gpu_strict.txt", "fast": WORK / "logs" / "cpu.txt", "fast_gpu": WORK / "logs" / "gpu.txt"}
PAIRS = [("strict", "fortran"), ("fast", "fortran"), ("strict", "fast"), ("fast_gpu", "fortran"), ("fast_gpu", "fast")]
VERDICT_J = (8, 10)


def pair_key(a, b):
    return f"{a}-{b}"


# ----------------------------------------------------------------------------- the dumps
def load_code(path):
    head, *rows = (path / "basis.txt").read_text().split("\n")
    n_rot, nphi, ntheta = (int(x) for x in head.split())
    basis = np.array([[int(x) for x in r.split()] for r in rows if r.strip()])
    assert basis.shape == (n_rot, 2)
    m0 = np.flatnonzero(basis[:, 1] == 0)
    rt = np.fromfile(path / "roundtrip.bin", dtype="<c16").reshape(len(m0), n_rot)
    matrix = np.fromfile(path / "matrix.bin", dtype="<f8")
    # matrix order (mJ block, theta, J): the mJ = 0 block follows the blocks of mJ < 0
    start = sum(ntheta * int(np.sum(basis[:, 1] == m)) for m in range(int(basis[:, 1].min()), 0))
    m0_block = matrix[start:start + ntheta * len(m0)].reshape(ntheta, len(m0))
    return {"n_rot": n_rot, "nphi": nphi, "ntheta": ntheta, "basis": basis, "m0": m0, "J": basis[m0, 0],
            "thetas": np.fromfile(path / "thetas.bin", dtype="<f8"), "weights": np.fromfile(path / "weights.bin", dtype="<f8"),
            "matrix": matrix, "m0_block": m0_block, "roundtrip": rt}


def ulp_compare(a, b):
    d = np.abs(a - b)
    return {"n": int(a.size), "n_identical": int(np.sum(a == b)), "max_abs": float(d.max()),
            "max_rel": float(np.max(d / np.abs(b), where=b != 0, initial=0.0)),
            "max_ulps": float(np.max(d / np.spacing(np.abs(b)), where=b != 0, initial=0.0))}


def round_trip_leak(c):
    """R_J: the (J, mJ = 0) components of the round trip of e_(J=0), minus e_(J=0)."""
    out = c["roundtrip"][0].copy()
    out[c["m0"][0]] -= 1.0
    leak_m0 = out[c["m0"]]
    other = np.delete(out, c["m0"])
    return leak_m0, float(np.abs(other).max())


def exact_products(c):
    """M_t0 M_tJ of the mJ = 0 block as exact rationals (the stored doubles), per (theta, J)."""
    M = c["m0_block"]
    col0 = [Fraction(float(x)) for x in M[:, 0]]
    return [[Fraction(float(M[t, j])) * col0[t] for j in range(M.shape[1])] for t in range(M.shape[0])]


def matrix_leak(c, prod):
    """R_J from the matrix values alone: nphi sum_t M_tJ M_t0 - delta_J0, exact, then rounded once."""
    nphi = c["nphi"]
    out = []
    for j in range(len(c["J"])):
        s = nphi * sum(prod[t][j] for t in range(len(prod))) - (1 if j == 0 else 0)
        out.append(float(s))
    return np.array(out)


def full_round_trip(c):
    """Every mJ = 0 basis function through the round trip: the largest deviation from the identity."""
    eye = np.zeros_like(c["roundtrip"])
    for k, i in enumerate(c["m0"]):
        eye[k, i] = 1.0
    d = np.abs(c["roundtrip"] - eye)
    off_m0 = np.delete(d, c["m0"], axis=1)
    return {"max_abs_deviation": float(d.max()), "max_outside_mJ0": float(off_m0.max()),
            "max_imag": float(np.abs(c["roundtrip"].imag).max())}


# ----------------------------------------------------------------------------- validity
def log_theta_table(path):
    lines = cpp_section(path.read_text(encoding="utf-8").splitlines(), "GRIDS")
    rows = [[float(x) for x in l.split()] for l in lines if re.match(r"^\s+\d+\s+\S+e[+-]\d+\s+\S+e[+-]\d+$", l)]
    return np.array([r[1] for r in rows]), np.array([r[2] for r in rows])


def validity(codes):
    """The dumped tables are those of the runs: Fortran theta nodes vs the run's .pes header; C++
    nodes and weights vs the run log's theta table (printed with 17 digits)."""
    pes = next(RUN_DIRS["fortran"].glob("*.pes"))
    out = {"fortran_thetas_equal_pes_header": bool(np.array_equal(codes["fortran"]["thetas"], read_pes(pes)["thetas"]))}
    for k, path in LOGS.items():
        th, w = log_theta_table(path)
        out[f"{k}_weights_equal_run_log"] = bool(np.array_equal(w, codes[k]["weights"]))
        out[f"{k}_thetas_vs_run_log"] = ulp_compare(codes[k]["thetas"], th)
    out["basis_identical"] = all(np.array_equal(codes[k]["basis"], codes["fortran"]["basis"]) for k in codes)
    out["fast_matrix_equals_fast_gpu_matrix"] = bool(np.array_equal(codes["fast"]["matrix"], codes["fast_gpu"]["matrix"]))
    out["strict_columns_identical"] = (CODES["strict"] / "columns_identical.txt").read_text().strip() == "1"
    # the strict and Fortran tables must be exactly the runs'; the fast-math helper's machine code
    # is generated at its own link (as the replay helper's), so it may differ by an ulp
    out["valid"] = (out["fortran_thetas_equal_pes_header"] and out["strict_weights_equal_run_log"]
                    and out["strict_thetas_vs_run_log"]["n_identical"] == out["strict_thetas_vs_run_log"]["n"]
                    and out["fast_weights_equal_run_log"] and out["basis_identical"] and out["strict_columns_identical"])
    return out


def fast_theta_sensitivity(codes, J_list):
    """Change of nphi M_tJ M_t0 if the fast helper's theta nodes were the fast run's (log table):
    w_t d/dtheta [P~_J P~_0] (theta_run - theta_helper), P~_J = sqrt(J + 1/2) P_J(cos theta)."""
    th_run, _ = log_theta_table(LOGS["fast"])
    c = codes["fast"]
    dth = th_run - c["thetas"]
    out = np.zeros(len(J_list))
    for i, j in enumerate(J_list):
        cj = np.zeros(j + 1); cj[j] = np.sqrt(j + 0.5)
        x = np.cos(c["thetas"])
        dPdx = np.polynomial.legendre.legval(x, np.polynomial.legendre.legder(cj))
        out[i] = np.sum(c["weights"] * np.sqrt(0.5) * (-np.sin(c["thetas"])) * dPdx * dth)
    return out, int(np.sum(dth != 0))


# ----------------------------------------------------------------------------- the BK family
def family_bk(J_list):
    readers = {k: (FortranBk(p) if k == "fortran" else NativeBk(p)) for k, p in BK_FILES.items()}
    n = readers["fortran"].n_records
    assert all(r.n_records == n for r in readers.values())
    lab = readers["fast"].labels
    idx = [int(np.flatnonzero((lab["J"] == J) & (lab["mJ"] == 0) & (lab["n"] == 0) & (lab["m"] == 0))[0]) for J in J_list]
    d = {pair_key(a, b): np.zeros((n, len(J_list)), dtype=complex) for a, b in PAIRS}
    raw = {pair_key(a, b): np.zeros((n, len(J_list)), dtype=complex) for a, b in PAIRS}
    cF = np.zeros((n, len(J_list)), dtype=complex)
    l2 = np.zeros(n)
    for k in range(n):
        rec = {c: r.record(k) for c, r in readers.items()}
        l2[k] = np.sqrt(np.sum(np.abs(rec["fortran"]) ** 2))
        cF[k] = rec["fortran"][0, idx]
        for a, b in PAIRS:
            diff = rec[a] - rec[b]
            den = np.sum(np.abs(rec[b]) ** 2)
            z = np.sum(np.conj(rec[b]) * diff) / den if den > 0 else 0.0
            raw[pair_key(a, b)][k] = diff[0, idx]
            d[pair_key(a, b)][k] = diff[0, idx] - z * rec[b][0, idx]
    window = l2 >= 1e-2 * l2.max()
    return {"d": d, "raw": raw, "cF": cF, "window": window, "times": readers["fast"].times()}


# ----------------------------------------------------------------------------- robustness: exp(-iV dt)
def potential_weights():
    """Mean over x, y and phi of exp(-iV dt) per (theta, Z, r) of the interaction region, the C++
    potential (item 4), and the selection V_mean < min V_asy + 1 eV."""
    pdir = WORK / "potential" / "build-fast"
    txt = (pdir / "grids.txt").read_text().split("\n")
    nzpot = int([l for l in txt if l.startswith("nz_potential")][0].split()[1])
    js = json.loads((INPUT_DIR / "input.json").read_text())
    nr, nxy, dt = js["nR"], js["nX"] * js["nY"], js["dtStep"]
    v = np.memmap(pdir / "v_full.bin", dtype="<f8", mode="r").reshape(nzpot * nr, nxy, -1)
    ntp = v.shape[2]
    nphi = 60
    ntheta = ntp // nphi
    fbar = np.zeros((nzpot * nr, ntheta), dtype=complex)
    vmean = np.zeros(nzpot * nr)
    for s in range(0, nzpot * nr, 250):
        blk = np.asarray(v[s:s + 250]).reshape(-1, nxy, ntheta, nphi)
        fbar[s:s + 250] = np.exp(-1j * blk * dt).mean(axis=(1, 3))
        vmean[s:s + 250] = blk.mean(axis=(1, 2, 3))
    vasy = np.array([float(l.split()[2]) for l in (pdir / "asymptotic_potential.dat").read_text().splitlines()
                     if l.strip() and not l.startswith("#")])
    sel = vmean < vasy.min() + 1.0 / CHAEV
    return fbar, sel, {"n_points": int(fbar.shape[0]), "n_selected": int(sel.sum()), "dt": dt}


# ----------------------------------------------------------------------------- main
def main():
    codes = {k: load_code(p) for k, p in CODES.items()}
    J = codes["fortran"]["J"]
    res = {"validity": validity(codes), "J": J}

    # grids and matrices
    res["thetas"] = {pair_key(a, b): ulp_compare(codes[a]["thetas"], codes[b]["thetas"]) for a, b in PAIRS[:3]}
    res["weights"] = {pair_key(a, b): ulp_compare(codes[a]["weights"], codes[b]["weights"]) for a, b in PAIRS[:3]}
    res["matrix_all"] = {pair_key(a, b): ulp_compare(codes[a]["matrix"], codes[b]["matrix"]) for a, b in PAIRS[:3]}
    res["matrix_mJ0"] = {pair_key(a, b): ulp_compare(codes[a]["m0_block"], codes[b]["m0_block"]) for a, b in PAIRS[:3]}

    # round trips: own arithmetic and matrix values alone
    R, Rmat, prods = {}, {}, {}
    res["round_trip"] = {}
    for k, c in codes.items():
        leak, other = round_trip_leak(c)
        prods[k] = exact_products(c)
        R[k], Rmat[k] = leak, matrix_leak(c, prods[k])
        res["round_trip"][k] = {"R_J": {int(j): float(leak[i].real) for i, j in enumerate(J)},
                                "R_J_imag_max": float(np.abs(leak.imag).max()),
                                "R_J_matrix_only": {int(j): float(Rmat[k][i]) for i, j in enumerate(J)},
                                "leak_outside_mJ0_max": other, "all_mJ0_functions": full_round_trip(c)}
    dR = {pair_key(a, b): R[a] - R[b] for a, b in PAIRS}
    dRmat = {pair_key(a, b): Rmat[a] - Rmat[b] for a, b in PAIRS}
    res["delta_R"] = {p: {int(j): float(v[i].real) for i, j in enumerate(J)} for p, v in dR.items()}
    res["delta_R_imag_max"] = {p: float(np.abs(v.imag).max()) for p, v in dR.items()}
    res["delta_R_matrix_only"] = {p: {int(j): float(v[i]) for i, j in enumerate(J)} for p, v in dRmat.items()}

    # observed family part of the BK differences
    lab = NativeBk(BK_FILES["fast"]).labels      # the analysed states: J = 2 ... JMax of the analysis
    fam = (lab["mJ"] == 0) & (lab["n"] == 0) & (lab["m"] == 0)
    Jf = [int(j) for j in J if j > 0 and j in set(lab["J"][fam].tolist())]
    fb = family_bk(Jf)
    w = fb["window"]
    norms = {p: np.sqrt(np.sum(np.abs(v[w]) ** 2, axis=0)) for p, v in fb["d"].items()}
    cnorm = np.sqrt(np.sum(np.abs(fb["cF"][w]) ** 2, axis=0))
    res["bk_family"] = {"window": {"n": int(w.sum()), "t_first": float(fb["times"][w][0]), "t_last": float(fb["times"][w][-1])},
                        "norm_c_fortran": {j: float(cnorm[i]) for i, j in enumerate(Jf)},
                        "norm_d": {p: {j: float(v[i]) for i, j in enumerate(Jf)} for p, v in norms.items()},
                        "norm_d_rel": {p: {j: float(v[i] / cnorm[i]) for i, j in enumerate(Jf)} for p, v in norms.items()},
                        "abs_d_max": {p: {j: float(np.abs(v[w, i]).max()) for i, j in enumerate(Jf)} for p, v in fb["d"].items()}}

    # predictions: T_J calibrated on strict - fast, prediction of strict - Fortran
    pred = {}
    for i, j in enumerate(Jf):
        k = list(J).index(j)
        r_sF, r_sf, r_fF = dR["strict-fortran"][k], dR["strict-fast"][k], dR["fast-fortran"][k]
        dsf, dsF, dfF = (fb["d"][p][w, i] for p in ("strict-fast", "strict-fortran", "fast-fortran"))
        lam = np.vdot(dsf, dsF) / np.vdot(dsf, dsf)
        coll = abs(np.vdot(dsf, dsF)) ** 2 / (np.vdot(dsf, dsf).real * np.vdot(dsF, dsF).real)
        ratio_R = (r_sF / r_sf) if r_sf != 0 else float("nan")
        predicted = norms["strict-fast"][i] * abs(r_sF) / abs(r_sf) if r_sf != 0 else float("nan")
        pred[j] = {"dR_strict_fortran": float(r_sF.real), "dR_strict_fast": float(r_sf.real), "dR_fast_fortran": float(r_fF.real),
                   "ratio_dR_strict_fortran_over_strict_fast": float(np.real(ratio_R)),
                   "observed_norm_strict_fortran": float(norms["strict-fortran"][i]),
                   "predicted_norm_strict_fortran": float(predicted),
                   "predicted_over_observed": float(predicted / norms["strict-fortran"][i]),
                   "transfer_T": float(norms["strict-fast"][i] / abs(r_sf)) if r_sf != 0 else float("nan"),
                   "collinearity": float(coll), "lambda_fit": [float(lam.real), float(lam.imag)],
                   "sign_matches": bool(np.sign(lam.real) == np.sign(np.real(ratio_R))),
                   "ratio_observed_fast_over_strict": float(norms["fast-fortran"][i] / norms["strict-fortran"][i]),
                   "ratio_predicted_fast_over_strict": float(abs(r_fF) / abs(r_sF)) if r_sF != 0 else float("nan"),
                   "collinearity_fast_fortran_vs_strict_fast": float(abs(np.vdot(dsf, dfF)) ** 2 /
                                                                     (np.vdot(dsf, dsf).real * np.vdot(dfF, dfF).real))}
    res["prediction"] = pred

    # verdict (criteria.md, 2026-11-08)
    def explains(j):
        p = pred[j]
        return 1 / 1.5 <= p["predicted_over_observed"] <= 1.5 and p["collinearity"] >= 0.9 and p["sign_matches"]
    within3 = all(1 / 3 <= pred[j]["predicted_over_observed"] <= 3 for j in VERDICT_J)
    zero = any(pred[j]["dR_strict_fortran"] == 0 for j in VERDICT_J)
    if not zero and all(explains(j) for j in VERDICT_J):
        verdict = "explains"
    elif not zero and (within3 or any(explains(j) for j in VERDICT_J)):
        verdict = "partly"
    else:
        verdict = "does not"
    res["criteria"] = {"factor_explains": 1.5, "factor_partly": 3.0, "collinearity_min": 0.9, "verdict_J": list(VERDICT_J)}
    res["summary"] = {
        "rounding_part_max_J_gt_0": float(max(np.abs(R[k][1:].real - Rmat[k][1:]).max() for k in codes)),
        "own_leak_max_J_gt_0": {k: float(np.abs(R[k][1:]).max()) for k in codes},
        "non_collinear_fraction": {j: float(1 - pred[j]["collinearity"]) for j in VERDICT_J}}
    res["verdict"] = {"value": verdict, "explains_by_J": {j: explains(j) for j in VERDICT_J}, "within_factor_3": within3,
                      "dR_strict_fortran_zero": zero}

    # robustness: the leak difference with the averaged exp(-iV dt) (matrix values, exact products)
    fbar, sel, info = potential_weights()
    rob = {"points": info}
    for j in VERDICT_J:
        k = list(J).index(j)
        nphi = codes["fortran"]["nphi"]
        def dG(a, b):
            delta = np.array([float(prods[a][t][k] - prods[b][t][k]) for t in range(len(prods[a]))])
            return nphi * fbar[sel] @ delta
        ratio = dG("strict", "fortran") / dG("strict", "fast")
        ref = dRmat["strict-fortran"][k] / dRmat["strict-fast"][k]
        a = np.abs(ratio)
        rob[j] = {"ratio_abs_min": float(a.min()), "ratio_abs_median": float(np.median(a)), "ratio_abs_max": float(a.max()),
                  "ratio_without_potential": float(ref),
                  "median_over_without": float(np.median(a) / abs(ref))}
    res["robustness_potential"] = rob

    # post hoc (criteria.md, 2026-11-08, after the results): the same ratio at every (Z, r) point of
    # the interaction region, by the anisotropy of the averaged exp(-iV dt); and the fast helper's
    # theta nodes (1 ulp from the run's at a few nodes) moved to the run's
    aniso = np.abs(fbar - fbar.mean(axis=1, keepdims=True)).max(axis=1)
    zr = {"Z": np.arange(fbar.shape[0]) // json.loads((INPUT_DIR / "input.json").read_text())["nR"]}
    dz = json.loads((INPUT_DIR / "input.json").read_text())["dZ"]
    post = {"all_points": {}}
    for j in VERDICT_J:
        k = list(J).index(j)
        nphi = codes["fortran"]["nphi"]
        delta = {p: np.array([float(prods[a][t][k] - prods[b][t][k]) for t in range(len(prods[a]))])
                 for p, (a, b) in {"strict-fortran": ("strict", "fortran"), "strict-fast": ("strict", "fast")}.items()}
        a = np.abs((fbar @ delta["strict-fortran"]) / (fbar @ delta["strict-fast"]))
        classes = {}
        for lo, hi, name in ((0, 1e-2, "anisotropy<1e-2"), (1e-2, 1e-1, "1e-2..1e-1"), (1e-1, 10, "anisotropy>1e-1")):
            m = (aniso >= lo) & (aniso < hi)
            classes[name] = {"n": int(m.sum()), "ratio_abs_median": float(np.median(a[m])),
                             "Z_max": float(zr["Z"][m].max() * dz)}
        post["all_points"][j] = {"ratio_abs_percentiles_1_50_90_99": [float(x) for x in np.percentile(a, [1, 50, 90, 99])],
                                 "by_anisotropy": classes,
                                 "nodes_strict_fortran": [int(t) for t in np.flatnonzero(delta["strict-fortran"])],
                                 "node0_share_strict_fast": float(delta["strict-fast"][0] / delta["strict-fast"].sum())}
    sens, n_nodes = fast_theta_sensitivity(codes, list(VERDICT_J))
    post["fast_helper_thetas"] = {"n_nodes_differing_from_run": n_nodes}
    for i, j in enumerate(VERDICT_J):
        k = list(J).index(j)
        r_sf_run = dR["strict-fast"][k].real - sens[i]       # strict - fast with the run's nodes
        p = pred[j]
        p_run = norms["strict-fast"][Jf.index(j)] * abs(p["dR_strict_fortran"]) / abs(r_sf_run)
        post["fast_helper_thetas"][j] = {"dR_change": float(-sens[i]), "dR_strict_fast_run_nodes": float(r_sf_run),
                                         "relative_change": float(-sens[i] / dR["strict-fast"][k].real),
                                         "predicted_over_observed_run_nodes": float(p_run / p["observed_norm_strict_fortran"])}
    res["post_hoc"] = post

    save_item("angular_matrices", res)
    print(json.dumps({"validity": res["validity"], "verdict": res["verdict"],
                      "prediction": {j: {k_: v for k_, v in pred[j].items() if k_ in ("dR_strict_fortran", "dR_strict_fast",
                       "observed_norm_strict_fortran", "predicted_norm_strict_fortran", "predicted_over_observed", "collinearity",
                       "ratio_observed_fast_over_strict", "ratio_predicted_fast_over_strict")} for j in Jf}}, indent=1))


if __name__ == "__main__":
    main()
