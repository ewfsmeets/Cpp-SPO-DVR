#!/usr/bin/env python3
"""Item 18b: the replay. The helper `replay` (C++ on-the-fly analysis, SMatrixAccumulator, linked
against the build that ran) re-accumulates the S-matrix from each run's BK file.

  self-check        run's own tables vs replay(run's own BK file)
  propagation-only  replay(fort.21) vs replay(C++ BK file), same helper binary: the only
                    difference is the BK coefficients, i.e. the propagation (and BK projection)
                    of the two codes, at full precision
  analysis-only     Fortran printed values vs replay(fort.21): the same BK coefficients through
                    the Fortran and the C++ analysis code, at the Fortran's printed precision
Writes data/items/item18b.json and data/replay_reaction_probability.tsv."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from common import (CHAEV, DATA, RUN_DIRS, WORK, FortranBk, NativeBk, clean_log_lines, parse_fortran_real,
                    printed_half_unit, save_item, stats)
from parse_logs import energy_blocks, fort25

REPLAY = WORK / "replay"
KEYS = ["E_coll", "v", "j", "mJ", "n", "m"]


def load(path, name):
    return pd.read_csv(path / name, sep="\t", comment="#")


def s_diff(a, b):
    assert (a[KEYS].to_numpy() == b[KEYS].to_numpy()).all()
    dS = np.abs((a["Re_S"] - b["Re_S"]).to_numpy() + 1j * (a["Im_S"] - b["Im_S"]).to_numpy())
    dP = np.abs((a["P"] - b["P"]).to_numpy())
    where = lambda i: {k: (float(b[k].iloc[i]) if k == "E_coll" else int(b[k].iloc[i])) for k in KEYS}
    sel = b["P"].to_numpy() > 1e-8
    dph = np.abs(np.angle(np.exp(1j * (a["phase"].to_numpy()[sel] - b["phase"].to_numpy()[sel]))))
    return {"S_abs": stats(dS, where), "P_abs": stats(dP, where), "phase_rad_P_gt_1e-8": stats(dph),
            "n_bitwise_equal_rows": int(np.sum(dS == 0) if False else np.sum((a["Re_S"].to_numpy() == b["Re_S"].to_numpy())
                                                                           & (a["Im_S"].to_numpy() == b["Im_S"].to_numpy())))}


def summary_diff(a, b):
    out = {}
    for col in ("P_scattered", "P_reaction", "P_v0", "P_v1", "E_vib_mean", "E_rot_mean", "E_dif_mean"):
        d = np.abs(a[col].to_numpy() - b[col].to_numpy())
        out[col] = stats(d, where=[float(x) for x in b["E_coll"]])
    return out


def rovib_diff(a, b):
    d = np.abs(a["P"].to_numpy() - b["P"].to_numpy())
    return stats(d, where=lambda i: {"E": float(b["E_coll"].iloc[i]), "v": int(b["v"].iloc[i]), "j": int(b["j"].iloc[i]),
                                     "mJ": int(b["mJ"].iloc[i])})


def files_identical(a, b):
    return all((a / f).read_bytes() == (b / f).read_bytes() for f in
               ("s_matrix.tsv", "scattering_summary.tsv", "probabilities_rovib.tsv", "probabilities_diffraction.tsv"))


def analysis_rounding_bound(M):
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    nb = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    lab = nb.labels
    index = {(int(lab["J"][i]), int(lab["mJ"][i]), int(lab["n"][i]), int(lab["m"][i])): i for i in range(nb.nstata)}
    is_ = np.array([index[(int(j), int(mj), int(n), int(m))] for j, mj, n, m in zip(M["j"], M["mJ"], M["n"], M["m"])])
    v = M["v"].to_numpy()
    e_init_cpp = float([l.split() for l in (RUN_DIRS["cpu"] / "gas_phase_states.dat").read_text().splitlines()
                        if not l.startswith("#")][0][2])
    dE = abs(fb.e_initial - e_init_cpp)
    t = fb.times()
    E_all = np.unique(M["E_total"].to_numpy()) / CHAEV                      # total energies (hartree), ascending
    # per (energy, t): dE t + ulp(|E| t) + 4 eps  (systematic drift of the phase reference, rounding
    # of the phase argument, rounding of the product and sum)
    w = dE * t[None, :] + np.spacing(np.abs(E_all[:, None] * t[None, :])) + 4 * np.finfo(float).eps
    acc = np.zeros((len(E_all), fb.nviba * fb.nstata))
    for k in range(fb.n_records):
        acc += w[:, k:k + 1] * np.abs(fb.record(k)).reshape(1, -1)
    ie = np.searchsorted(E_all, M["E_total"].to_numpy() / CHAEV)
    mass = 3.6742170840000003e+03
    e_tot = M["E_total"].to_numpy() / CHAEV
    kz = np.sqrt(2 * mass * (e_tot - nb.energies[v, is_]))
    k0 = np.sqrt(2 * mass * M["E_coll"].to_numpy() / CHAEV)
    from parse_logs import cpp_dotted, cpp_section
    cin = cpp_dotted(cpp_section(clean_log_lines("cpu"), "INPUT"))
    kz0, alpha = float(cin["Initial kz"].split()[0]), float(cin["Wave packet alpha"].split()[0])
    weight = np.exp(-((kz0 - k0) ** 2) / (2 * alpha))
    factor = (0.5 * alpha / np.pi) ** 0.25 / mass * np.sqrt(kz * k0 / weight)
    # set-up factor: kz = sqrt(2M(E - E_ch)) enters through exp(-i kz Zana) and sqrt(kz); the
    # Fortran's E - E_ch differs by d_rel = (eviroi_F - eviroi_C) - (enymo - E_ch,C)
    d_rel = np.abs((fb.e_initial - e_init_cpp) - (fb.enymo.reshape(fb.nviba, fb.nstata)[v, is_] - nb.energies[v, is_]))
    dkz = mass * d_rel / kz
    zana = 15.0
    S = np.abs(M["Re_S"].to_numpy() + 1j * M["Im_S"].to_numpy())
    factor_term = S * (zana * dkz + dkz / (2 * kz))
    analysis_rounding_bound.info = {"dE_eviroi_hartree": dE, "T": float(t[-1]),
                                    "phase_argument_max_rad": float(np.max(np.abs(E_all)) * t[-1]),
                                    "phase_argument_max_ulp_rad": float(np.spacing(np.max(np.abs(E_all)) * t[-1])),
                                    "bound": "|factor| dt sum_t |c_t| (dE t + ulp(|E| t) + 4 eps) + |S| (Zana dkz + dkz/(2 kz)),"
                                             " dkz = M |d(E - E_ch)| / kz"}
    return factor * fb.deltat * acc[ie, v * fb.nstata + is_] + factor_term


def analysis_setup_correction(M):
    """Replay S (C++ analysis of the Fortran BK file) with the C++ set-up scalars replaced by the
    Fortran's: channel energies (enymo), initial k_Z and alpha (fort.21 records 2, 10)."""
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    nb = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    lab = nb.labels
    index = {(int(lab["J"][i]), int(lab["mJ"][i]), int(lab["n"][i]), int(lab["m"][i])): i for i in range(nb.nstata)}
    is_ = np.array([index[(int(j), int(mj), int(n), int(m))] for j, mj, n, m in zip(M["j"], M["mJ"], M["n"], M["m"])])
    v = M["v"].to_numpy()
    e_c = nb.energies[v, is_]
    e_f = fb.enymo.reshape(nb.nviba, nb.nstata)[v, is_]
    mass = 3.6742170840000003e+03                     # total mass, C++ log (identical in both codes' input)
    e_tot = M["E_total"].to_numpy() / CHAEV
    kz_c = np.sqrt(2 * mass * (e_tot - e_c))
    kz_f = np.sqrt(2 * mass * (e_tot - e_f))
    k0 = np.sqrt(2 * mass * M["E_coll"].to_numpy() / CHAEV)
    import json
    from common import INPUT_DIR
    zana = json.loads((INPUT_DIR / "input.json").read_text())["Zana"]
    z0 = json.loads((INPUT_DIR / "input.json").read_text())["Z0"]
    from parse_logs import cpp_dotted, cpp_section
    cin = cpp_dotted(cpp_section(clean_log_lines("cpu"), "INPUT"))
    kz0_c, a_c = float(cin["Initial kz"].split()[0]), float(cin["Wave packet alpha"].split()[0])
    kz0_f, a_f = fb.kz, fb.alpha

    def factor(kz, kz0, alpha):
        weight = np.exp(-((kz0 - k0) ** 2) / (2 * alpha))
        pref = (0.5 * alpha / np.pi) ** 0.25 / mass
        return pref * np.exp(1j * (-kz * zana + (kz0 - k0) * z0 - np.pi)) * np.sqrt(kz * k0 / weight)

    s = M["Re_S"].to_numpy() + 1j * M["Im_S"].to_numpy()
    return s * (factor(kz_f, kz0_f, a_f) / factor(kz_c, kz0_c, a_c))


def decompose(A, B_):
    """Decomposition of the difference of two replayed S-matrices per energy: S_a = S_b (1 + alpha + i beta)
    on average (weighted by |S_b|^2): alpha = relative amplitude, beta = phase; the predicted
    dSum P = 2 alpha Sum P against the measured one; and the per-channel residual."""
    sa = A["Re_S"].to_numpy() + 1j * A["Im_S"].to_numpy()
    sb = B_["Re_S"].to_numpy() + 1j * B_["Im_S"].to_numpy()
    rows = []
    for e in np.unique(B_["E_coll"].to_numpy()):
        m = B_["E_coll"].to_numpy() == e
        z = np.sum(np.conj(sb[m]) * (sa[m] - sb[m])) / np.sum(np.abs(sb[m]) ** 2)
        dsum = float(np.sum(A["P"].to_numpy()[m] - B_["P"].to_numpy()[m]))
        psum = float(np.sum(B_["P"].to_numpy()[m]))
        rows.append({"E": float(e), "alpha": float(z.real), "beta": float(z.imag), "dsumP": dsum,
                     "predicted_2alpha_sumP": 2 * float(z.real) * psum})
    # per-channel residual after removing the uniform amplitude-and-phase factor of its energy:
    # a logic error that touches only some channels would show here
    zE = {r["E"]: complex(r["alpha"], r["beta"]) for r in rows}
    zz = np.array([zE[float(e)] for e in B_["E_coll"].to_numpy()])
    resid = np.abs((sa - sb) - zz * sb)
    selr = np.abs(sb) ** 2 > 1e-8
    rel = resid[selr] / np.abs(sb[selr])
    kmax = int(np.argmax(rel))
    idx = np.flatnonzero(selr)[kmax]
    # where the non-uniform part sits: the specular rotationally inelastic family
    fam = ((B_["n"] == 0) & (B_["m"] == 0) & (B_["mJ"] == 0) & (B_["j"] > 0) & (B_["v"] == 0)).to_numpy()
    tot2 = float(np.sum(np.abs(sa - sb) ** 2))
    res2 = float(np.sum(resid ** 2))
    fam_j = {}
    for jj in sorted(set(B_["j"].to_numpy()[fam].tolist())):
        mj = fam & (B_["j"].to_numpy() == jj)
        fam_j[int(jj)] = {"rel_residual_normweighted": float(np.sqrt(np.sum(resid[mj] ** 2) / np.sum(np.abs(sb[mj]) ** 2))),
                          "rel_residual_max_where_P_gt_1e-8": float(np.max(np.where(np.abs(sb[mj]) ** 2 > 1e-8,
                                                                                    resid[mj] / np.abs(sb[mj]), 0)))}
    per_el = {}
    for jj in sorted(set(B_["j"].to_numpy()[fam].tolist())):
        mj = fam & (B_["j"].to_numpy() == jj) & (np.abs(sb) > 0)
        r_el = np.abs(sa[mj] - sb[mj]) / np.abs(sb[mj])
        per_el[int(jj)] = {"median": float(np.median(r_el)), "max": float(r_el.max())}
    bins = {}
    P_ = np.abs(sb) ** 2
    for lo, hi, name in ((1e-2, 1.0, "P>1e-2"), (1e-4, 1e-2, "1e-4<P<1e-2"), (1e-6, 1e-4, "1e-6<P<1e-4"), (1e-8, 1e-6, "1e-8<P<1e-6")):
        mb = (P_ > lo) & (P_ <= hi)
        rb = resid[mb] / np.abs(sb[mb])
        bins[name] = {"n": int(mb.sum()), "p95": float(np.percentile(rb, 95)), "max": float(rb.max()),
                      "n_above_1e-10": int(np.sum(rb > 1e-10))}
    family = {"per_element_relative": per_el, "residual_by_size": bins,
              "residual_norm2_fraction_of_total": float(res2 / tot2),
              "residual_norm_fraction_of_total_dS": float(np.sqrt(res2 / tot2)),
              "family_fraction_of_residual_norm2": float(np.sum(resid[fam] ** 2) / res2),
              "residual_rel_normweighted_outside_family": float(np.sqrt(np.sum(resid[~fam] ** 2) / np.sum(np.abs(sb[~fam]) ** 2))),
              "per_J": fam_j}
    al = np.array([r["alpha"] for r in rows]); be = np.array([r["beta"] for r in rows])
    kappa_corr_inv_sqrtE = float(np.corrcoef(al, 1 / np.sqrt(np.array([r["E"] for r in rows])))[0, 1])
    ratio = np.array([r["dsumP"] / r["predicted_2alpha_sumP"] for r in rows])
    big = np.abs(sb) ** 2 > 1e-2
    return {"rel_dS_large_channels_median": float(np.median(np.abs(sa[big] - sb[big]) / np.abs(sb[big]))),
           "rel_dS_large_channels_max": float(np.max(np.abs(sa[big] - sb[big]) / np.abs(sb[big]))),
           "n_large_channels": int(big.sum()), "kappa_corr_inv_sqrtE": kappa_corr_inv_sqrtE, "max_P_single_channel": float(np.max(np.abs(sb) ** 2)),
           "specular_inelastic_family": family, "residual_rel_median": float(np.median(rel)), "residual_rel_p95": float(np.percentile(rel, 95)),
           "residual_rel_max": float(rel.max()), "residual_abs_max": float(resid.max()),
           "residual_n": int(selr.sum()),
           "residual_argmax": {k_: (float(B_[k_].iloc[idx]) if k_ == "E_coll" else int(B_[k_].iloc[idx]))
                               for k_ in KEYS} | {"abs_S": float(abs(sb[idx]))},
           "alpha_min": float(al.min()), "alpha_max": float(al.max()), "alpha_median": float(np.median(al)),
           "beta_abs_max": float(np.abs(be).max()), "beta_median": float(np.median(be)),
           "dsumP_over_2alpha_sumP_min": float(ratio.min()), "dsumP_over_2alpha_sumP_max": float(ratio.max()),
           "by_energy": rows}


def main():
    res = {"helper_builds_identical": {
        "fortran": files_identical(REPLAY / "build-fast" / "fortran", REPLAY / "build-cuda-fast" / "fortran"),
        "cpu": files_identical(REPLAY / "build-fast" / "cpu", REPLAY / "build-cuda-fast" / "cpu")}}
    rp = {k: REPLAY / "build-fast" / k for k in ("cpu", "gpu", "fortran")}
    sm = {f"replay_{k}": load(p, "s_matrix.tsv") for k, p in rp.items()}
    sm.update({k: load(RUN_DIRS[k], "s_matrix.tsv") for k in ("cpu", "gpu")})
    su = {f"replay_{k}": load(p, "scattering_summary.tsv") for k, p in rp.items()}
    su.update({k: load(RUN_DIRS[k], "scattering_summary.tsv") for k in ("cpu", "gpu")})
    rv = {f"replay_{k}": load(p, "probabilities_rovib.tsv") for k, p in rp.items()}
    rv.update({k: load(RUN_DIRS[k], "probabilities_rovib.tsv") for k in ("cpu", "gpu")})

    # self-check (run vs replay of its own file)
    sc = {}
    for k in ("cpu", "gpu"):
        sc[k] = {"bitwise": files_identical(RUN_DIRS[k], rp[k]), "s_matrix": s_diff(sm[k], sm[f"replay_{k}"]),
                 "summary": summary_diff(su[k], su[f"replay_{k}"]), "rovib": rovib_diff(rv[k], rv[f"replay_{k}"])}
    oc = (sm["cpu"]["Re_S"] - sm["replay_cpu"]["Re_S"]).to_numpy() + 1j * (sm["cpu"]["Im_S"] - sm["replay_cpu"]["Im_S"]).to_numpy()
    og = (sm["gpu"]["Re_S"] - sm["replay_gpu"]["Re_S"]).to_numpy() + 1j * (sm["gpu"]["Im_S"] - sm["replay_gpu"]["Im_S"]).to_numpy()
    sc["offset_same_in_both_runs_max_abs_diff"] = float(np.abs(oc - og).max())
    res["self_check"] = sc

    # propagation-only: replay(fort.21) vs replay(C++), and C++ vs C++
    po = {}
    for a, b in (("replay_cpu", "replay_fortran"), ("replay_gpu", "replay_fortran"), ("replay_gpu", "replay_cpu")):
        po[f"{a[7:]}-{b[7:]}"] = {"s_matrix": s_diff(sm[a], sm[b]), "summary": summary_diff(su[a], su[b]),
                                  "rovib": rovib_diff(rv[a], rv[b])}
    res["propagation_only"] = po
    # decomposition of the propagation-only difference (see decompose)
    dec = {f"{a}-{bref}": decompose(sm[f"replay_{a}"], sm[f"replay_{bref}"])
           for a, bref in (("cpu", "fortran"), ("gpu", "fortran"), ("gpu", "cpu"))}
    res["amplitude_decomposition"] = dec
    # overall norm-weighted relative difference of the S-matrices, every pair
    rel_all = {}
    for a, b in (("replay_cpu", "replay_fortran"), ("replay_gpu", "replay_fortran"), ("replay_gpu", "replay_cpu")):
        sa = sm[a]["Re_S"].to_numpy() + 1j * sm[a]["Im_S"].to_numpy()
        sb = sm[b]["Re_S"].to_numpy() + 1j * sm[b]["Im_S"].to_numpy()
        rel_all[f"{a[7:]}-{b[7:]}"] = float(np.sqrt(np.sum(np.abs(sa - sb) ** 2) / np.sum(np.abs(sb) ** 2)))
    res["S_relative_normweighted"] = rel_all
    # run tables vs replay(fort.21) (includes the run-vs-helper offset), used as dS for printed comparisons
    res["replay_vs_runs"] = {f"{k}-fortran": {"S_max_abs": s_diff(sm[k], sm["replay_fortran"])["S_abs"]["max"],
                                              "P_reaction_max_abs": summary_diff(su[k], su["replay_fortran"])["P_reaction"]["max"]}
                             for k in ("cpu", "gpu")}

    # analysis-only: Fortran printed (stdout, fort.25) vs replay(fort.21)
    blocks = energy_blocks(clean_log_lines("fortran"))
    f_su = su["replay_fortran"]
    f_rv = rv["replay_fortran"].set_index(["E_coll", "v", "j", "mJ"])["P"]
    f_df = load(rp["fortran"], "probabilities_diffraction.tsv").set_index(["E_coll", "v", "j", "mJ", "order"])["P"]
    E = f_su["E_coll"].to_numpy()
    ao = {}
    def tally(name, pairs, floor=0.0):
        ratios = []
        worst = {"ratio": -1}
        for c, text, where in pairs:
            tol = max(printed_half_unit(text), floor)
            r = abs(c - parse_fortran_real(text)) / tol
            ratios.append(r)
            if r > worst["ratio"]:
                worst = {"ratio": r, "where": where, "fortran": text, "replay": c}
        ratios = np.array(ratios)
        ao[name] = {"n": int(ratios.size), "n_within": int((ratios <= 1).sum()), "worst": worst,
                    "n_above_half_unit": int((ratios > 1).sum())}
    tally("P_reaction (half unit, no floor)", [(f_su["P_reaction"][i], b["p_reaction"], float(E[i])) for i, b in enumerate(blocks)])
    tally("P_sum", [(f_su["P_scattered"][i], b["p_sum"], float(E[i])) for i, b in enumerate(blocks)])
    tally("P_v", [(f_su[f"P_v{v}"][i], t, (float(E[i]), v)) for i, b in enumerate(blocks) for v, t in b["p_v"].items()])
    tally("P_vjm", [(float(f_rv.loc[(E[i],) + key]), t, (float(E[i]),) + key) for i, b in enumerate(blocks) for key, t in b["p_vjm"].items()])
    tally("P_order", [(float(f_df.loc[(E[i],) + key]), t, (float(E[i]),) + key) for i, b in enumerate(blocks) for key, t in b["p_order"].items()])
    tally("E_avg", [(f_su[col][i], b[key], (float(E[i]), col)) for i, b in enumerate(blocks)
                    for col, key in (("E_vib_mean", "e_vib"), ("E_rot_mean", "e_rot"), ("E_dif_mean", "e_dif"))])
    # fort.25 S-matrix vs replay(fort.21)
    F = pd.DataFrame(list(fort25(RUN_DIRS["fortran"] / "fort.25")),
                     columns=["ie", "v", "j", "mJ", "order", "n", "m", "ReS", "ImS", "P", "phase"])
    F["E_coll"] = E[F["ie"].to_numpy()]
    M = F.merge(sm["replay_fortran"], on=KEYS, how="inner", suffixes=("", "_r"))
    for col, fcol, rcol in (("fort25 Re S", "ReS", "Re_S"), ("fort25 Im S", "ImS", "Im_S"), ("fort25 |S|^2", "P", "P_r")):
        ft = M[fcol].to_numpy()
        f = np.array([parse_fortran_real(x) for x in ft])
        h = np.array([printed_half_unit(x) for x in ft])
        r = np.abs(M[rcol].to_numpy() - f) / h
        i = int(np.argmax(r))
        ao[col] = {"n": int(r.size), "n_within": int((r <= 1).sum()), "n_above_half_unit": int((r > 1).sum()),
                   "worst": {"ratio": float(r[i]), "fortran": ft[i], "replay": float(M[rcol].iloc[i]),
                             "where": {k: float(M[k].iloc[i]) for k in KEYS}}}
    # Re/Im S above half a printed unit: correct the replay for the Fortran's own set-up scalars
    # (channel energies enymo, initial k_Z and alpha from fort.21; s_matrix.cpp: S = sum * prefactor
    # * exp(i arg) * sqrt(kz k0 / weight), arg = -kz Zana + (kz0 - k0) Z0 - pi) and test the print.
    corr = analysis_setup_correction(M)
    for col, fcol, rcol in (("fort25 Re S (replay corrected for Fortran set-up)", "ReS", "Re_S"),
                            ("fort25 Im S (replay corrected for Fortran set-up)", "ImS", "Im_S")):
        ft = M[fcol].to_numpy()
        f = np.array([parse_fortran_real(x) for x in ft])
        h = np.array([printed_half_unit(x) for x in ft])
        raw_bad = np.abs(M[rcol].to_numpy() - f) > h
        val = corr.real if "Re" in col else corr.imag
        r = np.abs(val - f) / h
        ao[col] = {"n_raw_above_half_unit": int(raw_bad.sum()), "n_corrected_above_half_unit": int((r > 1).sum()),
                   "n_raw_above_now_within": int((raw_bad & (r <= 1)).sum()),
                   "worst_corrected_ratio_of_raw_outliers": float(r[raw_bad].max()) if raw_bad.any() else 0.0}
    # Bound for analysis-side differences (Fortran analysis vs C++ analysis of the same BK file):
    # (1) the analysis energies' reference eviroi differs (LAPACK vs Jacobi, item 6) by dE, so the
    #     phase factors exp(iEt) drift apart by dE t: |dS| <= |S| |dE| T;
    # (2) the absolute rounding floor of the 875-term sum: 4 eps |factor| dt sum_t |c_t|.
    bound = analysis_rounding_bound(M)
    res_rows = []
    for col, fcol, rcol in (("Re S", "ReS", "Re_S"), ("Im S", "ImS", "Im_S")):
        ft = M[fcol].to_numpy()
        f = np.array([parse_fortran_real(x) for x in ft])
        h = np.array([printed_half_unit(x) for x in ft])
        d = np.abs(M[rcol].to_numpy() - f)
        bad = d > h
        ok = d <= h + bound
        ao[f"fort25 {col}: above half unit, within half unit + bound"] = {
            "n_above_half_unit": int(bad.sum()), "n_explained": int((bad & ok).sum()), "n_unexplained": int((bad & ~ok).sum()),
            "max_excess_over_bound": float(((d - h) / bound)[bad].max()) if bad.any() else 0.0}
        for i in np.flatnonzero(bad):
            res_rows.append({"component": col, **{k: (float(M[k].iloc[i]) if k == "E_coll" else int(M[k].iloc[i])) for k in KEYS},
                             "fortran": ft[i], "replay": float(M[rcol].iloc[i]), "abs_diff": float(d[i]),
                             "half_unit": float(h[i]), "bound": float(bound[i]), "explained": bool(ok[i])})
    ao["analysis_side_S_outliers"] = res_rows
    sel = np.array([parse_fortran_real(x) for x in M["P"]]) > 1e-8
    dph = np.abs(np.angle(np.exp(1j * (M["phase_r"].to_numpy()[sel] - np.array([parse_fortran_real(x) for x in M["phase"][sel]])))))
    ao["fort25 phase (|S|^2 > 1e-8)"] = {"n": int(sel.sum()), "n_within": int((dph <= 5.05e-7).sum()), "max_abs_rad": float(dph.max())}
    res["analysis_only"] = ao
    res["analysis_bound_info"] = getattr(analysis_rounding_bound, "info", {})

    # global bias of the analysis: signed residuals replay(fort.21) - fort.25 printed values, in
    # half units; pure print rounding gives a uniform distribution (std 1/sqrt(3) = 0.577 half units)
    Pf = np.array([parse_fortran_real(x) for x in M["P"]])
    hP = np.array([printed_half_unit(x) for x in M["P"]]) / 1.01
    sel = Pf > 1e-6
    rr = (M["P_r"].to_numpy()[sel] - Pf[sel])
    sig = hP[sel] / np.sqrt(3.0)                       # rounding std of each printed value
    relw = 1.0 / (sig / Pf[sel]) ** 2
    rel_bias = float(np.sum(relw * rr / Pf[sel]) / np.sum(relw))
    rel_err = float(1.0 / np.sqrt(np.sum(relw)))
    # symmetry-equivalent channels print identical values: count each (energy, printed value) once
    keyv = (M["ie"].astype(str) + "|" + M["P"]).to_numpy()[sel]
    _, first = np.unique(keyv, return_index=True)
    rw, rr_, pf_ = relw[first], rr[first], Pf[sel][first]
    dedup_bias = float(np.sum(rw * rr_ / pf_) / np.sum(rw))
    dedup_err = float(1.0 / np.sqrt(np.sum(rw)))
    ph_f = np.array([parse_fortran_real(x) for x in M["phase"]])
    selp = (Pf > 1e-6)
    dph = np.angle(np.exp(1j * (M["phase_r"].to_numpy()[selp] - ph_f[selp])))
    sig_ph = 0.5e-6 / np.sqrt(3.0)
    ao["global_bias"] = {"n_P": int(sel.sum()), "P_relative_bias": rel_bias, "P_relative_bias_stderr": rel_err,
                         "residual_std_in_half_units": float(np.std(rr / hP[sel])),
                         "P_relative_bias_dedup": dedup_bias, "P_relative_bias_dedup_stderr": dedup_err, "n_P_dedup": int(first.size), "expected_std_pure_rounding": float(1 / np.sqrt(3)),
                         "n_phase": int(selp.sum()), "phase_bias_rad": float(np.mean(dph)),
                         "phase_bias_stderr_rad": float(sig_ph / np.sqrt(selp.sum())),
                         "phase_residual_std_in_half_units": float(np.std(dph / 0.5e-6))}
    # per-energy P_reaction table (ladder, Figure 1 bottom panel)
    lines = ["# P_reaction at full precision through the same C++ analysis (helper `replay`, build-fast): replay of "
             "each run's BK file; and the Fortran printed value.",
             "E_coll\tP_replay_fortran\tP_replay_cpu\tP_replay_gpu\tabs_cpu_minus_F\tabs_gpu_minus_F\tabs_gpu_minus_cpu\tP_fortran_text"
             "\tmax_dPvjm_cpu_F\tmax_dPvjm_gpu_F"]
    rc, rg, rf = su["replay_cpu"]["P_reaction"], su["replay_gpu"]["P_reaction"], su["replay_fortran"]["P_reaction"]
    vjm = {}
    for k in ("cpu", "gpu"):
        d = np.abs(rv[f"replay_{k}"]["P"].to_numpy() - rv["replay_fortran"]["P"].to_numpy())
        vjm[k] = pd.Series(d).groupby(rv["replay_fortran"]["E_coll"].to_numpy()).max().to_numpy()
    for i, b in enumerate(blocks):
        lines.append(f"{E[i]:.4f}\t{rf[i]:.16e}\t{rc[i]:.16e}\t{rg[i]:.16e}\t{abs(rc[i]-rf[i]):.3e}\t{abs(rg[i]-rf[i]):.3e}"
                     f"\t{abs(rg[i]-rc[i]):.3e}\t{b['p_reaction']}\t{vjm['cpu'][i]:.3e}\t{vjm['gpu'][i]:.3e}")
    (DATA / "replay_reaction_probability.tsv").write_text("\n".join(lines) + "\n")
    save_item("item18b", res)

    print("helper builds identical:", res["helper_builds_identical"])
    for k in ("cpu", "gpu"):
        s = sc[k]
        print(f"self-check {k}: bitwise {s['bitwise']}; |dS| max {s['s_matrix']['S_abs']['max']:.2e}, |dP| max "
              f"{s['s_matrix']['P_abs']['max']:.2e}, |dP_reaction| max {s['summary']['P_reaction']['max']:.2e}")
    print("offset identical in both runs to", sc["offset_same_in_both_runs_max_abs_diff"])
    for k, d in res["amplitude_decomposition"].items():
        print("amplitude", k, {kk: v for kk, v in d.items() if kk != "by_energy"})
    for pair, r in po.items():
        print(f"propagation-only {pair}: |dS| max {r['s_matrix']['S_abs']['max']:.2e} at {r['s_matrix']['S_abs']['argmax']}; "
              f"|dP| max {r['s_matrix']['P_abs']['max']:.2e}; P_reaction max {r['summary']['P_reaction']['max']:.2e} "
              f"(median {r['summary']['P_reaction']['median']:.2e}) at E={r['summary']['P_reaction']['argmax']}; "
              f"P(v,j,mJ) max {r['rovib']['max']:.2e}")
    for name, r in ao.items():
        if name == "analysis_side_S_outliers":
            continue
        print("analysis-only", name, {k: v for k, v in r.items() if k != "worst"}, "worst ratio", round(r.get("worst", {}).get("ratio", 0), 3))
    print("bound info", res["analysis_bound_info"])


if __name__ == "__main__":
    main()
