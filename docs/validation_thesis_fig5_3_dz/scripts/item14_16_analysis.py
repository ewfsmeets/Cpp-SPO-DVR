#!/usr/bin/env python3
"""Items 14-16: the on-the-fly analysis results of the three runs.

14  P_reaction, sum of P, P(v), P(v, j, mJ), average final energies (Fortran stdout per-energy
    blocks, printed) vs scattering_summary.tsv / probabilities_rovib.tsv (full precision)
15  P per diffraction order (stdout) vs probabilities_diffraction.tsv
16  S-matrix: fort.25 (Re S, Im S, |S|^2 as D13.6, phase F10.6) vs s_matrix.tsv, joined on
    (E, v, j, mJ, n, m); row-set differences reported
CPU <-> GPU: all tables at full precision.

Printed comparisons: |d| <= 1.01 x half a printed unit (P_reaction: max(that, 1e-11)). A value
outside that is 'explained' when the replay (item 18b) decomposes it into measured parts: the
Fortran's own BK coefficients through the C++ analysis, replay(fort.21), reproduce the printed
Fortran value within its half unit, and |C++ - replay(fort.21)| (the propagation difference at
full precision) is within the full-precision criterion (1e-10 absolute). Without the replay
files the predicted bound 2 sqrt(P) dS + dS^2 (dS = 1e-11, criteria.md) is used. Writes data/items/item14_16.json, data/reaction_probability.tsv and
data/state_resolved_dP.tsv."""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from common import DATA, RUN_DIRS, WORK, clean_log_lines, load_item, parse_fortran_real, printed_half_unit, save_item, stats
from parse_logs import energy_blocks, fort25

PRED_DS = 1e-11


def read_tsv(name, k):
    return pd.read_csv(RUN_DIRS[k] / name, sep="\t", comment="#")


def measured_ds():
    r = load_item("item18b")
    out = {}
    for pair in ("cpu-fortran", "gpu-fortran"):
        try:
            out[pair] = r["replay_vs_runs"][pair]["S_max_abs"]
        except (KeyError, TypeError):
            out[pair] = None
    return out


class Tally:
    def __init__(self, ds):
        self.ds = ds
        self.n = 0
        self.n_printed = 0
        self.explained = []
        self.outside = []
        self.worst = {"ratio": -1}
        self.ratios = []

    def add(self, c, ftext, where, floor=0.0, probability=True, replay=None):
        f = parse_fortran_real(ftext)
        tol = max(printed_half_unit(ftext), floor)
        d = abs(c - f)
        ratio = d / tol
        self.n += 1
        self.ratios.append(ratio)
        if ratio > self.worst["ratio"]:
            self.worst = {"ratio": ratio, "where": where, "fortran": ftext, "cpp": c, "abs_diff": d}
        if ratio <= 1:
            self.n_printed += 1
            return "printed"
        rec = {"where": where, "fortran": ftext, "cpp": c, "abs_diff": d, "half_unit": tol}
        if replay is not None:
            rec["replay_fortran"] = replay
            rec["propagation_diff"] = abs(c - replay)
            rec["fortran_print_rounding"] = abs(replay - f)
            ok = abs(replay - f) <= printed_half_unit(ftext) and abs(c - replay) <= 1e-10
            rec["method"] = "replay"
        else:
            bound = 2.0 * math.sqrt(max(f, 0.0)) * self.ds + self.ds ** 2 if probability else 0.0
            rec["bound"] = bound
            ok = d <= bound
            rec["method"] = "predicted bound"
        if ok:
            self.explained.append(rec)
            return "explained"
        self.outside.append(rec)
        return "OUTSIDE"

    def summary(self):
        r = np.array(self.ratios)
        return {"n": self.n, "n_within_printed": self.n_printed, "n_explained": len(self.explained),
                "n_outside": len(self.outside), "worst": self.worst, "dS_used": self.ds,
                "ratio_median": float(np.median(r)) if r.size else None,
                "explained_examples": sorted(self.explained, key=lambda x: -x["abs_diff"])[:5],
                "explained_max_abs_diff": max((x["abs_diff"] for x in self.explained), default=0.0),
                "explained_max_propagation_diff": max((x.get("propagation_diff", 0.0) for x in self.explained), default=0.0),
                "explained_methods": sorted({x["method"] for x in self.explained}),
                "explained_max_value": max((parse_fortran_real(x["fortran"]) for x in self.explained), default=0.0),
                "outside": self.outside[:20]}


def main():
    blocks = energy_blocks(clean_log_lines("fortran"))
    ds_meas = measured_ds()
    summ = {k: read_tsv("scattering_summary.tsv", k) for k in ("cpu", "gpu")}
    rovib = {k: read_tsv("probabilities_rovib.tsv", k) for k in ("cpu", "gpu")}
    diff = {k: read_tsv("probabilities_diffraction.tsv", k) for k in ("cpu", "gpu")}
    E = summ["cpu"]["E_coll"].to_numpy()
    assert len(E) == len(blocks) == 111
    rdir = WORK / "replay" / "build-fast" / "fortran"
    have_replay = (rdir / "s_matrix.tsv").exists()
    if have_replay:
        r_su = pd.read_csv(rdir / "scattering_summary.tsv", sep="\t", comment="#")
        r_rv = pd.read_csv(rdir / "probabilities_rovib.tsv", sep="\t", comment="#").set_index(["E_coll", "v", "j", "mJ"])["P"]
        r_df = pd.read_csv(rdir / "probabilities_diffraction.tsv", sep="\t", comment="#").set_index(["E_coll", "v", "j", "mJ", "order"])["P"]
    rv_ = (lambda col, ie: float(r_su[col][ie])) if have_replay else (lambda col, ie: None)
    res = {"n_energies": len(E), "dS": {"predicted": PRED_DS, "measured": ds_meas},
           "explanation_method": "replay" if have_replay else "predicted bound"}

    # ---------------- item 14
    out14 = {}
    reac_rows = []
    vjm_per_energy = {k: np.zeros(len(E)) for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu")}
    for k in ("cpu", "gpu"):
        ds = ds_meas.get(f"{k}-fortran") or PRED_DS
        T = {name: Tally(ds) for name in ("P_reaction", "P_sum", "P_v", "P_vjm", "E_avg")}
        s = summ[k]
        rv = rovib[k].set_index(["E_coll", "v", "j", "mJ"])["P"]
        for ie, b in enumerate(blocks):
            e = E[ie]
            T["P_reaction"].add(s["P_reaction"][ie], b["p_reaction"], {"E": e}, floor=1e-11, replay=rv_("P_reaction", ie))
            T["P_sum"].add(s["P_scattered"][ie], b["p_sum"], {"E": e}, replay=rv_("P_scattered", ie))
            for v, text in b["p_v"].items():
                T["P_v"].add(s[f"P_v{v}"][ie], text, {"E": e, "v": v}, replay=rv_(f"P_v{v}", ie))
            for col, key in (("E_vib_mean", "e_vib"), ("E_rot_mean", "e_rot"), ("E_dif_mean", "e_dif")):
                T["E_avg"].add(s[col][ie], b[key], {"E": e, "quantity": col}, probability=False, replay=rv_(col, ie))
            for (v, j, mj), text in b["p_vjm"].items():
                c = float(rv.loc[(e, v, j, mj)])
                T["P_vjm"].add(c, text, {"E": e, "v": v, "j": j, "mJ": mj},
                               replay=float(r_rv.loc[(e, v, j, mj)]) if have_replay else None)
                vjm_per_energy[f"{k}-fortran"][ie] = max(vjm_per_energy[f"{k}-fortran"][ie], abs(c - parse_fortran_real(text)))
        out14[f"{k}_vs_fortran"] = {name: t.summary() for name, t in T.items()}
    # per-v average energies printed by the Fortran have no C++ counterpart
    out14["n_per_v_average_energies_without_counterpart"] = sum(2 * len(b["edif_v"]) for b in blocks)
    # CPU <-> GPU, full precision
    gc = {}
    for col in ("P_scattered", "P_reaction", "P_v0", "P_v1", "E_vib_mean", "E_rot_mean", "E_dif_mean"):
        d = np.abs(summ["gpu"][col].to_numpy() - summ["cpu"][col].to_numpy())
        gc[col] = stats(d, where=[float(x) for x in E])
    a, b_ = rovib["gpu"], rovib["cpu"]
    assert (a[["E_coll", "v", "j", "mJ"]].to_numpy() == b_[["E_coll", "v", "j", "mJ"]].to_numpy()).all()
    d = np.abs(a["P"].to_numpy() - b_["P"].to_numpy())
    gc["P_vjm"] = stats(d, where=lambda i: {"E": float(b_["E_coll"][i]), "v": int(b_["v"][i]), "j": int(b_["j"][i]), "mJ": int(b_["mJ"][i])})
    vjm_per_energy["gpu-cpu"] = pd.Series(d).groupby(b_["E_coll"].to_numpy()).max().to_numpy()
    out14["gpu_vs_cpu"] = gc
    res["item14"] = out14

    # reaction probability table (Figure 1, ladder)
    lines = ["# Reaction probability per collision energy. Fortran: printed text (7 significant digits); C++: full "
             "precision (scattering_summary.tsv). half_unit: 1.01 x half a unit in the last printed digit.",
             "E_coll\tP_fortran_text\tP_fortran\thalf_unit\tP_cpu\tP_gpu\tabs_cpu_minus_F\tabs_gpu_minus_F\tabs_gpu_minus_cpu"
             "\tmax_dPvjm_cpu_F\tmax_dPvjm_gpu_F\tmax_dPvjm_gpu_cpu"]
    for ie, b in enumerate(blocks):
        f = parse_fortran_real(b["p_reaction"])
        pc, pg = summ["cpu"]["P_reaction"][ie], summ["gpu"]["P_reaction"][ie]
        lines.append("\t".join([f"{E[ie]:.4f}", b["p_reaction"], f"{f:.10e}", f"{printed_half_unit(b['p_reaction']):.3e}",
                                f"{pc:.16e}", f"{pg:.16e}", f"{abs(pc - f):.3e}", f"{abs(pg - f):.3e}", f"{abs(pg - pc):.3e}",
                                f"{vjm_per_energy['cpu-fortran'][ie]:.3e}", f"{vjm_per_energy['gpu-fortran'][ie]:.3e}",
                                f"{vjm_per_energy['gpu-cpu'][ie]:.3e}"]))
    (DATA / "reaction_probability.tsv").write_text("\n".join(lines) + "\n")
    # how many digits of P_reaction agree (all printed digits?)
    digits = []
    for ie, b in enumerate(blocks):
        f = parse_fortran_real(b["p_reaction"])
        for k in ("cpu", "gpu"):
            c = summ[k]["P_reaction"][ie]
            # C++ value rounded to the Fortran's printed format (D14.7: 7 significant digits)
            digits.append(f"{c:.6e}" == f"{f:.6e}")
    res["p_reaction_same_7_digits"] = {"n": len(digits), "n_same": int(sum(digits))}
    # where the 7 digits differ: distance of the C++ value to the Fortran's rounding boundary
    mism = []
    for ie, b in enumerate(blocks):
        f = parse_fortran_real(b["p_reaction"])
        h = printed_half_unit(b["p_reaction"]) / 1.01
        for k in ("cpu", "gpu"):
            c = summ[k]["P_reaction"][ie]
            if f"{c:.6e}" != f"{f:.6e}":
                mism.append({"E": float(E[ie]), "run": k, "fortran": b["p_reaction"], "cpp": float(c),
                             "distance_to_boundary": float(abs(c - f) - h)})
    res["p_reaction_7digit_mismatches"] = {"cases": mism, "max_distance_to_boundary": max(m["distance_to_boundary"] for m in mism)}
    res["p_reaction_range"] = [float(summ["cpu"]["P_reaction"].min()), float(summ["cpu"]["P_reaction"].max())]

    # ---------------- item 15: diffraction-order probabilities
    out15 = {}
    for k in ("cpu", "gpu"):
        ds = ds_meas.get(f"{k}-fortran") or PRED_DS
        T = Tally(ds)
        dd = diff[k].set_index(["E_coll", "v", "j", "mJ", "order"])["P"]
        missing = 0
        for ie, b in enumerate(blocks):
            for (v, j, mj, order), text in b["p_order"].items():
                key = (E[ie], v, j, mj, order)
                if key not in dd.index:
                    missing += 1
                    continue
                T.add(float(dd.loc[key]), text, {"E": float(E[ie]), "v": v, "j": j, "mJ": mj, "order": order},
                      replay=float(r_df.loc[key]) if have_replay else None)
        out15[f"{k}_vs_fortran"] = {**T.summary(), "fortran_rows_missing_in_cpp": missing,
                                    "cpp_rows": int(len(dd)), "fortran_rows": sum(len(b["p_order"]) for b in blocks)}
    d = np.abs(diff["gpu"]["P"].to_numpy() - diff["cpu"]["P"].to_numpy())
    out15["gpu_vs_cpu"] = stats(d)
    # cpp rows with P > 0 that the Fortran does not print
    fkeys = {(float(E[ie]),) + key for ie, b in enumerate(blocks) for key in b["p_order"]}
    c = diff["cpu"]
    extra = [(float(r.E_coll), int(r.v), int(r.j), int(r.mJ), int(r.order), float(r.P)) for r in c.itertuples()
             if (float(r.E_coll), int(r.v), int(r.j), int(r.mJ), int(r.order)) not in fkeys]
    out15["cpp_rows_not_printed_by_fortran"] = {"n": len(extra), "max_P": max((x[5] for x in extra), default=0.0)}
    res["item15"] = out15

    # ---------------- item 16: S-matrix
    out16 = {}
    smat = {k: read_tsv("s_matrix.tsv", k) for k in ("cpu", "gpu")}
    frows = list(fort25(RUN_DIRS["fortran"] / "fort.25"))
    F = pd.DataFrame(frows, columns=["ie", "v", "j", "mJ", "order", "n", "m", "ReS", "ImS", "P", "phase"])
    F["E_coll"] = E[F["ie"].to_numpy()]
    for col in ("ReS", "ImS", "P", "phase"):
        F[col + "_f"] = F[col].map(parse_fortran_real)
    F["half_P"] = F["P"].map(printed_half_unit)
    F["half_Re"] = F["ReS"].map(printed_half_unit)
    F["half_Im"] = F["ImS"].map(printed_half_unit)
    F["half_ph"] = F["phase"].map(printed_half_unit)
    keys = ["E_coll", "v", "j", "mJ", "n", "m"]
    out16["fortran_rows"] = int(len(F))
    out16["fortran_rows_with_S_zero"] = int(((F["ReS_f"] == 0) & (F["ImS_f"] == 0)).sum())
    for k in ("cpu", "gpu"):
        C = smat[k]
        M = F.merge(C, on=keys, how="outer", indicator=True, suffixes=("", "_c"))
        only_f = M[M["_merge"] == "left_only"]
        only_c = M[M["_merge"] == "right_only"]
        both = M[M["_merge"] == "both"].copy()
        r = {"cpp_rows": int(len(C)), "joined": int(len(both)), "fortran_only": int(len(only_f)),
             "fortran_only_all_S_zero": bool(((only_f["ReS_f"] == 0) & (only_f["ImS_f"] == 0)).all()),
             "cpp_only": int(len(only_c)), "cpp_only_max_P": float(only_c["P_c"].max()) if len(only_c) else 0.0,
             "order_label_equal": bool((both["order"] == both["order_c"]).all())}
        ds = ds_meas.get(f"{k}-fortran") or PRED_DS
        if have_replay:
            RS = pd.read_csv(rdir / "s_matrix.tsv", sep="\t", comment="#")[keys + ["Re_S", "Im_S", "P"]].rename(
                columns={"Re_S": "Re_S_r", "Im_S": "Im_S_r", "P": "P_r"})
            both = both.merge(RS, on=keys, how="left")
        for col, fcol, hcol, ccol in (("Re S", "ReS_f", "half_Re", "Re_S"), ("Im S", "ImS_f", "half_Im", "Im_S"),
                                       ("|S|^2", "P_f", "half_P", "P_c")):
            d = np.abs(both[ccol].to_numpy() - both[fcol].to_numpy())
            ratio = d / both[hcol].to_numpy()
            within = ratio <= 1
            if have_replay:
                rcol = {"Re S": "Re_S_r", "Im S": "Im_S_r", "|S|^2": "P_r"}[col]
                rvals = both[rcol].to_numpy()
                print_ok = np.abs(rvals - both[fcol].to_numpy()) <= both[hcol].to_numpy()
                prop = np.abs(both[ccol].to_numpy() - rvals)
                # rows where replay(fort.21) itself misses the print: analysis-side difference,
                # explained by item 18b's bound (eviroi phase reference, rounding, channel momentum)
                outl = {(o["component"], o["E_coll"], o["v"], o["j"], o["mJ"], o["n"], o["m"]): o
                        for o in load_item("item18b").get("analysis_only", {}).get("analysis_side_S_outliers", [])}
                comp = {"Re S": "Re S", "Im S": "Im S"}.get(col)
                ana_ok = np.array([
                    bool(comp) and (comp, float(e), int(vv), int(jj), int(mm), int(nn), int(m2)) in outl
                    and outl[(comp, float(e), int(vv), int(jj), int(mm), int(nn), int(m2))]["explained"]
                    for e, vv, jj, mm, nn, m2 in zip(both["E_coll"], both["v"], both["j"], both["mJ"], both["n"], both["m"])])
                expl = (~within) & (print_ok | ana_ok) & (prop <= 1e-10)
                unexplained_print = (~within) & ~print_ok & ~ana_ok
            else:
                bound = (2 * np.sqrt(np.maximum(both[fcol].to_numpy(), 0)) * ds + ds ** 2) if col == "|S|^2" else np.full(len(d), ds)
                expl = (~within) & (d <= bound)
                unexplained_print = np.zeros(len(d), dtype=bool)
            i = int(np.argmax(ratio))
            famm = ((both["n"] == 0) & (both["m"] == 0) & (both["mJ"] == 0) & (both["j"] > 0) & (both["v"] == 0)).to_numpy()
            r[col] = {"n": int(len(d)), "n_within_printed": int(within.sum()), "n_explained": int(expl.sum()),
                      "n_explained_in_specular_inelastic_family": int((expl & famm).sum()),
                      "n_outside": int((~within & ~expl).sum()), "worst_ratio": float(ratio[i]),
                      "worst": {kk: (float(both.iloc[i][kk]) if kk == "E_coll" else int(both.iloc[i][kk])) for kk in keys},
                      "worst_values": {"fortran": float(both.iloc[i][fcol]), "cpp": float(both.iloc[i][ccol])},
                      "max_abs_diff_explained": float(d[expl].max()) if expl.any() else 0.0,
                      "n_outside_where_replay_misses_fortran_print": int(unexplained_print.sum()),
                      "max_abs_diff": float(d.max())}
        sel = both["P_f"].to_numpy() > 1e-8
        dph = np.angle(np.exp(1j * (both["phase_c"].to_numpy()[sel] - both["phase_f"].to_numpy()[sel])))
        ratio = np.abs(dph) / both["half_ph"].to_numpy()[sel]
        i = int(np.argmax(ratio))
        r["phase"] = {"n": int(sel.sum()), "n_within_printed": int((ratio <= 1).sum()), "worst_ratio": float(ratio[i]),
                      "max_abs_diff_rad": float(np.abs(dph).max()),
                      "worst": {kk: float(both[sel].iloc[i][kk]) for kk in keys}}
        out16[f"{k}_vs_fortran"] = r
    a, b_ = smat["gpu"], smat["cpu"]
    assert (a[keys].to_numpy() == b_[keys].to_numpy()).all()
    dS = np.abs((a["Re_S"] - b_["Re_S"]).to_numpy() + 1j * (a["Im_S"] - b_["Im_S"]).to_numpy())
    dP = np.abs(a["P"].to_numpy() - b_["P"].to_numpy())
    sel = b_["P"].to_numpy() > 1e-8
    dph = np.abs(np.angle(np.exp(1j * (a["phase"].to_numpy()[sel] - b_["phase"].to_numpy()[sel]))))
    out16["gpu_vs_cpu"] = {"rows": int(len(b_)), "S_abs": stats(dS), "P_abs": stats(dP), "phase_rad_P_gt_1e-8": stats(dph),
                           "identical_rows": int(np.sum((dS == 0)))}
    res["item16"] = out16
    save_item("item14_16", res)

    # ---------------- summary print
    for k in ("cpu", "gpu"):
        o = out14[f"{k}_vs_fortran"]
        print(k, "item14:", {n: (v["n"], v["n_within_printed"], v["n_explained"], v["n_outside"], round(v["worst"]["ratio"], 3)) for n, v in o.items()})
        o = out15[f"{k}_vs_fortran"]
        print(k, "item15:", o["n"], o["n_within_printed"], o["n_explained"], o["n_outside"], round(o["worst"]["ratio"], 2), "missing", o["fortran_rows_missing_in_cpp"])
        o = out16[f"{k}_vs_fortran"]
        print(k, "item16 rows:", o["cpp_rows"], o["joined"], o["fortran_only"], o["fortran_only_all_S_zero"], o["cpp_only"])
        for col in ("Re S", "Im S", "|S|^2"):
            print("   ", col, o[col]["n"], o[col]["n_within_printed"], o[col]["n_explained"], o[col]["n_outside"], round(o[col]["worst_ratio"], 2))
        print("    phase", o["phase"])
    print("GPU-CPU item14:", {c: f"{v['max']:.2e}" for c, v in out14["gpu_vs_cpu"].items()})
    print("GPU-CPU item15:", out15["gpu_vs_cpu"]["max"], "item16 S:", out16["gpu_vs_cpu"]["S_abs"]["max"], "P:", out16["gpu_vs_cpu"]["P_abs"]["max"])
    print("P_reaction same 7 digits:", res["p_reaction_same_7_digits"], "range", res["p_reaction_range"])


if __name__ == "__main__":
    main()
