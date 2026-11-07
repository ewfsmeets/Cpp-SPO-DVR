#!/usr/bin/env python3
"""Item 13: the propagation. (a) averages every dtAve; (b) the step at which the specular
analysis line is fixed; (c) the BK coefficients of every analysis step, all pairs, in time;
(d, with --psi) the final wave function, all pairs.

Writes data/items/item13.json (and item13_psi.json with --psi) and data/bk_error_vs_time.tsv
(the time series behind Figure 3)."""
from __future__ import annotations

import struct
import sys

import numpy as np

from common import (DATA, PAIRS, RUN_DIRS, FortranBk, NativeBk, clean_log_lines, fortran_record_offsets,
                    parse_fortran_real, printed_half_unit, read_record, save_item, stats)
from parse_logs import averages_out, cpp_surface_fixed, fortran_surface_fixed

COLS = ["t", "norm", "Z", "r", "rest_norm"]


def averages():
    f = averages_out(RUN_DIRS["fortran"] / "averages.out")
    cpp = {k: [l.split() for l in (RUN_DIRS[k] / "averages.dat").read_text().splitlines() if not l.startswith("#")]
           for k in ("cpu", "gpu")}
    res = {"n_rows": {"fortran": len(f), "cpu": len(cpp["cpu"]), "gpu": len(cpp["gpu"])}}
    for k in ("cpu", "gpu"):
        per_col = {}
        worst = {"ratio": -1}
        n_out = 0
        for frow, crow in zip(f, cpp[k]):
            for j, col in enumerate(COLS):
                if col == "t":
                    continue      # the time stamps are the same by construction; not counted
                ft, c = frow[j], float(crow[j])
                tol = printed_half_unit(ft)
                ratio = abs(c - parse_fortran_real(ft)) / tol
                per_col.setdefault(col, []).append(ratio)
                if ratio > worst["ratio"]:
                    worst = {"ratio": ratio, "t": float(crow[0]), "column": col, "fortran": ft, "cpp": c}
                n_out += ratio > 1
        res[f"{k}_vs_fortran"] = {"n_values": len(f) * (len(COLS) - 1), "n_times": len(f), "n_outside": int(n_out), "worst": worst,
                                  "worst_ratio_by_column": {c: max(v) for c, v in per_col.items()}}
    a = np.array([[float(x) for x in r] for r in cpp["gpu"]])
    b = np.array([[float(x) for x in r] for r in cpp["cpu"]])
    gc = {}
    for j, col in enumerate(COLS):
        if col in ("rest_norm", "t"):
            d = np.abs(a[:, j] - b[:, j])
            gc[col] = {**stats(d, where=[float(t) for t in b[:, 0]]), "metric": "absolute"}
        else:
            d = np.abs(a[:, j] - b[:, j]) / np.abs(b[:, j])
            gc[col] = {**stats(d, where=[float(t) for t in b[:, 0]]), "metric": "relative"}
    res["gpu_vs_cpu"] = gc
    res["final"] = {"t": b[-1, 0], "norm_cpu": b[-1, 1], "norm_gpu": a[-1, 1], "norm_fortran_text": f[-1][1]}
    # norm absorbed vs left (for the ladder: norm left at T)
    res["norm_left_at_T"] = float(b[-1, 1])
    return res


def analysis_line(bk):
    cz, ct = cpp_surface_fixed(clean_log_lines("cpu"))
    gz, gt = cpp_surface_fixed(clean_log_lines("gpu"))
    fz = fortran_surface_fixed(clean_log_lines("fortran"))
    init = 4320  # initial rot-dif channel (0-based), v = 0; verified in setup.json (inchna)
    series = {k: np.abs(bk[k][:, 0, init]) for k in ("fortran", "cpu", "gpu")}
    t = bk["times"]
    def first_jump(s):
        # first record where the initial-channel coefficient exceeds 1e-10 (before: ~1e-19 noise)
        return int(np.argmax(s > 1e-10))
    jumps = {k: first_jump(series[k]) for k in series}
    j = jumps["cpu"]
    window = {k: [float(x) for x in series[k][j - 3:j + 3]] for k in series}
    return {"cpp_cpu": {"Z": cz, "t": ct}, "cpp_gpu": {"Z": gz, "t": gt}, "fortran": {"Z_text": fz},
            "first_record_with_specular_projection": {k: int(v) for k, v in jumps.items()},
            "time_of_that_record": {k: float(t[v]) for k, v in jumps.items()},
            "initial_channel_abs_around_fix": {"t": [float(x) for x in t[j - 3:j + 3]], **window},
            "max_abs_before_fix": {k: float(series[k][:j].max()) for k in series},
            "abs_at_fix": {k: float(series[k][j]) for k in series},
            "pair_diff_at_fix": {f"{a}-{b}": float(abs(bk[a][j, 0, init] - bk[b][j, 0, init])) for a, b in PAIRS}}


def load_bk():
    F = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    C = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    G = NativeBk(RUN_DIRS["gpu"] / "bk_coefficients.bin")
    assert F.n_records == C.n_records == G.n_records
    tc, tg = C.times(), G.times()
    assert np.array_equal(tc, tg) and np.allclose(tc, F.times(), rtol=0, atol=1e-9)
    bk = {"fortran": np.stack([F.record(k) for k in range(F.n_records)]),
          "cpu": np.stack([C.record(k) for k in range(C.n_records)]),
          "gpu": np.stack([G.record(k) for k in range(G.n_records)]),
          "times": tc, "labels": C.labels, "nstata": C.nstata, "times_fortran_minus_cpp_max": float(np.max(np.abs(F.times() - tc)))}
    return bk


def fit_growth(t, e, t_min):
    sel = (t >= t_min) & (e > 0)
    x, y = np.log(t[sel]), np.log(e[sel])
    b, a = np.polyfit(x, y, 1)
    return {"exponent": float(b), "prefactor": float(np.exp(a)), "t_min": float(t_min), "n": int(sel.sum())}


def bk_compare(bk, t_fix):
    labels = bk["labels"]
    ns = bk["nstata"]
    t = bk["times"]
    absmax = {k: float(np.abs(bk[k]).max()) for k in ("fortran", "cpu", "gpu")}
    scale = max(absmax["fortran"], 1e-4)
    inst = np.abs(bk["fortran"]).reshape(len(t), -1).max(axis=1)   # max |c(t)| per record
    res = {"n_records": int(len(t)), "n_coefficients_per_record": int(bk["fortran"][0].size), "max_abs_c": absmax,
           "scale": scale, "criterion": 1e-10 * scale, "times_fortran_minus_cpp_max": bk["times_fortran_minus_cpp_max"]}
    series = {"t": t, "max_abs_c_t": inst}
    l2 = np.sqrt((np.abs(bk["fortran"]) ** 2).reshape(len(t), -1).sum(axis=1))   # ||c_F(t)||_2 over all channels
    series["norm_c_t"] = l2
    window = l2 >= 1e-2 * l2.max()      # main scattered flux: the signal is not noise
    res["signal_window"] = {"criterion": "||c_F(t)|| >= 1e-2 max_t ||c_F||", "t_first": float(t[window][0]),
                            "t_last": float(t[window][-1]), "n": int(window.sum()), "max_norm_c": float(l2.max()),
                            "t_of_max": float(t[int(np.argmax(l2))])}
    chan_max = np.abs(bk["fortran"]).max(axis=0)       # per channel (v, is): max over time
    floor = 1e-6 * absmax["fortran"]
    big = chan_max > floor
    for a, b in PAIRS:
        d = np.abs(bk[a] - bk[b])
        k = int(np.argmax(d))
        rec, v, s = np.unravel_index(k, d.shape)
        e_t = d.reshape(len(t), -1).max(axis=1)
        series[f"{a}-{b}"] = e_t / scale
        series[f"{a}-{b}_rel_inst"] = np.where(inst > 0, e_t / np.where(inst > 0, inst, 1), 0)
        ch_rel = d.max(axis=0)[big] / chan_max[big]
        kk = int(np.argmax(ch_rel))
        vv, ss = np.unravel_index(np.flatnonzero(big)[kk], big.shape)
        flat = d.reshape(-1)
        dl2 = np.sqrt((d ** 2).reshape(len(t), -1).sum(axis=1))
        rel_l2 = dl2 / np.where(l2 > 0, l2, 1)
        series[f"{a}-{b}_l2"] = dl2
        series[f"{a}-{b}_rel_l2"] = rel_l2
        jfix = int(np.argmin(np.abs(t - t_fix)))
        steps = dl2[1:] / np.where(dl2[:-1] > 0, dl2[:-1], 1)
        bw, aw = np.polyfit(np.log(t[window]), np.log(rel_l2[window]), 1)
        res[f"{a}-{b}"] = {
            "max_abs": float(d.max()), "max_over_scale": float(d.max() / scale),
            "median_abs": float(np.median(flat)), "p95_abs": float(np.percentile(flat, 95)),
            "argmax": {"t": float(t[rec]), "v": int(v), "is": int(s), "J": int(labels["J"][s]), "mJ": int(labels["mJ"][s]),
                       "n": int(labels["n"][s]), "m": int(labels["m"][s]), "abs_c": float(abs(bk[b][rec, v, s]))},
            "within_criterion": bool(d.max() <= 1e-10 * scale),
            "per_channel_relative": {"floor": floor, "n_channels": int(big.sum()), "median": float(np.median(ch_rel)),
                                     "p95": float(np.percentile(ch_rel, 95)), "max": float(ch_rel.max()),
                                     "argmax": {"v": int(vv), "is": int(ss), "J": int(labels["J"][ss]),
                                                "mJ": int(labels["mJ"][ss]), "n": int(labels["n"][ss]),
                                                "m": int(labels["m"][ss]), "max_abs_c": float(chan_max[vv, ss])}},
            "growth_over_scale": fit_growth(t, e_t / scale, t_fix),
            "growth_relative_to_instantaneous_max": fit_growth(t, series[f"{a}-{b}_rel_inst"], t_fix),
            "identical": bool(np.array_equal(bk[a], bk[b])),
            "l2": {"rel_in_window": stats(rel_l2[window], where=[float(x) for x in t[window]]),
                   "abs_before_fix": stats(dl2[:jfix], where=[float(x) for x in t[:jfix]]),
                   "abs_after_window": stats(dl2[t > t[window][-1]], where=[float(x) for x in t[t > t[window][-1]]]),
                   "growth_rel_in_window": {"exponent": float(bw), "prefactor": float(np.exp(aw))},
                   "jump_at_fix": {"before": float(dl2[jfix - 1]), "at": float(dl2[jfix]), "ratio": float(dl2[jfix] / dl2[jfix - 1]),
                                   "ratio_max_metric": float(e_t[jfix] / e_t[jfix - 1])},
                   "max_step_ratio": float(steps.max()), "t_of_max_step_ratio": float(t[1:][int(np.argmax(steps))]),
                   "min_step_ratio": float(steps.min())},
        }
    # amplitude and phase parts of the difference per analysis step: c_a = c_F (1 + alpha + i beta)
    # on average (|c_F|^2 weighted). A per-step norm bias gives alpha(t) proportional to t.
    cf = bk["fortran"].reshape(len(t), -1)
    amp = {}
    for a, b in PAIRS:
        ref = bk[b].reshape(len(t), -1)
        da = bk[a].reshape(len(t), -1) - ref
        z = np.sum(np.conj(ref) * da, axis=1) / np.where(l2 > 0, np.sum(np.abs(ref) ** 2, axis=1), 1)
        series[f"{a}-{b}_alpha"] = z.real
        series[f"{a}-{b}_beta"] = z.imag
        w = window
        slope, icpt = np.polyfit(t[w], z.real[w], 1)
        y = z.real[w]
        def r2(pred):
            return float(1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2))
        ww = l2[w] ** 2                                    # |c_F(t)|^2 weights
        sw, iw = np.polyfit(t[w], y, 1, w=np.sqrt(ww))
        qs = np.array_split(y, 4)
        through0 = float(np.sum(t[w] * z.real[w]) / np.sum(t[w] ** 2))
        resid = z.real[w] - (icpt + slope * t[w])
        amp[f"{a}-{b}"] = {"alpha_window_median": float(np.median(z.real[w])), "alpha_window_min": float(z.real[w].min()),
                           "alpha_window_max": float(z.real[w].max()), "beta_window_abs_max": float(np.abs(z.imag[w]).max()),
                           "fit_slope_per_au": float(slope), "fit_intercept": float(icpt),
                           "fit_rms_residual": float(np.sqrt(np.mean(resid ** 2))),
                           "fit_r2": r2(icpt + slope * t[w]),
                           "fit_weighted_slope_per_au": float(sw), "fit_weighted_intercept": float(iw),
                           "fit_weighted_r2_unweighted_metric": r2(iw + sw * t[w]),
                           "fit_through_origin_r2": r2(through0 * t[w]),
                           "constant_r2": 0.0, "window_quartile_medians": [float(np.median(q)) for q in qs],
                           "slope_through_origin_per_au": through0,
                           "per_kinetic_block_amplitude_bias": through0 * 20.0 / 9.0,   # 9 kinetic blocks per 20 a.u.
                           "alpha_at_window_ends": [float(z.real[w][0]), float(z.real[w][-1])],
                           "t_window": [float(t[w][0]), float(t[w][-1])]}
    res["amplitude_vs_time"] = amp
    # channel families: the specular rotationally inelastic channels (v = 0, (n, m) = (0, 0),
    # mJ = 0, J > 0) carry the angular-transform image of the large specular component
    J, mJ, n_, m_ = (np.asarray(labels[k]) for k in ("J", "mJ", "n", "m"))
    fam = (n_ == 0) & (m_ == 0) & (mJ == 0) & (J > 0)
    jf = int(np.argmin(np.abs(t - t_fix)))
    famres = {}
    for a, b in PAIRS:
        d2 = np.abs(bk[a][:, 0, :] - bk[b][:, 0, :]) ** 2           # v = 0
        dall = (np.abs(bk[a] - bk[b]) ** 2).sum(axis=(1, 2))
        per_j = {}
        for jj in sorted(set(J[fam].tolist())):
            k = np.flatnonzero(fam & (J == jj))[0]
            dc = np.abs(bk[a][:, 0, k] - bk[b][:, 0, k])
            cmax = np.abs(bk[b][:, 0, k]).max()
            per_j[int(jj)] = {"max_abs_dc": float(dc.max()), "max_abs_c": float(cmax), "rel": float(dc.max() / cmax)}
        famres[f"{a}-{b}"] = {"fraction_of_sum_dc2_before_fix": float(d2[:jf][:, fam].sum() / dall[:jf].sum()),
                              "fraction_of_sum_dc2_all_t": float(d2[:, fam].sum() / dall.sum()),
                              "per_J": per_j}
    res["specular_inelastic_family"] = famres
    # before the fixing step vs after: ||dc|| of the pairs
    jf = int(np.argmin(np.abs(t - t_fix)))
    pre = series["cpu-fortran_l2"][:jf] / np.maximum(series["gpu-cpu_l2"][:jf], 1e-300)
    post = (t > t_fix) & (t < 5000)
    res["pre_post_fix"] = {"ratio_cpuF_over_gpuC_before_fix": {"min": float(pre[5:].min()), "max": float(pre[5:].max()),
                                                               "median": float(np.median(pre[5:]))},
                           "l2_2240_5000": {k: {"median": float(np.median(series[f"{k}_l2"][post])),
                                                "max": float(series[f"{k}_l2"][post].max())} for k in ("cpu-fortran", "gpu-cpu")}}
    res["n_channels_total"] = int(big.size)
    res["n_channels_above_floor"] = int(big.sum())
    # time series for Figure 3
    cols = list(series)
    lines = ["# BK coefficient differences per analysis step: max over channels |c_a - c_b| / scale (scale = "
             f"max(max|c_Fortran|, 1e-4) = {scale:.6e}); *_rel_inst: / max|c_Fortran(t)|; max_abs_c_t = max|c_Fortran(t)|",
             "\t".join(cols)]
    for i in range(len(t)):
        lines.append("\t".join(f"{series[c][i]:.6e}" for c in cols))
    (DATA / "bk_error_vs_time.tsv").write_text("\n".join(lines) + "\n")
    return res


# ----------------------------------------------------------------------------- final wave function
def checkpoint_layout(path):
    with open(path, "rb") as f:
        assert f.read(16).startswith(b"SPODVR-CHECKPNT")
        f.read(4)
        (n,) = struct.unpack("<Q", f.read(8)); f.read(n)
        (count,) = struct.unpack("<Q", f.read(8))
        for _ in range(count):
            (n,) = struct.unpack("<Q", f.read(8)); f.read(n + 8)
        (time,) = struct.unpack("<d", f.read(8))
        nrot, nxy, nzr, nr, nzsp = struct.unpack("<5i", f.read(20))
        idx, prev = struct.unpack("<2i", f.read(8))
        fixed = f.read(1)[0]
        (tol,) = struct.unpack("<d", f.read(8))
        spec_off = f.tell()
    scat_off = spec_off + 16 * nr * nzsp
    return {"time": time, "nrot": nrot, "nxy": nxy, "nzr": nzr, "nr": nr, "nzsp": nzsp, "surface_index": idx,
            "surface_previous": prev, "surface_fixed": fixed, "anly_tol": tol, "spec_off": spec_off, "scat_off": scat_off}


def psi_compare():
    lay = {k: checkpoint_layout(RUN_DIRS[k] / "checkpoint.bin") for k in ("cpu", "gpu")}
    L = lay["cpu"]
    nblock = L["nrot"] * L["nxy"]
    recs = fortran_record_offsets(RUN_DIRS["fortran"] / "fort.9")
    assert len(recs) == 3 + L["nzr"]
    nproj = struct.unpack("<i", read_record(RUN_DIRS["fortran"] / "fort.9", *recs[0]))[0]
    fixed, izan = struct.unpack("<2i", read_record(RUN_DIRS["fortran"] / "fort.9", *recs[2]))
    stride = recs[4][0] - recs[3][0]
    fbase = np.memmap(RUN_DIRS["fortran"] / "fort.9", dtype="<c16", mode="r", offset=recs[3][0],
                      shape=((L["nzr"] - 1) * stride // 16 + nblock,))
    fscat = np.lib.stride_tricks.as_strided(fbase, shape=(L["nzr"], nblock), strides=(stride, 16))
    fspec = np.frombuffer(read_record(RUN_DIRS["fortran"] / "fort.9", *recs[1]), dtype="<c16")
    psi = {"fortran": fscat}
    spec = {"fortran": fspec}
    for k in ("cpu", "gpu"):
        lk = lay[k]
        psi[k] = np.memmap(RUN_DIRS[k] / "checkpoint.bin", dtype="<c16", mode="r", offset=lk["scat_off"],
                           shape=(lk["nzr"], nblock))
        spec[k] = np.fromfile(RUN_DIRS[k] / "checkpoint.bin", dtype="<c16", count=lk["nr"] * lk["nzsp"], offset=lk["spec_off"])
    norms = {k: 0.0 for k in psi}
    acc = {f"{a}-{b}": {"d2": 0.0, "overlap": 0.0 + 0.0j, "max_block_rel": 0.0, "argmax_izr": -1} for a, b in PAIRS}
    nr = L["nr"]
    for izr in range(L["nzr"]):
        blocks = {k: np.asarray(psi[k][izr]) for k in psi}
        for k in blocks:
            norms[k] += float(np.vdot(blocks[k], blocks[k]).real)
        for a, b in PAIRS:
            d = blocks[a] - blocks[b]
            d2 = float(np.vdot(d, d).real)
            r = acc[f"{a}-{b}"]
            r["d2"] += d2
            r["overlap"] += np.vdot(blocks[b], blocks[a])
            nb2 = float(np.vdot(blocks[b], blocks[b]).real)
            if nb2 > 1e-30 and np.sqrt(d2 / nb2) > r["max_block_rel"] and nb2 > 1e-12 * 1e-6:
                r["max_block_rel"], r["argmax_izr"] = float(np.sqrt(d2 / nb2)), izr
    res = {"surface": {"cpu": {k: lay["cpu"][k] for k in ("surface_index", "surface_fixed", "time")},
                       "gpu": {k: lay["gpu"][k] for k in ("surface_index", "surface_fixed", "time")},
                       "fortran": {"fixed": fixed, "izanlySpecChan_1based": izan, "nproj": nproj}},
           "squared_norms_scattering_grid_sum": norms}
    for a, b in PAIRS:
        r = acc[f"{a}-{b}"]
        ref = np.sqrt(norms[b])
        ov = r["overlap"]
        izr = r["argmax_izr"]
        res[f"{a}-{b}"] = {"rel_l2": float(np.sqrt(r["d2"]) / ref), "global_phase_rad": float(np.angle(ov)),
                           "overlap_modulus_minus_1": float(abs(ov) / np.sqrt(norms[a] * norms[b]) - 1.0),
                           "max_block_rel": r["max_block_rel"], "argmax": {"iz": izr // nr, "ir": izr % nr} if izr >= 0 else None,
                           "specular_rel_l2": float(np.linalg.norm(spec[a] - spec[b]) / max(np.linalg.norm(spec[b]), 1e-300)),
                           "specular_norm": float(np.linalg.norm(spec[b]))}
    return res


def derive_psi():
    """Physical norms at T (from the averages and the specular vector) and the psi differences
    relative to the initial norm (1): the final psi is the small residual left on the grids."""
    from common import INPUT_DIR, load_item
    import json
    r = load_item("item13_psi")
    if not r:
        return None
    js = json.loads((INPUT_DIR / "input.json").read_text())
    last = [l.split() for l in (RUN_DIRS["cpu"] / "averages.dat").read_text().splitlines() if not l.startswith("#")][-1]
    total = float(last[1])
    spec = r["cpu-fortran"]["specular_norm"] ** 2 * js["dZ"] * js["dR"]
    scat = total - spec
    r["norm_at_T"] = {"total": total, "specular": spec, "scattering": scat}
    nscat = r["squared_norms_scattering_grid_sum"]
    spec_norm = {"fortran": r["cpu-fortran"]["specular_norm"], "cpu": r["gpu-cpu"]["specular_norm"]}
    for a, b in PAIRS:
        q = r[f"{a}-{b}"]
        q["scattering_amplitude_ratio_minus_1"] = float(np.sqrt(nscat[a] / nscat[b]) - 1)
        if a in spec_norm and b in spec_norm:
            q["specular_amplitude_ratio_minus_1"] = float(spec_norm[a] / spec_norm[b] - 1)
        q["scattering_abs_rel_initial_norm"] = q["rel_l2"] * np.sqrt(scat)
        q["specular_abs_rel_initial_norm"] = q["specular_rel_l2"] * np.sqrt(spec)
    save_item("item13_psi", r)
    return r


def main():
    if "--derive-psi" in sys.argv:
        r = derive_psi()
        print(r["norm_at_T"], {f"{a}-{b}": (r[f"{a}-{b}"]["scattering_abs_rel_initial_norm"], r[f"{a}-{b}"]["specular_abs_rel_initial_norm"]) for a, b in PAIRS})
        return
    if "--psi" in sys.argv:
        res = psi_compare()
        save_item("item13_psi", res)
        derive_psi()
        print({k: v for k, v in res.items() if k not in PAIRS and not isinstance(k, tuple)})
        for a, b in PAIRS:
            print(a, b, res[f"{a}-{b}"])
        return
    res = {"averages": averages()}
    bk = load_bk()
    res["analysis_line"] = analysis_line(bk)
    res["bk"] = bk_compare(bk, res["analysis_line"]["time_of_that_record"]["cpu"])
    save_item("item13", res)
    av = res["averages"]
    for k in ("cpu", "gpu"):
        print(k, "vs F averages:", av[f"{k}_vs_fortran"]["n_outside"], "outside of", av[f"{k}_vs_fortran"]["n_values"],
              "worst", {kk: round(vv, 3) for kk, vv in av[f"{k}_vs_fortran"]["worst_ratio_by_column"].items()})
    print("GPU vs CPU averages:", {c: (v["max"], v["metric"]) for c, v in av["gpu_vs_cpu"].items()})
    al = res["analysis_line"]
    print("line:", al["cpp_cpu"], al["cpp_gpu"], al["fortran"], al["first_record_with_specular_projection"],
          al["max_abs_before_fix"], al["abs_at_fix"])
    b = res["bk"]
    print("BK scale", b["scale"], "max|c|", b["max_abs_c"], "channels above floor", b["n_channels_above_floor"], "of", b["n_channels_total"])
    for a, bb in PAIRS:
        r = b[f"{a}-{bb}"]
        print(f" {a}-{bb}: max {r['max_abs']:.2e} = {r['max_over_scale']:.2e} of scale at {r['argmax']}; "
              f"per-channel rel max {r['per_channel_relative']['max']:.2e} median {r['per_channel_relative']['median']:.2e}; "
              f"growth exp {r['growth_over_scale']['exponent']:.2f}")
        l = r["l2"]
        print(f"    l2: rel in window median {l['rel_in_window']['median']:.2e} max {l['rel_in_window']['max']:.2e}; growth exp "
              f"{l['growth_rel_in_window']['exponent']:.2f}; abs before fix median {l['abs_before_fix']['median']:.1e}; "
              f"jump at fix x{l['jump_at_fix']['ratio']:.2f}; max step ratio {l['max_step_ratio']:.2f} at t={l['t_of_max_step_ratio']}")
    print("signal window", b["signal_window"])
    for k, v in b["amplitude_vs_time"].items():
        print("alpha(t)", k, {kk: (f"{vv:.3e}" if isinstance(vv, float) else vv) for kk, vv in v.items()})
    print("pre/post fix", b["pre_post_fix"])
    for k, v in b["specular_inelastic_family"].items():
        print("family", k, "before fix", round(v["fraction_of_sum_dc2_before_fix"], 4), "all t", round(v["fraction_of_sum_dc2_all_t"], 4),
              {j: f"{x['rel']:.1e} (|dc| {x['max_abs_dc']:.1e}, |c| {x['max_abs_c']:.1e})" for j, x in v["per_J"].items()})


if __name__ == "__main__":
    main()
