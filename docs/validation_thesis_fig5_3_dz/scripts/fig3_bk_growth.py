#!/usr/bin/env python3
"""Figure 3: BK coefficient differences against time, all code pairs. (a) the gate-4 metric
max over channels |c_a - c_b| / max(max|c_F|, 1e-4) with the criterion 1e-10; (b) the amplitude
part kappa(t) of the difference, c_a = c_b (1 + kappa + i beta) on average (|c_b|^2 weighted),
where ||c(t)|| >= 0.1% of its maximum, with its median per quarter of the signal window
(item13_propagation.py). The strict-build GPU run against the Fortran is added (2026-11-08,
data/strict_bk_error_vs_time.tsv from strict_gpu.py). Reads data/bk_error_vs_time.tsv,
data/strict_bk_error_vs_time.tsv, data/items/item13.json."""
from __future__ import annotations

import numpy as np

from common import DATA, FIGURES, load_item
from style import PAIR, INK_2, MUTED, plt, save, setup


def main():
    setup()
    lines = (DATA / "bk_error_vs_time.tsv").read_text().splitlines()
    head = lines[1].split("\t")
    a = np.array([[float(x) for x in l.split("\t")] for l in lines[2:]])
    col = {h: a[:, i] for i, h in enumerate(head)}
    t = col["t"]
    sl = (DATA / "strict_bk_error_vs_time.tsv").read_text().splitlines()
    sh = sl[1].split("\t")
    sa = np.array([[float(x) for x in l.split("\t")] for l in sl[2:]])
    scol = {h: sa[:, i] for i, h in enumerate(sh)}
    assert np.array_equal(scol["t"], t)
    it = load_item("item13")
    b = it["bk"]
    t_fix = it["analysis_line"]["time_of_that_record"]["fortran"]
    w0, w1 = b["signal_window"]["t_first"], b["signal_window"]["t_last"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.8, 7.4), sharex=True, gridspec_kw={"hspace": 0.1})
    for ax in (ax1, ax2):
        ax.axvspan(w0, w1, color="#f1f0ec", lw=0)
        ax.axvline(t_fix, color=MUTED, lw=1.0, ls=(0, (2, 2)))
    for key in ("cpu-fortran", "gpu-fortran", "gpu-cpu"):
        st = PAIR[key]
        ax1.plot(t, np.maximum(col[key], 1e-18), color=st["color"], ls=st["ls"], lw=1.3, label=st["label"])
        sel = col["norm_c_t"] >= 1e-3 * col["norm_c_t"].max()     # where the signal is not noise
        a_ = np.where(sel, 1e12 * col[f"{key}_alpha"], np.nan)
        ax2.plot(t, a_, color=st["color"], ls=st["ls"], lw=1.3, label=st["label"])
    st = PAIR["gpu_strict-fortran"]
    ax1.plot(t, np.maximum(scol["strict-fortran"], 1e-18), color=st["color"], ls=st["ls"], lw=1.3, label=st["label"] + " (added 2026-11-08)")
    ax2.plot(t, np.where(sel, 1e12 * scol["strict-fortran_alpha"], np.nan), color=st["color"], ls=st["ls"], lw=1.3,
             label=st["label"] + " (added 2026-11-08)")
    edges = np.linspace(w0, w1, 5)
    for key in ("cpu-fortran", "gpu-cpu"):
        q = b["amplitude_vs_time"][key]["window_quartile_medians"]
        for i in range(4):
            ax2.plot([edges[i], edges[i + 1]], [1e12 * q[i]] * 2, color=PAIR[key]["color"], lw=4, alpha=0.4,
                     label="median per quarter of the window" if (i == 0 and key == "cpu-fortran") else None)
    ax1.axhline(1e-10, color=INK_2, lw=1.0, ls=(0, (6, 3)))
    ax1.text(200, 1.4e-10, "criterion 1e-10 (gate 4)", fontsize=8, color=INK_2)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-17, 1e-8)
    ax1.set_ylabel("max |Δc(t)| / max|c|")
    ax1.legend(loc="lower right", fontsize=8.5)
    ax1.set_title("BK coefficients on the analysis line: differences between the codes in time", loc="left")
    ax1.text(t_fix + 120, 2e-9, f"line fixed, t = {t_fix:.0f}", fontsize=8, color=INK_2)
    ax1.text(w1 - 120, 2e-9, "signal window: ||c(t)|| ≥ 1% of max", fontsize=8, color=INK_2, ha="right")
    ax2.axhline(0, color=MUTED, lw=0.8)
    ax2.set_ylim(-4, 16)
    ax2.set_ylabel("amplitude part κ(t)  (×10⁻¹²)")
    ax2.text(0.01, 0.97, "shown where ||c(t)|| ≥ 0.1% of its maximum", transform=ax2.transAxes, fontsize=8, color=INK_2,
             va="top")
    ax2.set_xlabel("Time t (a.u.)")
    ax2.set_xlim(0, 17500)
    ax2.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92), fontsize=8)
    save(fig, FIGURES / "fig3_bk_error_growth.png")
    print("fig3 written")


if __name__ == "__main__":
    main()
