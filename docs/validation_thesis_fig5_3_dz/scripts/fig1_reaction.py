#!/usr/bin/env python3
"""Figure 1: reaction probability against collision energy, the three runs, the strict-build GPU run
(added 2026-11-08) and the thesis curve (top), and the differences between the codes on a log scale
(bottom). Reads only data/ extracts and data/items/*.json."""
from __future__ import annotations

import json

import numpy as np

from common import DATA, FIGURES, INPUT_DIR, load_item
from style import PAIR, PROGRAM, INK_2, MUTED, plt, save, setup


def read_tsv(path, skip):
    rows = [l.split("\t") for l in path.read_text().splitlines()[skip:]]
    return rows


def main():
    setup()
    rp = read_tsv(DATA / "reaction_probability.tsv", 2)
    E = np.array([float(r[0]) for r in rp])
    Pf = np.array([float(r[2]) for r in rp])
    half = np.array([float(r[3]) for r in rp])
    Pc = np.array([float(r[4]) for r in rp])
    Pg = np.array([float(r[5]) for r in rp])
    ftext = [r[1] for r in rp]
    rep = read_tsv(DATA / "replay_reaction_probability.tsv", 2)
    d_cf = np.array([float(r[4]) for r in rep])
    d_gf = np.array([float(r[5]) for r in rep])
    d_gc = np.abs(Pg - Pc)
    srp = read_tsv(DATA / "strict_reaction_probability.tsv", 2)
    Ps = np.array([float(r[1]) for r in srp])
    d_sf = np.array([float(r[4]) for r in srp])
    assert np.allclose([float(r[0]) for r in srp], E, rtol=0, atol=1e-9)
    th = np.loadtxt(DATA / "thesis_fig5_3_dz_digitised.tsv", comments="#", skiprows=3)
    i14 = load_item("item14_16")
    i17 = load_item("item17")
    p0 = load_item("phase0")
    js = json.loads((INPUT_DIR / "input.json").read_text())

    fig = plt.figure(figsize=(8.8, 10.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], width_ratios=[1.0, 0.36], hspace=0.12, wspace=0.03)
    ax = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[1, 0], sharex=ax)
    axt = fig.add_subplot(gs[:, 1])
    axt.set_axis_off()

    # ---- top: P_reaction
    ax.plot(th[:, 0], th[:, 1], color=PROGRAM["thesis"]["color"], lw=0, marker="D", ms=3.2, alpha=0.75,
            label=PROGRAM["thesis"]["label"])
    for k, P, off in (("fortran", Pf, 0), ("cpu", Pc, 2), ("gpu", Pg, 4), ("gpu_strict", Ps, 6)):
        st = PROGRAM[k]
        ax.plot(E, P, color=st["color"], ls=st["ls"], lw=1.6, marker=st["marker"], ms=5, markevery=(off, 8),
                mfc="white" if k != "fortran" else st["color"], label=st["label"])
    ax.set_ylabel("Reaction probability P_reaction")
    ax.set_xlim(0.395, 1.0)
    ax.set_ylim(-0.01, 0.68)
    ax.legend(loc="upper left", fontsize=8.5)
    nsame = i14["p_reaction_same_7_digits"]
    n_strict_same = int(sum(f"{s:.6e}" == f"{f:.6e}" for s, f in zip(Ps, Pf)))
    full = max(i14["item14"]["gpu_vs_cpu"]["P_reaction"]["max"], max(d_cf), max(d_gf))
    ax.text(0.015, 0.47,
            f"The curves coincide: |ΔP_reaction| ≤ {full:.1e} between any two\n"
            f"codes (propagation, full precision), far below one line width (~1e-3).\n"
            f"C++ rounded to the Fortran's 7 printed digits: {nsame['n_same']}/{nsame['n']} identical\n"
            f"(strict GPU build: {n_strict_same}/{E.size}, |ΔP_reaction| ≤ {d_sf.max():.1e} from the Fortran).\n"
            f"Thesis curve: agrees at {i17['cpu']['n_within']}/{i17['n_points_in_range']} digitised points "
            f"(max |ΔP| {i17['cpu']['max_abs_dP']:.3f}).",
            transform=ax.transAxes, ha="left", va="center", fontsize=8.4,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=MUTED, lw=0.6))
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title("H₂(v=0, j=0) + Cu(111), DZ PES: reaction probability, Fortran vs C++ CPU vs C++ GPU", loc="left")

    # ---- bottom: differences
    axb.fill_between(E, 1e-17, half, color="#eeede9", step="mid", lw=0)
    axb.step(E, half, where="mid", color=MUTED, lw=1.0, label="½ unit of the Fortran's 7th printed digit (print resolution; shaded below)")
    axb.axhline(1e-10, color=INK_2, ls=(0, (6, 3)), lw=1.0)
    axb.text(0.402, 1.35e-10, "full-precision criterion 1e-10", fontsize=8, color=INK_2)
    for key, d in (("gpu-cpu", d_gc), ("cpu-fortran", d_cf), ("gpu-fortran", d_gf), ("gpu_strict-fortran", d_sf)):
        st = PAIR[key]
        est = "direct (both on-the-fly tables)" if key == "gpu-cpu" else "propagation only: same C++ analysis of both BK files"
        axb.plot(E, np.maximum(d, 1e-17), color=st["color"], ls=st["ls"], marker=st["marker"], ms=3.5, lw=1.1,
                 label=f"{st['label']}, full precision ({est})")
    # printed comparison: both C++ runs vs the Fortran's printed digits
    for k, P, mk in (("cpu", Pc, "x"), ("gpu", Pg, "+")):
        eq = np.array([f"{c:.6e}" == f"{f:.6e}" for c, f in zip(P, Pf)])
        axb.plot(E[~eq], np.abs(P - Pf)[~eq], ls="none", marker=mk, ms=7, color="#0b0b0b",
                 label=f"{PROGRAM[k]['short']} − Fortran printed value where the 7 digits differ ({(~eq).sum()}; "
                       f"outside ½ unit, explained: item 14)")
        if k == "cpu":
            n_eq, eq_cpu = eq.sum(), eq
    axb.plot(E[eq_cpu], np.full(n_eq, 2.2e-16), ls="none", marker="|", ms=6, color=MUTED,
             label=f"bottom edge: energies where the C++ value has all 7 printed Fortran digits ({n_eq} of {E.size}, CPU)")
    axb.set_yscale("log")
    axb.set_ylim(1.5e-16, 3e-6)
    axb.set_xlabel("Collision (normal incidence) energy E_coll (eV)")
    axb.set_ylabel("|ΔP_reaction|")
    axb.legend(loc="upper left", bbox_to_anchor=(-0.12, -0.2), fontsize=8, ncol=1, borderaxespad=0.0)

    # ---- reproduction box
    c = p0["git_last_src_commit"]
    lines = ["To reproduce", "",
             "PES     Dai–Zhang LEPS", "        (H2Cu111DaiLight)",
             "state   H2 v=0, j=0, mJ=0",
             "        normal incidence",
             f"Z       {js['nZ']} pts, dZ {js['dZ']}",
             f"Zsp     {js['nZSpecular']} pts (specular)",
             f"r       {js['nR']} pts, dr {js['dR']}",
             f"x, y    {js['nX']} x {js['nY']} (60° cell)",
             f"jmax    {js['JMax']} (435 rot. fns)",
             f"dt      {js['dtStep']} a.u., T {js['tEnd']:.0f}",
             f"packet  {js['EMinWP']}–{js['EMaxWP']} eV",
             f"        {js['NrOfEs']} energies",
             "input   inputs/",
             "        thesis_fig5_3_dz/",
             f"commit  {c} (src)"]
    axt.text(0.04, 0.985, "\n".join(lines), family="DejaVu Sans Mono", fontsize=7.9, va="top", ha="left",
             transform=axt.transAxes, bbox=dict(boxstyle="round,pad=0.5", fc="#f6f6f4", ec=MUTED, lw=0.6))
    save(fig, FIGURES / "fig1_reaction_probability.png")
    print("written fig1; max full-precision |dP_reaction|", full)


if __name__ == "__main__":
    main()
