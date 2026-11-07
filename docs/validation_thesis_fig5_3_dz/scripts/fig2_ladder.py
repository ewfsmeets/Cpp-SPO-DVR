#!/usr/bin/env python3
"""Figure 2: error ladder 1, absolute differences in probability per collision energy, one row per
rung on a shared log axis: the CPU<->GPU rounding floor, C++<->Fortran at full precision (replay),
the Fortran's print resolution, and visible artefacts of the calculation that are identical in all
codes. Bars: min-max; dots: the individual values; gaps between adjacent rungs annotated with the
ratio of their medians. Writes data/items/ladder.json (the numbers on the figure)."""
from __future__ import annotations

import numpy as np

from common import DATA, FIGURES, load_item, save_item
from style import PAIR, INK_2, MUTED, plt, save, setup


def col(path, name, skip_comment=1):
    lines = path.read_text().splitlines()
    head = lines[skip_comment].split("\t")
    i = head.index(name)
    return np.array([float(l.split("\t")[i]) for l in lines[skip_comment + 1:]])


def main():
    setup()
    rp, rr = DATA / "reaction_probability.tsv", DATA / "replay_reaction_probability.tsv"
    E = col(rp, "E_coll")
    i19 = load_item("item19")
    i13 = load_item("item13")
    rungs = [
        ("C++ GPU vs C++ CPU,\nfull precision (same C++\ncode, two executions)", PAIR["gpu-cpu"]["color"],
         col(rp, "abs_gpu_minus_cpu"), col(rp, "max_dPvjm_gpu_cpu")),
        ("C++ CPU vs Fortran,\nfull precision (replay)", PAIR["cpu-fortran"]["color"],
         col(rr, "abs_cpu_minus_F"), col(rr, "max_dPvjm_cpu_F")),
        ("Fortran print resolution:\n½ unit of the 7th digit\nof P_reaction", MUTED, col(rp, "half_unit"), None),
    ]
    low = E <= 0.5 + 1e-9
    sym = i19["symmetry"]["cpu"]["diffraction_orbits"]
    artefacts = [
        ("Norm left on the grids\nat T = 17 500 a.u.", np.array([i13["averages"]["norm_left_at_T"]])),
        ("Reaction floor below\nthreshold (E ≤ 0.50 eV)", col(rp, "P_cpu")[low]),
        ("Spread between symmetry-\nequivalent diffraction channels\n(median, max)", np.array([sym["median_abs_spread"], sym["max_abs_spread"]])),
    ]
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    y = 0
    labels, ys = [], []
    medians = []
    rng = np.random.default_rng(1)
    for name, color, a, b in rungs:
        vals = a if b is None else np.concatenate([a, b])
        ax.plot([vals.min(), vals.max()], [y, y], color=color, lw=7, alpha=0.25, solid_capstyle="butt")
        ax.scatter(a, y + rng.uniform(-0.13, 0.13, a.size), s=10, color=color, zorder=3,
                   label=None)
        if b is not None:
            ax.scatter(b, y + rng.uniform(-0.13, 0.13, b.size), s=10, facecolors="white", edgecolors=color, zorder=3)
        medians.append(float(np.median(a)))
        labels.append(name)
        ys.append(y)
        y += 1
    for name, a in artefacts:
        ax.plot([a.min(), a.max()], [y, y], color="#e34948", lw=7, alpha=0.25, solid_capstyle="butt")
        ax.scatter(a, np.full(a.size, y), s=14, color="#e34948", zorder=3)
        medians.append(float(np.median(a)))
        labels.append(name)
        ys.append(y)
        y += 1
    ax.set_yticks(ys, labels, fontsize=8.6)
    ax.set_xscale("log")
    ax.set_xlim(1e-14, 3e-2)
    ax.set_ylim(-0.6, y - 0.4)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Absolute difference or size in probability (per collision energy)")
    ax.axhspan(2.5, y - 0.4, color="#fbeceb", lw=0, zorder=0)
    ax.text(1.5e-14, y - 0.55, "artefacts of this calculation, identical in all three codes", ha="left", va="top",
            fontsize=8.5, color=INK_2)
    # gaps between rungs 0-1, 1-2, 2-3 (medians)
    gaps = []
    for i in range(3):
        r = medians[i + 1] / medians[i]
        gaps.append(r)
        xm = np.sqrt(medians[i] * medians[i + 1])
        ax.annotate("", xy=(medians[i + 1], ys[i] + 0.5), xytext=(medians[i], ys[i] + 0.5),
                    arrowprops=dict(arrowstyle="<->", color=INK_2, lw=0.8))
        sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
        m_, e_ = f"{r:.0e}".split("e")
        label = f"×{r:.0f}" if r < 1000 else f"×{m_}×10{str(int(e_)).translate(sup)}"
        ax.text(xm, ys[i] + 0.56, label, ha="center", va="bottom", fontsize=8.5, color=INK_2)
    ax.scatter([], [], s=10, color=MUTED, label="P_reaction (one dot per energy)")
    ax.scatter([], [], s=10, facecolors="white", edgecolors=MUTED, label="largest |ΔP(v, j, mJ)| at that energy")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Error ladder: how the code differences compare with the scales of the calculation", loc="left")
    save(fig, FIGURES / "fig2_error_ladder.png")
    margin = col(rp, "half_unit") / col(rr, "abs_cpu_minus_F")
    margin_g = col(rp, "half_unit") / col(rr, "abs_gpu_minus_F")
    save_item("ladder", {"print_over_propagation_min": float(margin.min()), "print_over_propagation_min_gpu": float(margin_g.min()), "print_over_propagation_min_at_E": float(E[int(np.argmin(margin))]),
                         "print_over_propagation_E_le_053_max": float(margin[E <= 0.53 + 1e-9].max()),"medians": dict(zip([l.replace("\n", " ") for l in labels], medians)),
                         "gap_ratios_of_medians": gaps})
    print("medians", [f"{m:.2e}" for m in medians], "gaps", [f"{g:.1f}" for g in gaps])


if __name__ == "__main__":
    main()
