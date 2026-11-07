#!/usr/bin/env python3
"""Item 17: comparison with Figure 5.3 of M. F. Somers' thesis (dotted curve = DZ PES).

--digitise (extract mode; needs the output of thesis_extract.py): calibrate the axes on the tick
marks of the vector figure, take the dotted curve's dots (isolated small filled shapes inside the
frame: connected components of at most 2.2 pt), and freeze them in
data/thesis_fig5_3_dz_digitised.tsv plus the check image figures/thesis_digitisation_check.png
(the extracted points over the rendered original). Done before any run is drawn on it.

default: compare the runs with the frozen points (criteria.md, "Thesis"):
  agrees if |P_run(E) - P_thesis(E)| <= u_P + |dP/dE| u_E + 0.01 at every digitised point in
  0.40-0.95 eV, with u_P, u_E the largest tick calibration residual plus half the drawn dot size.
Writes data/items/item17.json."""
from __future__ import annotations

import json
import sys

import numpy as np

from common import DATA, FIGURES, WORK, save_item, load_item

THESIS = WORK / "thesis"
DIGITISED = DATA / "thesis_fig5_3_dz_digitised.tsv"


def components(rects, gap=0.25):
    n = len(rects)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            a, b = rects[i], rects[j]
            if a[0] - gap <= b[2] and b[0] - gap <= a[2] and a[1] - gap <= b[3] and b[1] - gap <= a[3]:
                parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def calibrate(centres, lo_edge, hi_edge, lo_val, hi_val, step):
    """Linear map from tick positions (pt) to values; ticks get the nearest multiple of `step`."""
    pos = np.array(list(centres) + [lo_edge, hi_edge])
    scale = (hi_val - lo_val) / (hi_edge - lo_edge)
    val = lo_val + np.round((pos - lo_edge) * scale / step) * step
    b, a = np.polyfit(val, pos, 1)          # pos = a + b * val
    resid = (pos - (a + b * val)) / b       # in data units
    return a, b, val, resid


def digitise():
    J = json.load(open(THESIS / "fig5_3_drawings.json"))
    meta, D = J["meta"], J["drawings"]
    fills = [d for d in D if d["type"] == "f"]
    longs = [d["rect"] for d in fills if max(d["rect"][2] - d["rect"][0], d["rect"][3] - d["rect"][1]) > 250]
    # frame: two horizontal and two vertical long thin rectangles (the page rule at the top is wider)
    horiz = sorted([r for r in longs if (r[2] - r[0]) > 250 and (r[2] - r[0]) < 300], key=lambda r: r[1])
    vert = sorted([r for r in longs if (r[3] - r[1]) > 250], key=lambda r: r[0])
    x_left, x_right = (vert[0][0] + vert[0][2]) / 2, (vert[-1][0] + vert[-1][2]) / 2
    y_top, y_bottom = (horiz[0][1] + horiz[0][3]) / 2, (horiz[-1][1] + horiz[-1][3]) / 2
    w = lambda r: r[2] - r[0]
    h = lambda r: r[3] - r[1]
    xt = [((r[0] + r[2]) / 2) for r in (d["rect"] for d in fills) if 0.5 < w(r) < 1.1 and abs(h(r) - 3.12) < 0.2]
    yt = [((r[1] + r[3]) / 2) for r in (d["rect"] for d in fills) if 0.5 < h(r) < 1.1 and abs(w(r) - 3.12) < 0.2]
    ax_, bx, xv, xres = calibrate(xt, x_left, x_right, 0.40, 1.00, 0.01)
    ay, by, yv, yres = calibrate(yt, y_bottom, y_top, 0.0, 0.8, 0.05)
    # dots of the dotted curve: small connected components strictly inside the frame
    inside = [d["rect"] for d in fills if d["rect"][0] > x_left + 0.4 and d["rect"][2] < x_right - 0.4
              and d["rect"][1] > y_top + 0.4 and d["rect"][3] < y_bottom - 0.4]
    groups = components(inside)
    dots = []
    for g in groups:
        x0 = min(inside[i][0] for i in g); y0 = min(inside[i][1] for i in g)
        x1 = max(inside[i][2] for i in g); y1 = max(inside[i][3] for i in g)
        ext = max(x1 - x0, y1 - y0)
        if ext <= 2.2:
            dots.append(((x0 + x1) / 2, (y0 + y1) / 2, ext))
    dots.sort()
    E = [(x - ax_) / bx for x, _, _ in dots]
    P = [(y - ay) / by for _, y, _ in dots]
    size = float(np.median([e for _, _, e in dots]))
    u_E = float(np.max(np.abs(xres)) + 0.5 * size / abs(bx))
    u_P = float(np.max(np.abs(yres)) + 0.5 * size / abs(by))
    info = {"source": "~/work/mf_somers_phd_thesis.pdf", "printed_page": meta["printed_page"], "pymupdf": meta["pymupdf"],
            "frame_pt": {"x_left": x_left, "x_right": x_right, "y_top": y_top, "y_bottom": y_bottom},
            "n_x_ticks": len(xt), "n_y_ticks": len(yt), "x_calibration_residual_max_eV": float(np.max(np.abs(xres))),
            "y_calibration_residual_max_P": float(np.max(np.abs(yres))), "pt_per_eV": float(bx), "pt_per_P": float(abs(by)),
            "dot_size_pt_median": size, "u_E_eV": u_E, "u_P": u_P, "n_points": len(dots),
            "E_range": [min(E), max(E)], "n_components_inside": len(groups)}
    lines = ["# Somers PhD thesis, Figure 5.3, dotted curve (DZ PES), digitised from the PDF's vector drawing "
             "(scripts/thesis_extract.py, item17_thesis.py --digitise).",
             f"# Axes calibrated on {len(xt)} x and {len(yt)} y tick marks: max residual "
             f"{info['x_calibration_residual_max_eV']:.1e} eV, {info['y_calibration_residual_max_P']:.1e} in P; "
             f"dot size {size:.2f} pt; uncertainty u_E = {u_E:.4f} eV, u_P = {u_P:.4f}.",
             "E_normal_eV\tP_reaction"]
    lines += [f"{e:.5f}\t{p:.5f}" for e, p in zip(E, P)]
    DATA.mkdir(parents=True, exist_ok=True)
    DIGITISED.write_text("\n".join(lines) + "\n")
    save_item("item17_digitisation", info)
    check_image(meta, dots, xt, yt, (x_left, x_right, y_top, y_bottom))
    print(json.dumps(info, indent=1))


def check_image(meta, dots, xt, yt, frame):
    import matplotlib.image as mpimg
    from style import plt, save, setup
    setup()
    img = mpimg.imread(THESIS / "fig5_3_crop.png")
    cx0, cy0 = meta["clip_pt"][0], meta["clip_pt"][1]
    s = meta["dpi"] / 72.0
    fig, ax = plt.subplots(figsize=(6.4, 6.9))
    ax.imshow(img, cmap="gray")
    ax.scatter([(x - cx0) * s for x, _, _ in dots], [(y - cy0) * s for _, y, _ in dots], s=14, facecolors="none",
               edgecolors="#e34948", linewidths=1.0, label=f"extracted points of the dotted curve ({len(dots)})")
    ax.scatter([(x - cx0) * s for x in xt], [(frame[3] + 4 - cy0) * s] * len(xt), marker="|", s=30, color="#2a78d6",
               label=f"x tick marks used ({len(xt)})")
    ax.scatter([(frame[0] - 4 - cx0) * s] * len(yt), [(y - cy0) * s for y in yt], marker="_", s=30, color="#2a78d6",
               label=f"y tick marks used ({len(yt)})")
    ax.set_axis_off()
    ax.legend(loc="lower right", bbox_to_anchor=(0.97, 0.16), fontsize=8, frameon=True, facecolor="white", framealpha=0.95)
    ax.set_title("Digitisation check: Somers thesis Fig. 5.3 (page 170), extracted points in red", fontsize=10)
    save(fig, FIGURES / "thesis_digitisation_check.png")


def compare():
    pts = np.loadtxt(DIGITISED, comments="#", skiprows=3)
    info = load_item("item17_digitisation")
    E_t, P_t = pts[:, 0], pts[:, 1]
    rows = [l.split("\t") for l in (DATA / "reaction_probability.tsv").read_text().splitlines()[2:]]
    E = np.array([float(r[0]) for r in rows])
    runs = {"fortran": np.array([float(r[2]) for r in rows]), "cpu": np.array([float(r[4]) for r in rows]),
            "gpu": np.array([float(r[5]) for r in rows])}
    sel = (E_t >= E[0]) & (E_t <= E[-1])
    res = {"n_points_total": int(len(E_t)), "n_points_in_range": int(sel.sum()), "range_eV": [float(E[0]), float(E[-1])],
           "u_E_eV": info["u_E_eV"], "u_P": info["u_P"], "allowance": 0.01,
           "x_axis": "normal incidence energy (thesis caption) = E_coll of the runs (normal incidence, EParallel = 0)"}
    for k, P in runs.items():
        Pi = np.interp(E_t[sel], E, P)
        slope = np.abs(np.gradient(P, E))
        sl = np.interp(E_t[sel], E, slope)
        tol = info["u_P"] + sl * info["u_E_eV"] + 0.01
        d = Pi - P_t[sel]
        i = int(np.argmax(np.abs(d)))
        res[k] = {"max_abs_dP": float(np.max(np.abs(d))), "at_E": float(E_t[sel][i]), "P_thesis_there": float(P_t[sel][i]),
                  "P_run_there": float(Pi[i]), "mean_dP": float(np.mean(d)), "median_abs_dP": float(np.median(np.abs(d))),
                  "n_within": int(np.sum(np.abs(d) <= tol)), "agrees": bool(np.all(np.abs(d) <= tol)),
                  "max_ratio_to_tolerance": float(np.max(np.abs(d) / tol))}
        # best E shift (diagnostic): shift the thesis points to minimise the rms difference
        shifts = np.linspace(-0.03, 0.03, 601)
        rms = [np.sqrt(np.mean((np.interp(E_t[sel] + s, E, P, left=np.nan, right=np.nan) - P_t[sel]) ** 2)) for s in shifts]
        j = int(np.nanargmin(rms))
        res[k]["best_shift_eV"] = float(shifts[j])
        res[k]["rms_at_best_shift"] = float(rms[j])
        res[k]["rms_unshifted"] = float(np.sqrt(np.mean(d ** 2)))
    save_item("item17", res)
    for k in runs:
        r = res[k]
        print(k, f"max |dP| {r['max_abs_dP']:.4f} at E={r['at_E']:.3f} (thesis {r['P_thesis_there']:.3f}, run {r['P_run_there']:.3f}); "
              f"within {r['n_within']}/{res['n_points_in_range']}; agrees {r['agrees']}; mean dP {r['mean_dP']:+.4f}; "
              f"best shift {r['best_shift_eV']*1000:+.1f} meV (rms {r['rms_at_best_shift']:.4f} vs {r['rms_unshifted']:.4f})")


if __name__ == "__main__":
    digitise() if "--digitise" in sys.argv else compare()
