#!/usr/bin/env python3
"""Item 19: physics checks, each evaluated for all three runs.

(a) channels open exactly at their thresholds: a channel is open (S listed / S != 0) iff
    E_total - E_channel > 0; the open sets of the codes are equal; the smallest margin
    |E_total - E_channel| over all decisions is reported.
(b) 0 <= P and sum P <= 1.
(c) symmetry, H2(v=0, j=0) at normal incidence on a C6v surface: P(v, j, mJ) = P(v, j, -mJ), and
    P(v, j, n, m) summed over mJ is equal on each C6v orbit of diffraction channels (orbits from
    the 60-degree hexagonal lattice; check: the diffraction order is constant on every orbit).
    Violations are quantified (relative spread on the orbit where max P > 1e-8).
(d) sub-threshold reaction floor (P_reaction at 0.40-0.50 eV, where the thesis curve is ~0).
(e) norm left on the grids at T.
Writes data/items/item19.json."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from common import CHAEV, RUN_DIRS, WORK, FortranBk, NativeBk, clean_log_lines, parse_fortran_real, save_item
from parse_logs import energy_blocks, fort25

KEYS = ["E_coll", "v", "j", "mJ", "n", "m"]


def c6v_orbits(nm_pairs):
    a1 = np.array([1.0, 0.0])
    a2 = np.array([math.cos(math.pi / 3), math.sin(math.pi / 3)])
    A = np.array([a1, a2])                      # rows: real-space basis
    B = np.linalg.inv(A).T                      # rows: reciprocal basis (a_i . b_j = delta_ij)
    ops = []
    for k in range(6):
        c, s = math.cos(k * math.pi / 3), math.sin(k * math.pi / 3)
        R = np.array([[c, -s], [s, c]])
        ops.append(R)
        ops.append(R @ np.array([[1.0, 0.0], [0.0, -1.0]]))
    Binv = np.linalg.inv(B)
    orbit_of = {}
    for n, m in nm_pairs:
        G = n * B[0] + m * B[1]
        images = []
        for R in ops:
            c = (R @ G) @ Binv
            ci = np.round(c)
            assert np.max(np.abs(c - ci)) < 1e-9
            images.append((int(ci[0]), int(ci[1])))
        orbit_of[(n, m)] = min(images)
    return orbit_of


def per_operation(s):
    """Relative difference P(g G) vs P(G), summed over mJ, for each of the 12 C6v operations."""
    g = s.groupby(["E_coll", "v", "j", "n", "m"])["P"].sum().to_dict()
    a1 = np.array([1.0, 0.0])
    a2 = np.array([math.cos(math.pi / 3), math.sin(math.pi / 3)])
    B = np.linalg.inv(np.array([a1, a2])).T
    Binv = np.linalg.inv(B)
    out = {}
    for k in range(6):
        for mir in (0, 1):
            c, si = math.cos(k * math.pi / 3), math.sin(k * math.pi / 3)
            R = np.array([[c, -si], [si, c]])
            if mir:
                R = R @ np.array([[1.0, 0.0], [0.0, -1.0]])
            rel = []
            for (E, v, j, n, m), p in g.items():
                cc = ((n * B[0] + m * B[1]) @ R.T) @ Binv
                n2, m2 = int(round(cc[0])), int(round(cc[1]))
                q = g.get((E, v, j, n2, m2))
                if q is not None and max(p, q) > 1e-8:
                    rel.append(abs(p - q) / max(p, q))
            rel = np.array(rel)
            out[f"rot{60 * k}{'_mirror' if mir else ''}"] = {"n": int(rel.size), "median": float(np.median(rel)),
                                                              "p95": float(np.percentile(rel, 95))}
    return out


def symmetry(s):
    """s: DataFrame with E_coll, v, j, mJ, n, m, P (and order if present)."""
    out = {}
    pairs = sorted(set(zip(s["n"], s["m"])))
    orbit = c6v_orbits(pairs)
    s = s.assign(orbit=[orbit[(n, m)] for n, m in zip(s["n"], s["m"])])
    if "order" in s:
        out["order_constant_on_orbits"] = bool((s.groupby("orbit")["order"].nunique() == 1).all())
    out["n_orbits"] = len(set(orbit.values()))
    # mJ <-> -mJ (summed over diffraction)
    pv = s.groupby(["E_coll", "v", "j", "mJ"])["P"].sum()
    pm = pv.copy()
    pm.index = pm.index.set_levels(pm.index.levels[3], level=3)
    neg = pv.rename(index=lambda x: x)
    d = []
    for (E, v, j, mj), p in pv.items():
        q = pv.get((E, v, j, -mj), 0.0)
        mx = max(p, q)
        if mx > 1e-8 and mj > 0:
            d.append((abs(p - q) / mx, abs(p - q), E, v, j, mj))
    d.sort(reverse=True)
    out["mJ_symmetry"] = {"n_pairs": len(d), "max_relative": d[0][0] if d else 0.0, "max_abs": max(x[1] for x in d) if d else 0.0,
                          "where": {"E": d[0][2], "v": int(d[0][3]), "j": int(d[0][4]), "mJ": int(d[0][5])} if d else None}
    # the diagonal mirror (n, m) -> (m, n), summed over mJ: the one operation that also maps the
    # parallelogram of FFT momenta in (x, y) (n, m = -nx/2 .. nx/2 - 1) onto itself
    pd_ = s.groupby(["E_coll", "v", "j", "n", "m"])["P"].sum()
    dm = []
    for (E, v, j, n, m), p in pd_.items():
        if n < m:
            q = pd_.get((E, v, j, m, n), np.nan)
            mx = max(p, q)
            if mx > 1e-8:
                dm.append((abs(p - q) / mx, E, v, j, n, m))
    dm.sort(reverse=True)
    out["diagonal_mirror"] = {"n_pairs": len(dm), "max_relative": float(dm[0][0]),
                              "where": {"E": float(dm[0][1]), "v": int(dm[0][2]), "j": int(dm[0][3]), "n": int(dm[0][4]), "m": int(dm[0][5])}}
    # orbits of diffraction channels, summed over mJ
    g = s.groupby(["E_coll", "v", "j", "n", "m", "orbit"])["P"].sum().reset_index()
    spread = []
    for key, grp in g.groupby(["E_coll", "v", "j", "orbit"]):
        mx = grp["P"].max()
        if mx > 1e-8 and len(grp) > 1:
            spread.append(((grp["P"].max() - grp["P"].min()) / mx, grp["P"].max() - grp["P"].min(), key, len(grp)))
    spread.sort(key=lambda x: -x[0])
    rel = np.array([x[0] for x in spread])
    out["diffraction_orbits"] = {"n_orbits_compared": len(spread), "max_relative_spread": float(rel.max()),
                                 "median_relative_spread": float(np.median(rel)), "p95_relative_spread": float(np.percentile(rel, 95)),
                                 "max_abs_spread": float(max(x[1] for x in spread)),
                                 "median_abs_spread": float(np.median([x[1] for x in spread])),
                                 "where": {"E": float(spread[0][2][0]), "v": int(spread[0][2][1]), "j": int(spread[0][2][2]),
                                           "orbit": list(spread[0][2][3]), "members": spread[0][3]}}
    return out


def main():
    res = {}
    E = pd.read_csv(RUN_DIRS["cpu"] / "scattering_summary.tsv", sep="\t", comment="#")["E_coll"].to_numpy()
    tables = {k: pd.read_csv(RUN_DIRS[k] / "s_matrix.tsv", sep="\t", comment="#") for k in ("cpu", "gpu")}
    F = pd.DataFrame(list(fort25(RUN_DIRS["fortran"] / "fort.25")), columns=["ie", "v", "j", "mJ", "order", "n", "m", "ReS", "ImS", "P", "phase"])
    F["E_coll"] = E[F["ie"].to_numpy()]
    F["P"] = F["P"].map(parse_fortran_real)
    F["S0"] = (F["ReS"].map(parse_fortran_real) == 0) & (F["ImS"].map(parse_fortran_real) == 0)
    tables["fortran"] = F[~F["S0"]][KEYS + ["order", "P"]].reset_index(drop=True)

    # (a) thresholds: C++ open = listed rows; Fortran open = rows with S != 0
    nb = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    lab = nb.labels
    e_init_c = float([l.split() for l in (RUN_DIRS["cpu"] / "gas_phase_states.dat").read_text().splitlines() if not l.startswith("#")][0][2])
    open_sets = {}
    margins = {}
    for k, ech, einit in (("cpp", nb.energies, e_init_c), ("fortran", fb.enymo.reshape(fb.nviba, fb.nstata), fb.e_initial)):
        total = E / CHAEV + einit                      # hartree (both codes: E_coll/chaev + eviroi)
        diff = total[:, None, None] - ech[None, :, :]  # (energy, v, is)
        open_sets[k] = diff > 0
        margins[k] = float(np.min(np.abs(diff)))
    idx = {(int(lab["J"][i]), int(lab["mJ"][i]), int(lab["n"][i]), int(lab["m"][i])): i for i in range(nb.nstata)}
    ie_of = {float(e): i for i, e in enumerate(E)}
    listed = {}
    for k, t in tables.items():
        arr = np.zeros_like(open_sets["cpp"])
        for e, v, j, mj, n, m in zip(t["E_coll"], t["v"], t["j"], t["mJ"], t["n"], t["m"]):
            arr[ie_of[float(e)], int(v), idx[(int(j), int(mj), int(n), int(m))]] = True
        listed[k] = arr
    res["thresholds"] = {
        "n_open_cpp_by_energy_rule": int(open_sets["cpp"].sum()), "n_open_fortran_by_energy_rule": int(open_sets["fortran"].sum()),
        "rule_sets_equal": bool(np.array_equal(open_sets["cpp"], open_sets["fortran"])),
        "cpu_listed_equals_rule": bool(np.array_equal(listed["cpu"], open_sets["cpp"])),
        "gpu_listed_equals_rule": bool(np.array_equal(listed["gpu"], open_sets["cpp"])),
        "fortran_nonzero_S_subset_of_rule": bool(np.all(~listed["fortran"] | open_sets["fortran"])),
        "fortran_open_with_S_exactly_zero": int(np.sum(open_sets["fortran"] & ~listed["fortran"])),
        "smallest_margin_hartree": margins, "smallest_margin_eV": {k: v * CHAEV for k, v in margins.items()}}

    # (b) bounds
    blocks = energy_blocks(clean_log_lines("fortran"))
    bounds = {}
    for k in ("cpu", "gpu"):
        su = pd.read_csv(RUN_DIRS[k] / "scattering_summary.tsv", sep="\t", comment="#")
        bounds[k] = {"min_P": float(tables[k]["P"].min()), "max_sum_P": float(su["P_scattered"].max()),
                     "min_P_reaction": float(su["P_reaction"].min())}
    bounds["fortran"] = {"min_P": float(tables["fortran"]["P"].min()), "max_sum_P": max(parse_fortran_real(b["p_sum"]) for b in blocks),
                         "min_P_reaction": min(parse_fortran_real(b["p_reaction"]) for b in blocks)}
    res["bounds"] = bounds

    # (c) symmetry
    res["symmetry"] = {k: symmetry(t) for k, t in tables.items()}
    ops = per_operation(tables["cpu"])
    res["per_operation_cpu"] = ops
    broken = {k: v for k, v in ops.items() if k not in ("rot0", "rot60_mirror")}
    res["per_operation_summary"] = {"broken_median_min": min(v["median"] for v in broken.values()),
                                    "broken_median_max": max(v["median"] for v in broken.values()),
                                    "broken_p95_min": min(v["p95"] for v in broken.values()),
                                    "broken_p95_max": max(v["p95"] for v in broken.values()),
                                    "diagonal_mirror_median": ops["rot60_mirror"]["median"],
                                    "diagonal_mirror_p95": ops["rot60_mirror"]["p95"]}
    # the orbit spread is a property of the calculation: compare it between the codes
    sp = {k: v["diffraction_orbits"]["max_relative_spread"] for k, v in res["symmetry"].items()}
    res["symmetry_spread_codes_max_difference"] = float(max(sp.values()) - min(sp.values()))
    # first-order shell at two energies (v = 0, j = 0, summed over mJ) for the report
    shell = {}
    t = tables["cpu"]
    for E_ in (0.6, 0.865):
        g = t[np.isclose(t["E_coll"], E_) & (t["v"] == 0) & (t["j"] == 0) & (t["order"] == 1)].groupby(["n", "m"])["P"].sum()
        shell[str(E_)] = {f"({n},{m})": float(p) for (n, m), p in g.items()}
    res["first_order_shell_v0_j0"] = shell

    # (d) sub-threshold reaction floor and (e) norm left
    su = pd.read_csv(RUN_DIRS["cpu"] / "scattering_summary.tsv", sep="\t", comment="#")
    low = su[su["E_coll"] <= 0.5 + 1e-9]
    res["subthreshold"] = {"E_max": 0.5, "n": int(len(low)), "P_reaction_min": float(low["P_reaction"].min()),
                           "P_reaction_max": float(low["P_reaction"].max()), "P_reaction_median": float(low["P_reaction"].median()),
                           "fortran_same_range": [min(parse_fortran_real(b["p_reaction"]) for b, e in zip(blocks, E) if e <= 0.5 + 1e-9),
                                                  max(parse_fortran_real(b["p_reaction"]) for b, e in zip(blocks, E) if e <= 0.5 + 1e-9)]}
    av = [l.split() for l in (RUN_DIRS["cpu"] / "averages.dat").read_text().splitlines() if not l.startswith("#")]
    res["norm_left_at_T"] = {"cpu": float(av[-1][1]),
                             "gpu": float([l.split() for l in (RUN_DIRS["gpu"] / "averages.dat").read_text().splitlines() if not l.startswith("#")][-1][1]),
                             "fortran_text": [l.split() for l in open(RUN_DIRS["fortran"] / "averages.out") if len(l.split()) == 5][-1][1]}
    save_item("item19", res)
    print("thresholds:", res["thresholds"])
    print("bounds:", res["bounds"])
    for k, v in res["symmetry"].items():
        print("symmetry", k, v)
    print("subthreshold:", res["subthreshold"])
    print("norm left:", res["norm_left_at_T"])


if __name__ == "__main__":
    main()
