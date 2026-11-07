#!/usr/bin/env python3
"""Item 18a: the C++ fine analysis (`onlyAnalysis`) of each run's BK file: fort.21 (Fortran run),
bk_coefficients.bin (C++ CPU and GPU runs). A different estimator from the on-the-fly results
(windowed FFT of the BK series, interpolated to the analysis energies), so the three are compared
only with each other, at full precision (criterion 1e-10 absolute for P, S, P_reaction).
Writes data/items/item18a.json and data/fine_analysis_reaction_probability.tsv."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import DATA, PAIRS, RUNS, WORK, save_item, stats

DIRS = {k: RUNS / f"thesis_fig5_3_dz_oa_{k}" / "fine_analysis" for k in ("cpu", "gpu", "fortran")}
KEYS = ["E_coll", "v", "j", "mJ", "n", "m"]


def load(k, name, sub=""):
    return pd.read_csv(DIRS[k] / sub / name, sep="\t", comment="#")


def main():
    res = {"grid": {k: json.loads((DIRS[k] / "grid.json").read_text()) for k in DIRS}}
    res["grid_identical"] = all(res["grid"][k] == res["grid"]["fortran"] for k in DIRS)
    tables = {k: {"summary": load(k, "scattering_summary.tsv"), "rovib": load(k, "probabilities_rovib.tsv"),
                  "diffraction": load(k, "probabilities_diffraction.tsv"), "s": load(k, "s_matrix.tsv")} for k in DIRS}
    spectrum = {k: {"summary": load(k, "scattering_summary.tsv", "spectrum"),
                    "rovib": load(k, "probabilities_rovib.tsv", "spectrum")} for k in DIRS}
    res["n_energies"] = int(len(tables["fortran"]["summary"]))
    res["n_spectrum_points"] = int(len(spectrum["fortran"]["summary"]))
    for a, b in PAIRS:
        ta, tb = tables[a], tables[b]
        r = {}
        for col in ("P_scattered", "P_reaction", "P_v0", "P_v1"):
            d = np.abs(ta["summary"][col].to_numpy() - tb["summary"][col].to_numpy())
            r[col] = stats(d, where=[float(x) for x in tb["summary"]["E_coll"]])
        assert (ta["rovib"][["E_coll", "v", "j", "mJ"]].to_numpy() == tb["rovib"][["E_coll", "v", "j", "mJ"]].to_numpy()).all()
        r["P_vjm"] = stats(np.abs(ta["rovib"]["P"].to_numpy() - tb["rovib"]["P"].to_numpy()))
        r["P_order"] = stats(np.abs(ta["diffraction"]["P"].to_numpy() - tb["diffraction"]["P"].to_numpy()))
        sa, sb = ta["s"], tb["s"]
        same_rows = len(sa) == len(sb) and (sa[KEYS].to_numpy() == sb[KEYS].to_numpy()).all()
        r["s_rows_equal"] = bool(same_rows)
        if same_rows:
            dS = np.abs((sa["Re_S"] - sb["Re_S"]).to_numpy() + 1j * (sa["Im_S"] - sb["Im_S"]).to_numpy())
            r["S"] = stats(dS)
            r["P"] = stats(np.abs(sa["P"].to_numpy() - sb["P"].to_numpy()))
            sel = sb["P"].to_numpy() > 1e-8
            r["phase_rad_P_gt_1e-8"] = stats(np.abs(np.angle(np.exp(1j * (sa["phase"].to_numpy()[sel] - sb["phase"].to_numpy()[sel])))))
        sp = np.abs(spectrum[a]["summary"]["P_reaction"].to_numpy() - spectrum[b]["summary"]["P_reaction"].to_numpy())
        r["spectrum_P_reaction"] = stats(sp, where=[float(x) for x in spectrum[b]["summary"]["E_coll"]])
        r["within_1e-10"] = bool(max(r["P_reaction"]["max"], r["P_vjm"]["max"], r.get("S", {"max": 0})["max"]) <= 1e-10)
        res[f"{a}-{b}"] = r
    # per-energy table (the fine-analysis P_reaction of the three runs)
    lines = ["# Fine analysis (C++ onlyAnalysis) of each run's BK file: P_reaction at the analysis energies (interpolated).",
             "E_coll\tP_reaction_fortran\tP_reaction_cpu\tP_reaction_gpu"]
    for i in range(res["n_energies"]):
        lines.append(f"{tables['fortran']['summary']['E_coll'][i]:.4f}\t{tables['fortran']['summary']['P_reaction'][i]:.16e}"
                     f"\t{tables['cpu']['summary']['P_reaction'][i]:.16e}\t{tables['gpu']['summary']['P_reaction'][i]:.16e}")
    (DATA / "fine_analysis_reaction_probability.tsv").write_text("\n".join(lines) + "\n")
    save_item("item18a", res)
    print("grid identical:", res["grid_identical"], "energies", res["n_energies"], "spectrum points", res["n_spectrum_points"])
    for a, b in PAIRS:
        r = res[f"{a}-{b}"]
        print(f"{a}-{b}: P_reaction max {r['P_reaction']['max']:.2e} (median {r['P_reaction']['median']:.2e}); P(v,j,mJ) max "
              f"{r['P_vjm']['max']:.2e}; S max {r.get('S', {}).get('max', float('nan')):.2e}; spectrum P_r max "
              f"{r['spectrum_P_reaction']['max']:.2e}; within 1e-10: {r['within_1e-10']}")


if __name__ == "__main__":
    main()
