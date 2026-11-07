#!/usr/bin/env python3
"""Negative control (heavy: writes 2 x 266 MB copies of fort.21 and runs the replay helper):
inject a known error into one channel of the Fortran's BK coefficients, as if the propagation had
a channel-local error, replay it through the same C++ analysis, and run the comparison of item 18b
(per-energy uniform factor removed, per-channel residual) against replay(C++ CPU BK file).
Shows the size of a single-channel error the comparison resolves. Writes data/items/negative_control.json."""
from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pandas as pd

from common import INPUT_DIR, RUN_DIRS, WORK, FortranBk, NativeBk, save_item

KEYS = ["E_coll", "v", "j", "mJ", "n", "m"]
CHANNEL = {"v": 0, "J": 0, "mJ": 0, "n": -1, "m": -1}   # first diffraction order, |S|^2 ~ 1e-3
HELPER = WORK / "helpers" / "build-fast" / "replay"


def per_channel_residual(sa_df, sb_df):
    sa = sa_df["Re_S"].to_numpy() + 1j * sa_df["Im_S"].to_numpy()
    sb = sb_df["Re_S"].to_numpy() + 1j * sb_df["Im_S"].to_numpy()
    E = sb_df["E_coll"].to_numpy()
    resid = np.empty(len(sb))
    for e in np.unique(E):
        m = E == e
        z = np.sum(np.conj(sb[m]) * (sa[m] - sb[m])) / np.sum(np.abs(sb[m]) ** 2)
        resid[m] = np.abs((sa[m] - sb[m]) - z * sb[m])
    sel = np.abs(sb) ** 2 > 1e-8
    rel = np.where(sel, resid / np.where(sel, np.abs(sb), 1), np.nan)
    return rel, sel


def main():
    out_dir = WORK / "negative_control"
    out_dir.mkdir(parents=True, exist_ok=True)
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    nb = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    lab = nb.labels
    k = int(np.flatnonzero((lab["J"] == CHANNEL["J"]) & (lab["mJ"] == CHANNEL["mJ"]) & (lab["n"] == CHANNEL["n"])
                           & (lab["m"] == CHANNEL["m"]))[0]) + CHANNEL["v"] * fb.nstata
    base = pd.read_csv(WORK / "replay" / "build-fast" / "cpu" / "s_matrix.tsv", sep="\t", comment="#")
    ref = pd.read_csv(WORK / "replay" / "build-fast" / "fortran" / "s_matrix.tsv", sep="\t", comment="#")
    rel0, sel = per_channel_residual(ref, base)
    is_ch = ((base["v"] == CHANNEL["v"]) & (base["j"] == CHANNEL["J"]) & (base["mJ"] == CHANNEL["mJ"])
             & (base["n"] == CHANNEL["n"]) & (base["m"] == CHANNEL["m"])).to_numpy()
    res = {"channel": CHANNEL, "unperturbed": {"channel_residual_median": float(np.nanmedian(rel0[is_ch & sel])),
                                              "all_channels_p95": float(np.nanpercentile(rel0[sel], 95)),
                                              "all_channels_p99": float(np.nanpercentile(rel0[sel], 99))}}
    for delta in (1e-11, 1e-10, 1e-9):
        path = out_dir / f"fort.21.delta{delta:.0e}"
        shutil.copyfile(RUN_DIRS["fortran"] / "fort.21", path)
        mm = np.memmap(path, dtype=np.uint8, mode="r+")
        for r in range(fb.n_records):
            o = fb.first_offset + r * fb.stride + 16 * k
            c = np.frombuffer(mm[o:o + 16].tobytes(), dtype="<c16")[0] * (1.0 + delta)
            mm[o:o + 16] = np.frombuffer(np.array([c], dtype="<c16").tobytes(), dtype=np.uint8)
        mm.flush()
        del mm
        rdir = out_dir / f"replay_delta{delta:.0e}"
        subprocess.run([str(HELPER), str(INPUT_DIR / "input.json"), str(path), "fortran", str(rdir)], check=True,
                       capture_output=True)
        pert = pd.read_csv(rdir / "s_matrix.tsv", sep="\t", comment="#")
        rel, _ = per_channel_residual(pert, base)
        ch = rel[is_ch & sel]
        p95 = float(np.nanpercentile(rel[sel], 95))
        p_r = pd.read_csv(rdir / "scattering_summary.tsv", sep="\t", comment="#")["P_reaction"].to_numpy()
        p_r0 = pd.read_csv(WORK / "replay" / "build-fast" / "fortran" / "scattering_summary.tsv", sep="\t",
                           comment="#")["P_reaction"].to_numpy()
        res[f"{delta:.0e}"] = {"delta": delta, "channel_residual_median": float(np.nanmedian(ch)),
                               "channel_residual_min": float(np.nanmin(ch)), "all_channels_p95": p95,
                               "channel_median_over_p95": float(np.nanmedian(ch) / p95),
                               "energies_where_channel_above_p99_of_unperturbed": int(np.sum(ch > res["unperturbed"]["all_channels_p99"])),
                               "n_energies": int(ch.size),
                               "max_change_P_reaction": float(np.max(np.abs(p_r - p_r0)))}
        path.unlink()
    save_item("negative_control", res)
    for key, v in res.items():
        print(key, v)


if __name__ == "__main__":
    main()
