"""Figure of docs/flux_analysis_findings.md: the flux of the initial channel at r_c.

Reads initial_channel_rc2.0.tsv and initial_channel_rc4.1.tsv (written by specular_flux.cpp) and
writes initial_channel.png next to this script: per flux surface r_c, the flux analysis's P_scat
(scattering grid only, what the Fortran computes), the flux P_spec of the initial channel that it
leaves out, their sum P_total, and the reaction probability 1 - sum R of the S-matrix analysis with
a band of +-(its largest negative value), the noise of 1 - sum R on this small grid.

Usage: python3 -I -B docs/flux_analysis_findings/make_figure.py   (run time < 5 s; needs matplotlib)
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
CASES = (("initial_channel_rc2.0.tsv", "r_c = 2.0 bohr (near the barrier)"),
         ("initial_channel_rc4.1.tsv", "r_c = 4.1 bohr (beyond the barrier)"))
# Categorical slots 1-3 of the reference palette; the reference curve in a neutral ink.
SERIES = (("P_scat", "P_scat (flux analysis, scattering grid only)", "#2a78d6"),
          ("P_spec", "P_spec (initial channel, left out)", "#eb6834"),
          ("P_total", "P_total = P_scat + P_spec", "#1baf7a"))
REFERENCE_INK = "#52514e"
GRID_INK = "#e4e3df"


def read_columns(path):
    """The TSV's columns as lists of floats, keyed by header name."""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {key: [float(r[key]) for r in rows] for key in rows[0]}


def plot_case(ax, data, title):
    energy = data["E_flux_eV"]
    reference = data["one_minus_sumR"]
    noise = max(-x for x in reference)  # largest negative 1 - sum R
    ax.fill_between(energy, [x - noise for x in reference], [x + noise for x in reference],
                    color=GRID_INK, linewidth=0, label=f"1 - ΣR ± {noise:.1e} (its noise)")
    ax.plot(energy, reference, color=REFERENCE_INK, linewidth=1.5, linestyle=(0, (4, 2)),
            label="1 - ΣR (S-matrix analysis)")
    # P_total underneath and wider: at r_c = 4.1 bohr it lies within 9 % of P_scat.
    for key, label, color in SERIES:
        wide = key == "P_total"
        ax.plot(energy, data[key], color=color, linewidth=4 if wide else 2, zorder=2 if wide else 3,
                solid_capstyle="round", label=label)
    ax.axhline(0.0, color=REFERENCE_INK, linewidth=0.8)
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_xlabel("collision energy E (eV)")
    ax.grid(color=GRID_INK, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True, constrained_layout=True)
    for ax, (name, title) in zip(axes, CASES):
        plot_case(ax, read_columns(HERE / name), title)
    axes[0].set_ylabel("probability")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle("cu111_dailight_small + flux analysis, 7000 steps: the flux of the initial channel",
                 fontsize=12, x=0.01, ha="left")
    fig.savefig(HERE / "initial_channel.png", dpi=110)


if __name__ == "__main__":
    main()
