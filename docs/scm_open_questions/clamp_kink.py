"""Which bond-length clamp produced the published SCM potential? (docs/scm_open_questions.md, Q1)

Clamping r freezes the pair parameters, so ΔV(r) has a kink (a jump of ∂ΔV/∂r) at each clamp limit:
at 1.847 bohr with the clamp as coded, at 1.370 and 3.490 bohr (0.725 and 1.847 Å) with the published
one. This script finds the kinks in

  - the published elbow cuts of Smits & Somers, J. Chem. Phys. 158, 014704 (2023), data folder FIG23
    (V_SCM of three production slabs and V_BOSS on the same (r, Z) grid; ΔV = V_SCM − V_BOSS), and
  - ΔV of the example slab on the same grid, with each clamp (the controls).

A kink between grid points r_i and r_i+1 shows as a spike of the second difference of ΔV along r; the
score of a point sums |d2(r) − (d2(r − 2h) + d2(r + 2h))/2| / max|d2| over the rows Z ≤ 3 Å.
Writes clamp_kink.png next to this script and prints the strongest spikes of every cut.

Usage: python3 -I -B docs/scm_open_questions/clamp_kink.py DIR   (run time ≈ 2 s; needs matplotlib)
DIR holds the published files Cu111_BOSS_briT90P90.2d.elbow, Cu111_Surf252_briT90P90.2d.elbow, …, from
https://pubs.tc.lic.leidenuniv.nl/public/bauke.smits/2023_JCP_QDEAMSCM_Quantum/FIG23/
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # python3 -I leaves the script's folder off the path
import make_figures as mf  # noqa: E402  (the ΔV transcription and the example data)
import matplotlib.pyplot as plt  # noqa: E402

PUBLISHED = [("Cu111_BOSS_briT90P90", "Cu111_Surf252_briT90P90"),
             ("Cu111_BOSS_briT90P90", "Cu111_Surf556_briT90P90"),
             ("Cu111_BOSS_briT90P90", "Cu111_Surf1012_briT90P90"),
             ("Cu111_BOSS_t2fAT90P30", "Cu111_Surf252_t2fAT90P30")]
R_GRID = [0.4 + 0.019 * i for i in range(101)]    # Å, as in the published files
Z_GRID = [0.25 + 0.0625 * j for j in range(45)]   # Å, the rows up to Z = 3 Å
B = mf.BOHR_IN_ANGSTROM


def read_elbow(path):
    """{(r, Z): V} of a published elbow file (columns r (Å), Z (Å), V (eV))."""
    grid = {}
    for line in Path(path).read_text().splitlines():
        f = line.split()
        if len(f) >= 3:
            grid[(round(float(f[0]), 6), round(float(f[1]), 6))] = float(f[2])
    return grid


def published_delta_v(directory, boss_name, slab_name):
    """ΔV(r, Z) = V_SCM − V_BOSS of one published cut, rows Z ≤ 3 Å, as {Z: [ΔV at each r]}."""
    boss = read_elbow(Path(directory) / f"{boss_name}.2d.elbow")
    slab = read_elbow(Path(directory) / f"{slab_name}.2d.elbow")
    return {z: [slab[(round(r, 6), round(z, 6))] - boss[(round(r, 6), round(z, 6))] for r in R_GRID] for z in Z_GRID}


def example_delta_v(clamp):
    """ΔV of the example slab on the elbow grid at the bridge (1/2, 0)·a_T, θ = φ = 90°."""
    P = mf.read_parameters()
    atoms, a_t = mf.read_atoms()
    kept, s = mf.kept_atoms(atoms), a_t / mf.A0
    return {z: [mf.delta_v(P, kept, s, 0.5 * a_t, 0.0, z / B, r / B, math.pi / 2, math.pi / 2, clamp(P))
                for r in R_GRID] for z in Z_GRID}


def spike_scores(rows):
    """The kink score of every interior r (see the module docstring)."""
    scores = [0.0] * len(R_GRID)
    for dv in rows.values():
        d2 = [dv[i + 1] - 2 * dv[i] + dv[i - 1] for i in range(1, len(dv) - 1)]   # d2[k] belongs to r[k+1]
        scale = max(abs(x) for x in d2) or 1.0
        for k in range(2, len(d2) - 2):
            scores[k + 1] += abs(d2[k] - 0.5 * (d2[k - 2] + d2[k + 2])) / scale
    return scores


def report(label, rows):
    scores = spike_scores(rows)
    top = sorted(range(len(scores)), key=lambda i: -scores[i])[:2]
    median = sorted(scores)[len(scores) // 2]
    spikes = ", ".join(f"r = {R_GRID[i]:.3f} Å ({R_GRID[i] / B:.3f} bohr): {scores[i]:.1f}" for i in top)
    print(f"{label:34s} strongest spikes {spikes}; median {median:.2f}")


def plot(published, coded, clamped):
    """ΔV(r) at three heights: the published slab-252 bridge cut, and the example slab with each clamp."""
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0), dpi=100)
    for z_angstrom, colour in ((1.0, "#1f77b4"), (1.1875, "#d62728"), (1.5, "#2ca02c")):
        r_bohr = [r / B for r in R_GRID]
        label = f"Z = {z_angstrom / B:.2f} bohr"
        axes[0].plot(r_bohr, published[z_angstrom], color=colour, label=label)
        axes[1].plot(r_bohr, coded[z_angstrom], color=colour, label=label + ", as coded")
        axes[1].plot(r_bohr, clamped[z_angstrom], color=colour, ls="--")
    axes[1].plot([], [], color="0.3", ls="--", label="published clamp")
    for ax in axes:
        ax.axvline(1.847, color="0.2", lw=0.9)
        for limit in (0.725 / B, 1.847 / B):
            ax.axvline(limit, color="0.55", lw=0.9, ls=":")
        ax.set_xlim(1.0, 4.35)
        ax.set_xlabel("bond length r (bohr)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("ΔV (eV)")
    axes[1].set_ylabel("ΔV (eV)")
    axes[0].set_title("Published: slab 252, bridge, θ = φ = 90°", fontsize=10)
    axes[1].set_title("Example slab, same site and grid", fontsize=10)
    fig.tight_layout()
    fig.savefig(mf.HERE / "clamp_kink.png", facecolor="white")


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    coded = example_delta_v(mf.coded_clamp)
    clamped = example_delta_v(mf.published_clamp)
    report("example slab, clamp as coded", coded)
    report("example slab, published clamp", clamped)
    published = {}
    for boss_name, slab_name in PUBLISHED:
        rows = published_delta_v(directory, boss_name, slab_name)
        report(slab_name, rows)
        published.setdefault(slab_name, rows)
    plot(published["Cu111_Surf252_briT90P90"], coded, clamped)


if __name__ == "__main__":
    main()
