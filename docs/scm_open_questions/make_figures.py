"""Figures of docs/scm_open_questions.md from the example SCM data (potentials/H2Cu111SRP_C3v.dir).

A short Python transcription of the 5.24 SCM pair term (DistortedSurfacePES.dir: SCM_Potential.f90
SCM_SurfInteract_3Body and Switch_Rydberg_Potential, SCM_Types.f90 Setup_3Body_Rydberg_Par,
SCM_Wrapper.f90 ConvertCoMToXYZ), accurate to well below 1 meV; for illustration only, not the
implementation. Writes pair_function.png and edge_jump.png next to this script.

Usage: python3 -I -B docs/scm_open_questions/make_figures.py   (run time ≈ 20 s; needs matplotlib)
"""
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "potentials" / "H2Cu111SRP_C3v.dir"
HARTREE_TO_EV = 1 / 0.036749308286  # the precise factor of GetCRPDistSurfPotential
BOHR_IN_ANGSTROM = 0.529177249
A0 = 4.916154  # nearest-neighbour distance of the CRP set H2Cu111SRP_C3v, bohr (2.601516559 Å)


def read_parameters():
    """The 18 parameters of layer 1 (params.dat VERSION 2; both layers are identical)."""
    return [float(x) for x in (DATA / "params.dat").read_text().split()[2:20]]


def read_atoms():
    """(displaced, ideal, layer) of every atom of the example slab, and its lattice constant a_T."""
    lines = (DATA / "examplePosfile.dat").read_text().splitlines()
    a_t = float(lines[2])
    atoms = []
    for line in lines[3:]:
        f = line.split()
        if f:
            atoms.append((tuple(map(float, f[0:3])), tuple(map(float, f[3:6])), int(f[6])))
    return atoms, a_t


def pair_parameters(P, t):
    """p1..p8 at clamped bond length t, and the constant C beyond p8."""
    p = [P[8 + k] * t + P[k] for k in range(8)]
    q8 = p[7] - p[0]
    c = -math.exp(-p[1] * q8) * (p[2] + q8 * (p[3] + q8 * (p[4] + q8 * p[5])))
    return p, c


def pair_function(R, p, c):
    """f(R) in hartree: modified Rydberg, switched to C between p7 and p8."""
    if R > p[7]:
        return c
    q = R - p[0]
    v = -(p[2] + q * (p[3] + q * (p[4] + q * p[5]))) * math.exp(-p[1] * q)
    if R > p[6]:
        rho = 0.5 * math.cos((R - p[7]) * math.pi / (p[7] - p[6])) + 0.5
        v = (1 - rho) * v + rho * c
    return v


def delta_v(P, atoms, s, X, Y, Z, r, theta, phi):
    """ΔV in eV at a Cartesian centre (X, Y, Z), as coded (bohr clamp, rescaling on)."""
    p, c = pair_parameters(P, min(max(r, P[16]), P[17]))
    d = (0.5 * r * math.sin(theta) * math.cos(phi), 0.5 * r * math.sin(theta) * math.sin(phi),
         0.5 * r * math.cos(theta))
    h = [(X - d[0], Y - d[1], Z - d[2]), (X + d[0], Y + d[1], Z + d[2])]
    centre = [(X / s), (Y / s), Z]
    h_ideal = [tuple(centre[k] - d[k] for k in range(3)), tuple(centre[k] + d[k] for k in range(3))]
    v = 0.0
    for displaced, ideal, _ in atoms:
        for a in range(2):
            v = v - pair_function(math.dist(h_ideal[a], ideal), p, c) + pair_function(math.dist(h[a], displaced), p, c)
    return v * HARTREE_TO_EV


def plot_pair_function(P):
    """f(R) for the bohr reading of the clamp at r = 1.4 and r ≥ 1.847 bohr, and the Å reading at r = 3.49 bohr."""
    fig, ax = plt.subplots(figsize=(8.8, 4.4), dpi=100)
    radii = [1.5 + 0.01 * i for i in range(951)]
    curves = [(1.4, "t = 1.4 bohr: H₂ at its bond length (both readings)"),
              (P[17], "t = 1.847: every r ≥ 1.847 bohr, as coded (clamp in bohr)"),
              (P[17] / BOHR_IN_ANGSTROM, "t = 3.49: r = 3.49 bohr if the clamp is in Å")]
    for t, label in curves:
        p, c = pair_parameters(P, t)
        ax.plot(radii, [pair_function(R, p, c) * HARTREE_TO_EV for R in radii], label=label)
    p, _ = pair_parameters(P, 1.4)
    for x, text in ((p[6], "p7"), (p[7], "p8")):
        ax.axvline(x, color="0.6", lw=0.8, ls="--")
        ax.text(x + 0.08, 0.9, text, color="0.4")
    ax.axvline(3.27, color="0.85", lw=4, zorder=0)
    ax.text(3.35, -0.45, "R ≈ 3.3 bohr: H–Cu distance of H in the fcc hollow", color="0.4", fontsize=9)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlim(1.5, 11)
    ax.set_ylim(-0.5, 2.0)
    ax.set_xlabel("H–Cu distance R (bohr)")
    ax.set_ylabel("pair function f(R) (eV)")
    ax.set_title("SCM pair function of the example params.dat for three values of the clamped bond length t")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "pair_function.png", facecolor="white")


def plot_edge_jump(P, atoms, a_t):
    """ΔV along a line across the cell edge v = 0 at u = (2/3)·a, where the jump at this height is largest (0.68 eV),
    as the periodic grid sees it and as the snapshot is."""
    s = a_t / A0
    kept = [a for a in atoms if a[2] <= 2 and sum(x * x for x in a[0]) < 16.0 ** 2]
    u = 2.0 / 3.0 * a_t
    Z, r, theta = 2.25, 1.95, math.pi / 2
    fractions = [-0.5 + 0.01 * i for i in range(101)]

    def at(v):
        return delta_v(P, kept, s, u + 0.5 * v, math.sqrt(3) / 2 * v, Z, r, theta, 0.0)

    periodic = [at((f % 1.0) * a_t) for f in fractions]   # what the periodic x, y grid repeats
    snapshot = [at(f * a_t) for f in fractions]            # the snapshot itself, around the edge
    fig, ax = plt.subplots(figsize=(8.8, 4.4), dpi=100)
    ax.plot(fractions, snapshot, label="the snapshot (atoms around the actual position)", color="0.55", ls="--")
    ax.plot(fractions, periodic, label="what the run uses: the origin cell, repeated periodically", lw=2)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.text(0.01, max(periodic) - 0.02, "cell edge v = 0", color="0.3", fontsize=9)
    ax.set_xlabel("v / a (lattice coordinate across the edge, at u = (2/3)·a)")
    ax.set_ylabel("ΔV (eV)")
    ax.set_title("ΔV across a cell edge near the barrier: Z = 2.25, r = 1.95 bohr, θ = 90°, φ = 0")
    ax.legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / "edge_jump.png", facecolor="white")
    jump = abs(periodic[50] - at(a_t))
    print(f"edge_jump.png: jump at v = 0 {jump:.3f} eV")


def main():
    P = read_parameters()
    atoms, a_t = read_atoms()
    plot_pair_function(P)
    plot_edge_jump(P, atoms, a_t)


if __name__ == "__main__":
    main()
