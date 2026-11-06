#!/usr/bin/env python3
"""Checks the Z grids, wave packet and absorbers of an input.json against the rules in inputs/README.md.

Usage: uv run -q --with numpy python -I tools/wave_packet_check.py <input.json> --zasymp <Z_asymp (bohr)>

Prints the packet widths, the momentum margins, the packet position checks, and the 1D reflection R
and transmission T of the Z, specular-Z and r absorbers (exact quadratic absorber and strength formula
of the code). Units: legacy factors of the code (1 amu = 1822.89 me, 1 Ha = 27.2114 eV).
"""
import argparse
import json
from math import erfc, sqrt

import numpy as np

AMU, HARTREE_EV = 1822.89, 27.2114

# Vibok-Balint-Kurti table (Fortran: MyFunctions.f90 / C++: optical_potential_strength).
XTAB = np.arange(1.0, 15.01, 0.5)
FTAB = np.array([2.750, 2.640, 2.510, 2.380, 2.230, 2.090, 1.940, 1.800, 1.670, 1.550, 1.440, 1.340, 1.250, 1.170,
                 1.100, 1.040, 0.986, 0.949, 0.903, 0.870, 0.840, 0.813, 0.788, 0.764, 0.741, 0.718, 0.697, 0.678,
                 0.663])
ENERGIES_EV = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5]


def lengths_per_wavelength(target_ev, length, mass):
    return length * sqrt(2.0 * mass * target_ev / HARTREE_EV) / (2.0 * np.pi)


def strength(target_ev, length, mass):
    """Absorber strength (Ha); clamped to the table ends like the code."""
    x = min(max(lengths_per_wavelength(target_ev, length, mass), XTAB[0]), XTAB[-1])
    return float(np.interp(x, XTAB, FTAB)) * target_ev / HARTREE_EV


def reflection_transmission(e_ev, start, end, target_ev, mass, h=0.002):
    """1D R and T of V = -i A (3/2) ((q-start)/L)^2 on [start, end] (transfer matrix, step h)."""
    a = strength(target_ev, end - start, mass)
    e = e_ev / HARTREE_EV
    k0 = sqrt(2.0 * mass * e)
    q = np.arange(start, end, h) + h / 2
    v = -1j * a * 1.5 * ((q - start) / (end - start)) ** 2
    psi, dpsi = 1.0 + 0j, 1j * k0  # transmitted wave e^{ik0 q} beyond the end, integrated backwards
    for vi in v[::-1]:
        k = np.sqrt(2.0 * mass * (e - vi) + 0j)
        c, s = np.cos(k * h), np.sin(k * h)
        psi, dpsi = c * psi - s / k * dpsi, k * s * psi + c * dpsi
    incoming = (psi + dpsi / (1j * k0)) / 2
    reflected = (psi - dpsi / (1j * k0)) / 2
    return abs(reflected / incoming) ** 2, 1.0 / abs(incoming) ** 2


def tail(n_sigma):
    """One-sided Gaussian tail beyond n_sigma standard deviations."""
    return 0.5 * erfc(n_sigma / sqrt(2.0))


def absorber_table(name, start, end, target, mass, grid_end):
    ratio = lengths_per_wavelength(target, end - start, mass)
    flag = "" if 1.0 <= ratio <= 15.0 else "  OUTSIDE the table [1, 15]: strength clamped"
    print(f"\n{name}: {start:.3f}-{end:.3f} bohr, target {target} eV, L/lambda = {ratio:.1f}{flag}")
    if abs(end - grid_end) > 1e-9:
        print(f"  warning: absorber end {end:.3f} != last grid point {grid_end:.3f} (waves beyond it wrap around)")
    print("  E (eV)   " + " ".join(f"{e:8.2f}" for e in ENERGIES_EV))
    rt = [reflection_transmission(e, start, end, target, mass) for e in ENERGIES_EV]
    print("  R        " + " ".join(f"{r:8.0e}" for r, _ in rt))
    print("  T        " + " ".join(f"{t:8.0e}" for _, t in rt))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input")
    parser.add_argument("--zasymp", type=float, required=True, help="Z_asymptotic of the PES (bohr)")
    args = parser.parse_args()
    p = json.load(open(args.input))

    mass = (p["Mass1"] + p["Mass2"]) * AMU
    reduced = p["Mass1"] * p["Mass2"] / (p["Mass1"] + p["Mass2"]) * AMU
    k_min = sqrt(2 * mass * p["EMinWP"] / HARTREE_EV)
    k_max = sqrt(2 * mass * p["EMaxWP"] / HARTREE_EV)
    k0 = (k_min + k_max) / 2
    sigma_k = (k0 - k_min) / 2  # = sqrt(alpha), density width
    sigma_z = 1 / (2 * sigma_k)
    dz, z_start, tol = p["dZ"], p.get("ZStart", 0.0), p.get("TolSp", 1e-7)
    k_grid = np.pi / dz
    z_end = z_start + (p["nZ"] - 1) * dz
    zsp_end = z_start + (p["nZSpecular"] - 1) * dz

    print(f"packet {p['EMinWP']}-{p['EMaxWP']} eV, M = {mass:.1f} me: k0 = {k0:.3f}, sigma_k = {sigma_k:.3f}, "
          f"sigma_Z = {sigma_z:.3f} bohr")
    n_grid = (k_grid - k0) / sigma_k
    print(f"grid momentum limit pi/dZ = {k_grid:.2f}: (k_grid - k0)/sigma_k = {n_grid:.1f} "
          f"-> wrapped outward {tail(n_grid):.1e}  [need >= 6 and << TolSp = {tol:g}]")
    n_out = k0 / sigma_k
    print(f"own outward tail: k0/sigma_k = {n_out:.1f} -> {tail(n_out):.1e}  [need >= 6, i.e. EMax <= 4 EMin]")
    print(f"scattering grid {z_start}-{z_end:.3f} bohr, specular grid {z_start}-{zsp_end:.3f} bohr")
    print(f"Zana = {p['Zana']} (Z_asymp {args.zasymp}); Z0 = {p['Z0']}: "
          f"{(p['Z0'] - p['Zana']) / sigma_z:.1f} sigma_Z above Zana [>= 8 recommended]")
    outer = p["Z0"] + 5.2 * sigma_z
    print(f"packet outer tail (1e-7) at {outer:.2f}; gap to ZSpOptStart {p['ZSpOptStart']}: "
          f"{p['ZSpOptStart'] - outer:.2f} bohr [>= ~2]")
    print(f"slowest part: {mass / k_min:.0f} a.u. per bohr")

    absorber_table("Z absorber (scattering grid)", p["ZOptStart"], p["ZOptEnd"], p["ZOptEnergy"], mass, z_end)
    absorber_table("Z absorber (specular grid)", p["ZSpOptStart"], p["ZSpOptEnd"], p["ZSpOptEnergy"], mass, zsp_end)
    r_end = p["RStart"] + (p["nR"] - 1) * p["dR"]
    absorber_table("r absorber (reduced mass)", p["ROptStart"], p["ROptEnd"], p["ROptEnergy"], reduced, r_end)


if __name__ == "__main__":
    main()
