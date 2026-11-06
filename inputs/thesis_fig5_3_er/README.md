# Somers PhD thesis, Figure 5.3: H₂(v=0, j=0) + Cu(111), ER PES

The solid curve of Figure 5.3: the ER PES is M. Persson et al., J. Chem. Phys. 110, 2240 (1999),
Fortran module `PESForCu111Persson` (SystemID `H2Cu111Persson`, a = 2.5173 Å = 4.757 a₀,
Z_asymptotic = 8.0 a₀). The reasoning (why ΔZ = 0.12, the wave-packet width, the absorbers) is in
`../thesis_fig5_3_dz/README.md`; this case follows the same rules (`../README.md`; check with `tools/wave_packet_check.py`).

## Taken from Table 5.1 (column "ER PES, H₂")
r: 0.5 + 32 × 0.25, N_x = N_y = 16, j_max = 28, T = 15 000 a.u. Table 5.1's Z grid (0.0 + 90 × 0.15,
N_tot = 435 × 90·32·16·16 = 321 million) is **not used**: ΔZ = 0.15 cannot hold the momentum spread
of the packet (8e-6 of the norm wraps around to outward momentum for 0.35–1.05 eV, see the DZ
README). Here N_tot = 435 × 112·32·16·16 = 399 million.

## Deduced or chosen
- Time step 2.5 a.u. (6000 steps, as in chapter 6), `dtAna` 20 a.u. (chapter 4).
- Wave packet 0.40–0.95 eV, 111 energies (0.005 eV apart), `Z0` = 11.0 (8σ_Z above `Zana`; the
  1e-7 tail ends at 12.85, 3.2 a₀ below the specular absorber).
- Z grids (ΔZ 0.12, Z from 0.0): scattering 112 points (to 13.32 a₀), specular 210 points
  (to 25.08 a₀).
- `Zana` 8.04 (first grid point at or beyond Z_asymptotic).
- Absorbers (1D R + T for these exact absorbers):

  | Absorber | Range (a₀) | Target | R + T | R at 0.02 eV |
  |---|---|---|---|---|
  | Z, scattering | 8.04–13.32 | 0.2 eV | ≤ 7e-6 for 0.05–0.95 eV (4e-6 at 0.95) | 4e-3 |
  | Z, specular | 16.08–25.08 (L/λ = 14.9) | 0.4 eV | ≤ 3e-9 for 0.2–1.0 eV; 2e-8 at 1.2 eV | 4e-5 |
  | r | 4.0–8.25 | 0.3 eV | ≤ 2e-5 above 0.2 eV; 1e-3 at 0.1 eV | 0.13 |

- No `UseSymmetry`: the Persson PES is not periodic on the grid (it converts coordinates with
  0.529167 Å/bohr but the lattice constant with 0.529177), so the symmetry self-check fails.

## Cost (this machine)
CPU (8 threads): 5.9 GiB wave function plus 3.7 GiB potential, 9.7 GiB in total. exp(−iV dt)
exceeds the 4 GiB table limit and is computed in every step. With the 90-point grid a step took
about 5.8 s (static schedule); this grid has 1.24× the points in Z.
GPU: does not fit a 12 GB card (estimate 14.8 GiB, of which 7.4 GiB is the exp(−iV dt) table:
without symmetry it holds all 256 (x, y) points). The 90-point grid did not fit either (~12.3 GiB).
