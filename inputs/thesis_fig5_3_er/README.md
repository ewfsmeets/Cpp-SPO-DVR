# Somers PhD thesis, Figure 5.3: H₂(v=0, j=0) + Cu(111), ER PES

The solid curve of Figure 5.3: the ER PES is M. Persson et al., J. Chem. Phys. 110, 2240 (1999),
Fortran module `PESForCu111Persson` (SystemID `H2Cu111Persson`, a = 2.5173 Å = 4.757 a₀,
Z_asymptotic = 8.0 a₀). The reasoning and the absorber analysis are in `../thesis_fig5_3_dz/README.md`.

## Taken from Table 5.1 (column "ER PES, H₂")
Z: 0.0 + 90 × 0.15 (`nZ` 90), r: 0.5 + 32 × 0.25, N_x = N_y = 16, j_max = 28, T = 15 000 a.u.
N_tot = 435 × 90·32·16·16 = 320.7 million (the table's "321 million").

## Deduced or chosen
- Time step 2.5 a.u. (6000 steps, as in chapter 6), `dtAna` 20 a.u. (chapter 4). Wave packet 0.35–1.05 eV at
  `Z0` = 10.0, 141 energies.
- `Zana` 8.1 (first grid point beyond Z_asymptotic). Specular grid: 160 points (to 23.85 a₀).
- Absorbers (1D R + T over 0.05–1.0 eV):

  | Absorber | Range (a₀) | Target | R + T, 0.05–1.0 eV | R at 0.02 eV |
  |---|---|---|---|---|
  | Z, scattering | 8.1–13.35 | 0.2 eV | ≤ 7e-6 | 5e-3 |
  | Z, specular | 12.0–23.85 | 0.2 eV | ≤ 1.3e-6 | 3e-6 |
  | r | 4.0–8.25 | 0.3 eV | ≤ 2e-5 above 0.2 eV; 1e-3 at 0.1 eV | 0.13 |

- No `UseSymmetry`: the Persson PES is not periodic on the grid (it converts coordinates with
  0.529167 Å/bohr but the lattice constant with 0.529177), so the symmetry self-check fails.

## Cost (this machine, 8 threads)
4.8 GiB wave function plus 3.0 GiB potential, 7.8 GiB in total. exp(−iV dt) exceeds the 4 GiB table
limit and is computed in every step: about 5.8 s per step (static schedule), so the 6000 steps take
about 10 h.
