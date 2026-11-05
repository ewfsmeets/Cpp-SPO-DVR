# Somers PhD thesis, Figure 5.3: H₂(v=0, j=0) + Cu(111), DZ PES

Reproduces the dotted curve of Figure 5.3 (M. F. Somers, PhD thesis, Leiden, chapter 5): reaction
probability against normal incidence energy (0.4–1.0 eV) on the Dai–Zhang LEPS PES. The solid curve
(ER/Persson PES) is `../thesis_fig5_3_er/`.

"DZ PES" in the thesis is J. Dai and J. Z. H. Zhang, J. Chem. Phys. 102, 6280 (1995). That is the
Fortran module `PESForCu111DaiLight` (SystemID `H2Cu111DaiLight`, lattice constant 4.84 a₀). The
module's header comment attributes it to "Dai and Light", but the reference is Dai and Zhang.

## Taken from Table 5.1 (column "DZ PES, H₂")
| Table 5.1 | Value | Input key |
|---|---|---|
| Z₀, ΔZ, N_Z | 0.0, 0.15, 128 | `ZStart`, `dZ`, `nZ` (scattering grid) |
| r₀, Δr, N_r | 0.5, 0.20, 40 | `RStart`, `dR`, `nR` |
| N_x (N_y) | 12 | `nX`, `nY` |
| a | 4.84 a₀ | fixed by the PES |
| j_max | 28 | `JMax`, `MJMax` |
| T | 25 000 a.u. | `tEnd` |

Check: N_tot = 435 rotational functions (even j ≤ 28, all m_j) × 128·40·12·12 = 320.7 million. This
matches "321 million" in the table, so N_Z is the scattering grid. The program reports the same
size.

## Deduced or chosen (not in the thesis)
- **Time step** 2.5 a.u. (10 000 steps), as in chapter 6 (Table 6.1). Chapter 4 (Table 4.5) used
  5.0 a.u. with the same SPO code; the smaller step is the safer choice. `dtAna` 20 a.u. (chapter 4).
- **Wave packet** 0.35–1.05 eV, a margin around the 0.40–1.00 eV of the figure; 141 analysis
  energies. `Z0` = 17.0: the packet must start beyond Z_asymptotic = 15 (the DZ PES is cut off there),
  which is why this grid reaches 19 a₀ and T is longer than for the ER PES.
- **Analysis line** `Zana` = 15.0 (= Z_asymptotic, the lowest allowed value).
- **Specular grid** 192 points (to 28.65 a₀; not in the thesis; its memory cost is negligible).
- **Absorbers** (quadratic, strength from the Vibók–Balint-Kurti table at the target energy
  `…OptEnergy`). The target energy is a tuning knob, not the packet energy. A low target leaves the
  fast part of the wave transmitted through the absorber, and the FFT grid then wraps it around to
  Z = 0. A high target reflects the slow part. 1D reflection R and transmission T (computed for these
  exact absorbers):

  | Absorber | Range (a₀) | Target | R + T, 0.05–1.0 eV | R at 0.02 eV |
  |---|---|---|---|---|
  | Z, scattering | 15.0–19.05 | 0.2 eV | ≤ 5e-5 (1e-5 at 1 eV) | 2e-2 |
  | Z, specular | 19.05–28.65 | 0.2 eV | ≤ 2e-6 | 6e-6 |
  | r | 4.1–8.3 | 0.3 eV | ≤ 2e-5 above 0.2 eV; 1e-3 at 0.1 eV | 0.13 |

  For comparison, the 0.05 eV target of chapter 4 transmits 2% at 1 eV on this short (4 a₀)
  absorber. Waves below ~0.03 eV normal energy (diffraction or rotational channels just above
  threshold) are reflected whatever the choice, but they move slowly and barely return within T.

## Cost (this machine, 8 threads)
4.8 GiB wave function, 6.3 GiB in total (`UseSymmetry`, C6v: 19 of 144 (x, y) points). About
4.4 s per step (measured before the switch to a dynamic schedule), so the 10 000 steps take about 12 h.
