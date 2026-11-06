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
| r₀, Δr, N_r | 0.5, 0.20, 40 | `RStart`, `dR`, `nR` |
| N_x (N_y) | 12 | `nX`, `nY` |
| a | 4.84 a₀ | fixed by the PES |
| j_max | 28 | `JMax`, `MJMax` |
| T | 25 000 a.u. | `tEnd` |
| Z₀, ΔZ, N_Z | 0.0, 0.15, 128 | **not used**: ΔZ = 0.12, N_Z = 160 (below) |

Table 5.1 gives N_tot = 321 million = 435 rotational functions (even j ≤ 28, all m_j) × 128·40·12·12,
so its N_Z is the scattering grid. Here N_tot = 435 × 160·40·12·12 = 400.9 million.

## Z grid: why ΔZ = 0.12 instead of 0.15
General rules: `../README.md`; check with `tools/wave_packet_check.py`. The packet width is not an input: it follows from `EMinWP`/`EMaxWP` (α = (k₀ − k_min)²/4, Fortran
`getAlphaCoeffWavePacket`), so its momentum spread σ_k = √α is fixed by the energy range. A grid
only holds momenta up to π/ΔZ. Whatever lies above wraps around in the FFT to large *outward*
momentum: it runs into the specular absorber, and it stops the dynamical specular analysis line
(`DynamZSpAna`) at once, because that line is fixed as soon as norm moves outward past it.

| ΔZ | Energy range (eV) | σ_Z (a₀) | k_max = π/ΔZ | Norm wrapped to outward momentum |
|---|---|---|---|---|
| 0.15 | 0.35–1.05 (first version of this input) | 0.28 | 20.9 | 8e-6 |
| 0.15 | 0.40–1.00 | 0.33 | 20.9 | 3e-7 |
| 0.14 | 0.35–1.05 | 0.28 | 22.4 | 1.3e-7 (above `TolSp` = 1e-7) |
| 0.12 | 0.40–0.95 (**this input**) | 0.36 | 26.2 | 1e-20 |

The first version (ΔZ 0.15, 0.35–1.05 eV, `Z0` 17) lost 3.3e-6 of norm in the first 500 a.u. and
fixed the specular analysis line at 18.75 a₀, behind the packet. With this input the loss after
500 a.u. is 4.7e-8 and the line keeps moving.

## Deduced or chosen (not in the thesis)
- **Time step** 2.5 a.u. (10 000 steps), as in chapter 6 (Table 6.1). Chapter 4 (Table 4.5) used
  5.0 a.u. with the same SPO code; the smaller step is the safer choice. `dtAna` 20 a.u. (chapter 4).
- **Wave packet** 0.40–0.95 eV (the thesis figure starts at 0.4 eV), 111 analysis energies
  (0.005 eV apart). The packet's momentum density at `EMinWP`/`EMaxWP` is e⁻² ≈ 14% of its peak.
- **Grids in Z** (ΔZ 0.12, Z from 0.0): scattering grid 160 points (to 19.08 a₀), specular grid
  256 points (to 30.60 a₀; 2D, so its cost is negligible). Both are valid FFT sizes.
- **Analysis line** `Zana` = 15.0 = Z_asymptotic (the DZ PES is cut off there; the lowest allowed
  value), grid point 125.
- **Initial position** `Z0` = 18.0: 8σ_Z above `Zana`; the tail at the 1e-7 level ends at
  Z0 + 1.85 = 19.85, 2.2 a₀ below the specular absorber. Further out only costs time (the
  slowest analysed part, 0.40 eV, needs ~350 a.u. per a₀).
- **Absorbers** (quadratic, strength from the Vibók–Balint-Kurti table at the target energy
  `…OptEnergy`). The target energy is a tuning knob, not the packet energy. A low target leaves the
  fast part of the wave transmitted through the absorber, and the FFT grid then wraps it around to
  Z = 0. A high target reflects the slow part. The absorbers end on the last grid point (the
  absorbing potential is zero beyond its end). 1D reflection R and transmission T, computed for
  these exact absorbers:

  | Absorber | Range (a₀) | Target | R + T | R at 0.02 eV |
  |---|---|---|---|---|
  | Z, scattering | 15.0–19.08 | 0.2 eV | ≤ 4e-5 for 0.05–0.95 eV (8e-6 at 0.95) | 2e-2 |
  | Z, specular | 22.08–30.60 | 0.4 eV | ≤ 4e-9 for 0.2–0.95 eV; 2e-8 at 1.2 eV | 1e-4 |
  | r | 4.1–8.3 | 0.3 eV | ≤ 2e-5 above 0.2 eV; 1e-3 at 0.1 eV | 0.13 |

  The scattering absorber sees every channel, down to the thresholds (hence the low target). The
  specular grid carries only the initial channel, elastically reflected at the packet energies,
  so its absorber is tuned higher. Waves below ~0.03 eV normal energy (diffraction or rotational
  channels just above threshold) are reflected whatever the choice, but they move slowly and
  barely return within T.

## Cost (this machine)
6.0 GiB wave function, 1.3 GiB exp(−iV dt) table (`UseSymmetry`, C6v: 19 of 144 (x, y) points).
CPU (8 threads): 7.3 GiB during the propagation, about 1.25× the time of the 128-point grid
(~5.5 s per step, ~15 h). GPU (`--device cuda`): 8.8 GiB estimated plus the CUDA context; fits
the 12 GB RTX 4070 Ti.

## Fortran reference run
`DVR.inp` is the same input for the Fortran executable (keep it in sync with `input.json`; the
title is shorter because Fortran strings hold 80 characters). The executable reads `DVR.inp`
from the current directory and writes there: copy it into an empty run directory and start the
executable from that directory.
