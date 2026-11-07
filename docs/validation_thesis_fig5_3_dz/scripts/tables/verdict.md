| # | Item | Verdict | C++ CPU − Fortran | C++ GPU − Fortran | C++ GPU − C++ CPU | Criterion |
|---|---|---|---|---|---|---|
| [1](stages.md#item-1) | Input read, decisions taken | identical / printed | 54/54 echoed values = input (p); input.json ≡ DVR.inp 55/56 keys (Title by design) | same | shared code | exact; ½ printed unit |
| [2](stages.md#item-2) | Grids (Z, r, x, y, k, θ, φ, J) | identical / rounding | grids, φ, counts identical; θ nodes ≤ 6 ulp (f) | same | shared code | exact; few ulp |
| [3](stages.md#item-3) | MISC SETUP lines | identical / printed | identical 8, n/a 1, printed 19 | same | shared code | exact; ½ printed unit |
| [4](stages.md#item-4) | 6D potential, full grid | rounding | 2.8e-14 at the same point, 7.5e-14 as used (f) | same (V computed on the host) | identical (bitwise) | 1e-13 · max(\|V\|, 1 eV) |
| [5](stages.md#item-5) | Gas-phase V_asy(r) | rounding | 5.1e-16 (f) | same | identical | 1e-14 · max(\|V\|, 1 eV) |
| [6](stages.md#item-6) | Rovibrational states | printed / rounding | 348 values ≤ 0.99 ½ unit (p); E(0,0) 7.1e-15 rel. (f) | same | identical | ½ printed unit; 1e-12 rel. |
| [7](stages.md#item-7) | Diffraction basis (incl. warning) | identical | identical 5 | same | shared code | exact |
| [8](stages.md#item-8) | Rotation channels | identical | identical 4 | same | shared code | exact |
| [9](stages.md#item-9) | Rotation-diffraction channels | identical / rounding | labels of 9504 states identical; channel energies ≤ 1.3e-14 rel. (f) | same | identical | exact; 1e-13 rel. |
| [10](stages.md#item-10) | Optical potentials | printed | printed 6 | same | shared code | ½ printed unit |
| [11](stages.md#item-11) | Initial wave packet | printed / rounding / explained | explained 2, identical 5, printed 3, rounding 3 | same | shared code | ½ printed unit; 1e-13 rel. |
| [12](stages.md#item-12) | Analysis set-up, energy grid | identical / printed | identical 5, printed 1; 111/111 energies (p) | same | identical | exact; ½ printed unit |
| [13a](stages.md#item-13a) | Averages every 500 a.u. | printed | 144/144 within ½ unit (p) | 144/144 (p) | norm 5.0e-12 rel. (f) | ½ printed unit; 1e-10 |
| [13b](stages.md#item-13b) | Specular analysis line fixed | identical | same step (t = 2220 a.u.) | same step | same step | exact |
| [13c](stages.md#item-13c) | BK coefficients, 875 steps | rounding · open (channel-local excess) | 6.1e-12 (f) | 7.3e-12 (f) | 1.2e-12 (f) | 1e-10 · max\|c\| |
| [13d](stages.md#item-13d) | Final wave function | rounding | ‖Δψ‖/‖ψ_T‖ 2.0e-10, 1.2e-13 of the initial norm (f) | 1.2e-13 (f) | 1.9e-14 (f) | reported (no threshold) |
| [14](stages.md#item-14) | P_reaction, ΣP, P(v), P(v, j, mJ) | printed / explained / rounding | P_reaction 109/111 within ½ unit, 2 explained (p); full precision ≤ 1.4e-11 (f) | ≤ 1.6e-11 (f) | ≤ 1.9e-12 (f) | ½ printed unit, floor 1e-11 for P_reaction; 1e-10 |
| [15](stages.md#item-15) | Diffraction-order probabilities | printed / explained | 82432/82488 within ½ unit, 56 explained (p) | 82430/82488 within ½ unit, 58 explained (p) | 6.4e-13 (f) | ½ printed unit; 1e-10 |
| [16](stages.md#item-16) | S-matrix elements | printed / explained / rounding · open (channel-local excess) | \|S\|² 636327/636373 within ½ unit (p); \|ΔS\| ≤ 1.7e-12 (f) | \|ΔS\| ≤ 2.2e-12 (f) | \|ΔS\| ≤ 5.5e-13 (f) | ½ printed unit; 1e-10 |
| [17](stages.md#item-17) | Thesis Fig. 5.3 (dotted) | agrees | 67/67 points; max \|ΔP\| 0.014 | 67/67 | — | digitisation + 0.01 |
| [18a](stages.md#item-18a) | Fine analysis (onlyAnalysis) of each BK file | rounding | P_reaction ≤ 1.4e-11 (f) | ≤ 1.5e-11 (f) | ≤ 1.8e-12 (f) | 1e-10 |
| [18b](stages.md#item-18b) | Replay through one analysis code | rounding / printed | propagation only: \|ΔP\| ≤ 9.5e-13 (f); analysis only: all printed P within ½ unit (p) | \|ΔP\| ≤ 1.2e-12 (f) | \|ΔP\| ≤ 2.7e-13 (f) | 1e-10; ½ printed unit |
| [19](stages.md#item-19) | Physics checks | same in all codes | thresholds, bounds, P(mJ) = P(−mJ) ≤ 4.7e-12 rel. | same | same | property shared by all three |
