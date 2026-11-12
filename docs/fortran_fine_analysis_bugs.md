# Fortran fine analysis: interpolation bugs and their fix

The fine analysis (Fortran `nread = 2`, input `onlyAnalysis` / `alsoAnalysis`) has two steps:
1. It computes the S-matrix on the spectral grid E_k = k·π/T of the BK time series
   (`gtaife`, `hlftrf`, `probaf`).
2. It interpolates the grid values linearly to the NrOfEs analysis energies (`StateToStateProbs`,
   `smatrix_calc`).

The interpolation in the Fortran code (`Analyzor.f90 smatrix_calc`, nread = 2 branch without NAG,
lines 732–790) has two bugs, and its sign handling causes a third class of errors. The C++ code
fixes all three. The probabilities are not affected; the S-matrix elements and phases are.
Line numbers refer to version 5.19, and so do the counts of wrong values in the last table
(the 5.19 model with all three bugs). Upstream fixed bug 1 on 2026-11-09
(svn r1713, after version 5.24); that is the reference of the port
([validation.md](validation.md#reference-data)). Bugs 2 and 3 are still present there.

## The Fortran code
```fortran
interpolation_factor = ( etransZ - ftrans(idat-1) ) / ( ftrans(idat) - ftrans(idat-1) )
...
prob1 = s1_real**2 + s1_imag**2
prob2 = s2_real**2 + s2_imag**2
prob_interpolated = prob1 + (prob2 - prob1) * interpolation_factor

if ( s1_real .RNE. 0.0 ) then
  phase1 = ATAN(s1_imag/s1_real)                     ! (2) ATAN, not ATAN2
else
  phase1 = 0.0
endif
(same for phase2)
phase_interpolated = phase1 + (phase2 - phase1) * interpolation_factor

sreal_interpolated = s1_real + (s2_real - s1_real) * interpolation_factor
simag_interpolated = s1_imag + (s2_imag - s1_real) * interpolation_factor   ! (1) s1_real
sreal = sqrt(prob_interpolated) * cos(phase_interpolated)
simag = sqrt(prob_interpolated) * sin(phase_interpolated)
sreal = sign(sreal, sreal_interpolated)              ! (3) signs taken separately
simag = sign(simag, simag_interpolated)
smtrx(is,iv,knt) = cmplx( sreal, simag )
```
s1 and s2 are the S-matrix elements at the grid points below and above E_z. The intent is:
- interpolate |S|² linearly;
- interpolate the phase linearly;
- rebuild S from both.

## What goes wrong
1. **Re S in place of Im S (line 778).** `simag_interpolated` should be
   `s1_imag + (s2_imag - s1_imag)·f`.
   - The value only sets the sign of Im S, so the result has the wrong sign of Im S whenever the
     two expressions have different signs.
   - This happens even on a grid point (f = 1). The sign then comes from s2_imag + s1_imag − s1_real.
2. **ATAN instead of ATAN2 (lines 761, 768).**
   - ATAN gives phases in (−π/2, π/2) only. Where Re S changes sign between the two grid points,
     φ1 and φ2 jump by about π, so the interpolated phase is wrong by O(1).
   - The split of √P between |Re S| and |Im S| is then wrong. The sign correction of (3) cannot
     repair a wrong magnitude.
   - The `.RNE. 0.0` test is the absolute 1e-12 comparison of `MyOperators.f90`.
   - ATAN2 alone would not be enough: the phase would still need unwrapping across ±π.
3. **Signs taken separately (lines 785–786).** The magnitude |sin φ_int| and the sign of the
   linear interpolation of Im S have their zeros at different f whenever |s1| ≠ |s2|. So where
   Im S (or Re S) changes sign between the grid points, there is an interval of f with the wrong
   sign.

**Not affected:**
- P = sreal² + simag² = `prob_interpolated` up to rounding. Every probability is therefore
  unaffected: the stdout per-energy blocks, P(v), P(v, j, mJ), the diffraction orders, the
  reaction probability, and all average energies.
- The spectrum on the grid itself (probaf: fil10, unit 3) is unaffected; probaf uses ATAN2.

**Affected:** Re S, Im S and the phase at the analysis energies: Fortran fort.25, C++
`fine_analysis/s_matrix.tsv`.

## The C++ interpolation
`interpolate_s_matrix` (`src/analysis/fine_analysis.cpp`):
- P = P1 + (P2 − P1)·f, as in the Fortran, so the probabilities stay comparable.
- φ = arg(s1 + (s2 − s1)·f): the phase of the linearly interpolated S.
  - It is continuous across ±π and needs no special case for a zero end point (a closed channel).
  - It fails only when the phase changes by about π between neighbouring grid points.
- S = √P·e^{iφ}.

If a populated channel (|S|² > 1e-8 of the largest |S|²) changes its phase by more than π/2
between neighbouring grid points, the program prints a warning. The interpolation is then
unreliable, and the run should be propagated longer.

## How the fix is verified (the Fortran code is not modified)
The reference data comes from the unmodified Fortran binary. It runs as an analysis-only run on
the BK coefficients of its own full hardwall_small run:
`tools/fortran_reference/run_fine_analysis_reference.sh`, variants `full`, `truncated`, `nen95`,
`small` and `empty`. The C++ analyses the same file, so both programs start from identical
input.

1. **The fine-grid S agrees with the Fortran, and the bug model is what the Fortran does**
   (`tests/test_regression_fine_analysis.cpp`, gate 5).
   - The test helper `fortran_interpolation_with_bugs` emulates the reference: lines 739–787
     literally, with 5.24's `s1_imag` in place of the typo (bugs 2 and 3 remain). The unit tests
     of item 2 keep the 5.19 model with all three bugs (`fortran519_interpolation_with_bugs`).
     Applied to the C++ grid values, it reproduces every Re S, Im S and phase in the Fortran
     fort.25 to within half a unit in the last printed digit (counts printed by
     `ctest -R Gate5 -V`, 5.24 r1713 reference data, 2026-11-12):

     | variant | values compared |
     |---|---|
     | full | 144 |
     | nen95 | 243 |
     | small | 297 |
     | Dai–Light full (`cu111_dailight_small`) | 105303 |

     Values with |S| < 1e-8 of the largest are skipped.
   - In the nen95 variant the analysis energies fall on grid points (f = 0 or 1). There the
     Fortran's quad-precision `ensctr` picks the neighbouring interval for 7 values; the test
     accepts either interval.
   - All probabilities agree with the fixed C++ directly, to within half a printed unit:
     - full: 1365 printed in the per-energy blocks, 9904 in fort.25;
     - plus every value of the spectrum (fil10, unit 3, probaf table).
2. **The differences occur only where predicted**
   (`tests/test_fine_analysis.cpp`, `FortranInterpolationBugs.*`).
   - The emulated 5.19 result (`fortran519_interpolation_with_bugs`, all three bugs) is compared
     with the intended one (P and the unwrapped phase interpolated linearly).
   - Every difference falls into one of three classes:
     - Re S changes sign between the grid points (2);
     - Im S changes sign (3);
     - the typo changes the sign of the interpolated Im S (1).
   - Everywhere else the two agree to rounding.
3. **The C++ result is close to the exact S; the Fortran result is not.** The exact S at an
   analysis energy is the time-to-energy transform (trapezoid sum) evaluated at that energy, with
   no interpolation. The table gives the largest relative error |S − S_exact|/|S_exact|:

| case | pairs | Fortran wrong (class 1 / 2 / 3) | C++ error | Fortran error |
|---|---|---|---|---|
| hardwall_small, Fortran BK file, 39 analysis energies | 39 | 2 (2 / 0 / 0) | 7.4e-4 | 0.31 (phase 0.31 rad) |
| synthetic S(E), phase step 0.9 rad per grid point | 200 | 96 (36 / 60 / 36) | 0.039 | O(1) |

The hardwall_small row counts only pairs with |S| > 1e-3 of the largest at both grid points.
The C++ error there is the interpolation error of |S|² (9.8e-4 relative), which the Fortran
shares.

For hardwall_small the phase of the elastic channel advances only about 0.13 rad per grid point
(6.17 meV). Re S and Im S therefore rarely change sign between neighbours, and only the typo
(class 1) occurs. The synthetic case has a fast phase, so all three classes occur.
