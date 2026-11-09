# Validation against the Fortran code

The C++ program is checked against the original Fortran program (`~/work/SPO-DVR`, version 5.19,
gfortran build) by five acceptance gates. Each gate is a regression test in `tests/` that compares
the C++ results with Fortran output stored in `tests/reference/<case>/`; see the README section
"Testing" for how to run them, and `tools/fortran_reference/` for regenerating the Fortran data.

**Tolerances.** Values that the Fortran prints in E/F format (e.g. E15.8, mantissa 0.dddddddd) are
compared to half a unit in the last printed digit, not with a fixed relative tolerance (for
0.147E+01 the last digit is worth 3.4e-8 relative). Tables are read from the dedicated output files
(e.g. `averages.out`): the stdout contains other tables in the same format. Where a tolerance is
max(half a printed unit, floor), raising the floor changes nothing while the printed unit binds.

Known differences between the two programs (Fortran bugs fixed in C++, quirks kept on purpose,
quad-precision sites) are listed in [fortran_deviations.md](fortran_deviations.md).

**Production grid.** The gates use small grids. A full thesis calculation (`inputs/thesis_fig5_3_dz`,
Dai–Zhang PES, 7000 steps, 435 rotational functions) was run with the Fortran, the C++ CPU and the
C++ GPU code and compared state by state, with criteria written before the comparisons:
[validation_thesis_fig5_3_dz](validation_thesis_fig5_3_dz/README.md). With the default (strict)
GPU build the S matrix agrees with the Fortran to 3.5e-13 (norm-weighted relative) and all 111
reaction probabilities match the Fortran's 7 printed digits.

## Gates
Main reference case: `tests/reference/hardwall_small/DVR.inp` (Fortran input), with the Fortran
output next to it; the corrugated Cu(111) cases are in their own directories.
**All gate values below are specific to the parameters in that `DVR.inp`** (grids, masses, JMax, …).
If `DVR.inp` changes, the Fortran reference must be regenerated and the gate values updated; the
C++ regression tests must use the equivalent JSON input (`tests/reference/hardwall_small/input.json`)
derived from exactly the same parameters.
1. **Hard wall PES, gas-phase potential** (`tests/test_regression_hardwall.cpp`): V(r) at Z = zasymp = 7.0, θ = π/2 (Fortran
   "INFO FROM SETTING UP THE ROVIBRATIONAL ASYMPTOTIC POTENTIAL" table, r = 0.8 … 9.05, 56 points)
   reproduced to ~1e-14 relative (full 6D potential ~1e-13, see the quad-precision sites in
   fortran_deviations.md). Sanity: V(r = 1.4) = 0 exactly. With GCC 15 the fast-math builds
   (`SPODVR_FAST_MATH`) fail the V(r) table by 1.6e-14 (at r index 5); the default build passes.
2. **Gas-phase rovibrational states** (same test; Fortran `GasPhase.f90` / `FourierGridHamiltonian.f90`): energies
   e(v,j) for v = 0,1 and j = 0…4 (Fortran "INFO FROM SETTING UP THE ASYMPTOTIC ROVIBRATIONAL
   STATES"), e.g. e(0,0) = 0.258378 eV, e(1,0) = 0.762842 eV; also ⟨r⟩, psi(1), psi(nr), norms.
   The printed Fortran values have only 6 digits — dump full precision from Fortran when a tighter
   gate is needed.
3. **Propagation** (`tests/test_regression_propagation.cpp`): initial wave-packet block
   (norm, |ψ(k)| values, iLow/iHigh, norm inside) and the averages table (time, norm, ⟨Z⟩, ⟨r⟩,
   rest norm) against the Fortran output (stdout; averages from averages.out), to all printed
   digits: |C++ − Fortran| ≤ half a unit in the 8th mantissa digit of E15.8 (the Fortran is not
   patched). Passes for hardwall_short (21 rows) and the full hardwall_small run (7000 steps, 141
   rows). Default case `tests/reference/hardwall_short` (same as
   hardwall_small but 1000 steps); the full 7000-step hardwall_small run only with
   `SPODVR_LONG_TESTS=1`. The hard wall has no θ/φ/x/y coupling, so gate 3 only checks that the
   angular transform and x/y FFTs act as the identity there; their correctness rests on the
   analytic unit tests (test_angular, test_propagation). Corrugated cases (60° cell, C6v):
   `cu111_dailight_short` (1000 steps, default), `cu111_dailight_short_sym` (C6v UseSymmetry,
   nPhisOptimal: nφ 12), `cu111_dailight_small` (7000 steps, long; also fine analysis, gate 5
   `Gate5DaiLight`). Status: all pass (gates 3, 4, 5). Full run: averages worst 0.99 of half a
   unit, BK (1400 records) 6.4e-13 of max|c| = 0.26, probabilities within half a unit except the
   reaction probability (see gate 4). Run times of this case: performance.md "Fortran vs C++".
   fast-math (`SPODVR_FAST_MATH`) moved the reaction
   probability to 2.0e-12 / 5.3e-12 from the Fortran (old floor 1e-12): floor raised to 1e-11;
   not the default build.
   Setup gates (gates 1/2 type, `tests/test_regression_pes.cpp`) for cu111_dailight_small (with
   switching 6..8), cu111_persson_small and cu111_dailight_short_sym (the Fortran symmetric .pes,
   expanded through symmetryIndex by extract_pes.py, vs the C++ gather): all pass, vasy ≤ 9e-16,
   6D ≤ 3.6e-14 of max(|V|, 1 eV); computed (x, y) points 10 = 10 (C6v triangle 1).
   The same gates for the CRP library (SystemID `CRPlibrary`) on the synthetic data sets of its
   tests (pairwise model; `tests/crplib/reference/<sym>_pairwise/dataset`):
   `crp_c6v_synthetic_small`, `crp_c3v_synthetic_small` (6 × 6) and `crp_c4v_synthetic_small`
   (8 × 8), with r and Z partly outside the data (clamped) and switching 6..8: all pass, vasy
   bitwise, 6D ≤ 3.0e-14 of max(|V|, 1 eV) (the grids' quad-precision ulps, as above).
   With `UseSymmetry` (C++ only for C3v and C4v; the subgroup usable with nφ = 10) the same cases
   against the full Fortran potential: 6D ≤ 5.7e-14.
4. **On-the-fly analysis** (same test as gate 3, same runs): every probability printed in the
   Fortran per-energy blocks (total energy, sum, reaction probability, P(v), P(v, j, mJ); 1435
   values) to half a unit in the last printed digit (values below 1e-20 are noise and compared
   absolutely), every BK coefficient of every analysis step (`fortran_fort21.bin`) to 1e-10 of
   max(max|c|, 1e-4) (exception: the reaction probability = 1 − sum is compared to
   max(half a printed unit, 1e-11) absolute, the accuracy of the sum it is derived from; decided
   after cu111_dailight_small: 2.9e-6 printed to 7 digits differed by 3.6e-13; floor
   raised from 1e-12 to 1e-11 for the fast-math build (5.3e-12), see
   tests/analysis_reference.hpp), and
   the probabilities of a Fortran run split by a restart
   (`fortran_restart_stdout.txt`, tools/fortran_reference/run_restart_reference.sh).
   Status: all pass. Full run: probabilities worst 0.99 of half a unit, BK (1400 records) 6.4e-13
   of max|c| = 0.41; Fortran restart run worst 0.99. In hardwall_short all BK coefficients are
   rounding noise (≤ 2e-15, the packet has not reached Zana), hence the scale floor 1e-4.
5. **Fine analysis** (`tests/test_regression_fine_analysis.cpp`): the C++ analyses the Fortran's
   own BK file (hardwall_small `fortran_fort21.bin`); Fortran reference from unmodified
   onlyAnalysis runs (`tools/fortran_reference/run_fine_analysis_reference.sh`, variants full,
   truncated (tEnd 13860), nen95 (NrOfEs -95), small (tEnd 800), empty (tEnd 100, nen95 = 0) in
   `hardwall_small/fine_<variant>/`). Compared to half a printed unit: AsympInfo numbers,
   probaf table, all of fil10 and unit 3 (spectrum), the interpolated per-energy blocks and
   fort.25 probabilities; fort.25 Re S / Im S / phase against an emulation of the buggy Fortran
   interpolation (see fortran_deviations.md). Status: all pass (worst 0.99 of half a unit).
   Note: the associated Gauss-Legendre quadrature for |mJ| > 0 is dead code in the Fortran
   (`setupdvr.f90` forces `SymmetryInPhi = .FALSE.`, so only mJ = 0 is used). The rotational
   basis uses `legendre.f90 AssociatedLegendre` instead (C++ `associated_legendre`), tested
   analytically (closed forms, orthonormality, cos²θ matrix elements, transform round trip).
