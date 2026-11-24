# Validation against the Fortran code

The C++ program is checked against the original Fortran program (upstream SPO-DVR 5.24 at svn
r1727, built with `-O3`: [Reference data](#reference-data)) by five acceptance gates. Each gate is
a regression test in `tests/` that compares the C++ results with Fortran output stored in
`tests/reference/<case>/`; see the README section "Testing" and [CONTRIBUTING.md](CONTRIBUTING.md#testing) for how to run them.

**Tolerances.** Values that the Fortran prints in E/F format (e.g. E15.8, mantissa 0.dddddddd) are
compared to half a unit in the last printed digit, not with a fixed relative tolerance (for
0.147E+01 the last digit is worth 3.4e-8 relative). Tables are read from the dedicated output files
(e.g. `averages.out`): the stdout contains other tables in the same format. Where a tolerance is
max(half a printed unit, floor), raising the floor changes nothing while the printed unit binds.

**NaN.** A NaN on either side (C++ or Fortran) fails every gate. The largest deviations and the
reference scales are reduced with a maximum that keeps NaN (`nan_aware_max` in
`tests/comparison.hpp`; `std::max(x, NaN)` returns x), and each check on them fails for NaN
(`NaN ≤ tolerance` is false). Negative tests: `tests/test_comparison_helpers.cpp`.

Known differences between the two programs (Fortran bugs fixed in C++, quirks kept on purpose,
quad-precision sites) are listed in [fortran_deviations.md](fortran_deviations.md).

<a id="production-grid"></a>**Production grid.** The gates use small grids. A full thesis calculation (`inputs/thesis_fig5_3_dz`,
Dai–Zhang PES, 7000 steps, 435 rotational functions) was run with the Fortran, the C++ CPU and the
C++ GPU code and compared state by state, with criteria written before the comparisons:
[validation_thesis_fig5_3_dz](validation_thesis_fig5_3_dz/README.md). With the default (strict)
GPU build the S matrix of its propagation agrees with the Fortran's to 3.5e-13 (norm-weighted
relative; both runs' BK coefficients through the same C++ analysis), and all 111 reaction
probabilities lie within half a unit of the Fortran's 7th printed digit. That Fortran run used 5.19; on the
same grid, the potential and the first 40 steps of 5.24 built with `-O3` are bitwise identical to
5.19's (checked 2026-11-11), so the comparison holds for the current reference too.

## Reference data

The Fortran output in `tests/reference/`, `tests/crplib/reference/` and `tests/scm/reference/` (the
CRP and SCM driver outputs) comes from one build, the
reference install `r1727-3d00331`: upstream 5.24 plus the fine-analysis sign fix of svn r1713 and the
input-check and analysis fixes of r1714–r1727 (git `3d00331`, svn r1727), the commit the port is up to speed
with ([upstream_sync.md](upstream_sync.md)).

- **Why `-O3`:** the upstream default (`-Ofast -mtune=native -march=native`) brings fast-math and
  FMA (CPU-dependent results; the computed machine epsilon underflows to 0, so the code falls back
  to 1e-12). With `-O3`, 5.24 is bitwise identical to 5.19 on every gate case of 5.19 features
  (checked 2026-11-11).
- **Provenance:** the install's `BUILDINFO.tsv` (copies in `tools/fortran_reference/buildinfo/`)
  and a `MANIFEST.tsv` per case (`ctest -R ReferenceManifest`). Build, regenerate, compare, bump:
  [tools/fortran_reference/README.md](../tools/fortran_reference/README.md).
- **5.19** (install `v5.19`: symlinks into the unmodified 5.19 distribution) is the archive, kept
  for the thesis record.

### Reference history

| Period | Fortran build | Reference data |
|---|---|---|
| until 2026-11-06 | 5.19 (svn r1572), gfortran 13 build of 2026-10-26, `-O3` | generated, mostly at 16 threads |
| from the GCC 15 switch (2026-11-06) | 5.19 rebuilt with gfortran 15 | kept; the rebuild reproduces them except the Dai–Light BK files (≤ 8.3e-13 of the gate-4 scale) and last digits of fort.25 |
| since 2026-11-12 | 5.24 r1713 `-O3` (install `r1713-dce0ecc`), gfortran 15, 4 threads | all regenerated (`regenerate_all.sh`) and compared with `refdata.py compare … --expected tools/fortran_reference/expected/r1713.tsv`. Identical: the averages of all 5 propagation cases, the potentials and gas-phase states, the hard-wall BK files, the CRP driver outputs, fine-analysis fort.10 and fort.3. Changed by the toolchain (gfortran 13 → 15 and its run-time libraries): the Dai–Light BK files by ≤ 8.3e-13 (short), 4.1e-13 (short_sym), 1.1e-14 (small) of the gate-4 scale, the propagation fort.25 (read by no test), stdout probabilities below 1e-20 by ≤ 1e-21. Changed by the r1713 fix: the signs of Im S and of the phase in the fine-analysis fort.25. Stdout labels |
| since 2026-11-23 | 5.24 r1727 `-O3` (install `r1727-3d00331`), gfortran 15, 4 threads | all regenerated and compared with `refdata.py compare … --expected tools/fortran_reference/expected/r1727.tsv`: 993 of 1022 files identical (every BK file, averages, potential, gas phase, fine-analysis output, CRP and SCM driver output). Changed by r1720: in the propagation `fortran_fort25.txt` of the six cases whose specular analysis line becomes fixed, the initial state's rows (every v) rotate by −k_v (Z_sp − Zana), \|S\| and P unchanged (class `initial_phase`; hardwall_short identical, its line never fixes). Changed by r1725: the stdout label `(nproj = 1)`. Removed before the switch: the gate-5 variant `fine_nen95` (NrOfEs −95 is an input error since r1719) |

With the switch the C++ follows 5.24 twice: `RStart ≤ 0` is rejected (5.19: only `< 0`), and
gate 5 emulates the interpolation with 5.24's `s1_imag`; the unit tests `FortranInterpolationBugs.*`
keep the 5.19 model with its `s1_real` typo ([fortran_fine_analysis_bugs.md](fortran_fine_analysis_bugs.md)).

## Gates

Status on the 5.24 r1713 reference data (measured 2026-11-15 at 8b3df8d with 8 threads; the margins of
gates 3–5 below were measured at the switch, and the setup data of gates 1 and 2 did not change
with it): the full ctest, 809 tests = 803 passed + 6 skipped (the five long tests and
`CrpLib.Cases/FortranMoleculeSurface.GasPhase/"c6v_pairwise_angstrom"`, which has no gas-phase
reference); the long tests (`SPODVR_LONG_TESTS=1 ctest -R 'Long($|\.)'`), 5 of 5; the CUDA tests
including the long ones, 24 of 24 (gates 3 and 4 on the GPU; 7000-step BK within 8.5e-13 (hard
wall) and 9.3e-13 (Dai–Light) of the scale, measured 2026-11-12).

Main reference case: `tests/reference/hardwall_small/DVR.inp` (Fortran input), with the Fortran
output next to it; the corrugated Cu(111) cases are in their own directories.
**All gate values below are specific to the parameters in that `DVR.inp`** (grids, masses, JMax, …).
If `DVR.inp` changes, the Fortran reference must be regenerated and the gate values updated; the
C++ regression tests must use the equivalent JSON input (`tests/reference/hardwall_small/input.json`)
derived from exactly the same parameters.

1. **Potential** (`tests/test_regression_hardwall.cpp`, `tests/test_regression_pes.cpp`). Gate 1a:
   the gas-phase potential V(r) at Z = zasymp, θ = π/2 (the Fortran table "INFO FROM SETTING UP THE
   ROVIBRATIONAL ASYMPTOTIC POTENTIAL"). Gate 1b: the 6D potential (hartree, before the −dt
   premultiplication) against the Fortran's `.pes` file (`pes_reference.txt`: V(Z, r) for the hard
   wall, which depends on nothing else; a subsample of the 6D grid for the other cases), and the
   angular grid. Gate 1c (`RegressionPesCut`, in `tests/test_regression_pes.cpp`): the write2DPES cut,
   V(r, Z) of the PES itself (eV, no switching, no VMin/VMax) at a fixed x, y, θ, φ over every r and
   scattering-grid Z point, against the Fortran's 2dpes.dat (`pes_cut_<k>.dat`,
   `run_pes_cut_reference.sh`), for each cut of the case's `pes_cuts.json`.
   - Hard wall: V(r) at Z = zasymp = 7.0 (r = 0.8 … 9.05, 56 points) reproduced to ~1e-14 relative
     (full 6D potential ~1e-13, see the quad-precision sites in fortran_deviations.md). Sanity:
     V(r = 1.4) = 0 exactly. With GCC 15 and plain `-ffast-math` the fast-math builds
     (`SPODVR_FAST_MATH`) failed the V(r) table by 1.6e-14 (at r index 5); since NaN and inf are kept
     (`-fno-finite-math-only`, 2026-11-15) `build-fast` passes it, as the default build does.
   - Corrugated Cu(111) cases (`tests/test_regression_pes.cpp`, gates 1 and 2):
     cu111_dailight_small (with switching 6..8), cu111_persson_small and cu111_dailight_short_sym (the Fortran symmetric .pes,
     expanded through symmetryIndex by extract_pes.py, vs the C++ gather): all pass, vasy ≤ 9e-16,
     6D ≤ 3.6e-14 of max(|V|, 1 eV); computed (x, y) points 10 = 10 (C6v triangle 1).
   - The same gates for the CRP library (SystemID `CRPlibrary`) on the synthetic data sets of its
     tests (pairwise model; `tests/crplib/reference/<sym>_pairwise/dataset`):
     `crp_c6v_synthetic_small`, `crp_c3v_synthetic_small` (6 × 6) and `crp_c4v_synthetic_small`
     (8 × 8), with r and Z partly outside the data (clamped) and switching 6..8: all pass, vasy
     bitwise, 6D ≤ 3.0e-14 of max(|V|, 1 eV) (the grids' quad-precision ulps, as above).
     With `UseSymmetry` (C++ only for C3v and C4v; the subgroup usable with nφ = 10) the same cases
     against the full Fortran potential: 6D ≤ 5.7e-14.
     On the real CRP data set (H₂/Cu(111) SRP48, `potentials/H2Cu111SRP_C3v.dir`):
     `crp_cu111_srp48_small` (the upstream input on coarse grids, VMax imposed): vasy 7.7e-16, 6D
     1.5e-13 of max(|V|, 1 eV) with and without `UseSymmetry`. This case allows 3e-13: CRP D4
     on data up to 77 Hartree (1.5e-14 with gfortran's libmvec cos, checked 2026-11-11). Gate 1b
     clips the reference at VMin/VMax: the Fortran writes its .pes file before it imposes them.
   - SCM distorted-surface PES (SystemID `CRPDistrtSurface`): `scm_cu111_srp48_small` (gates 1–4;
     `crp_cu111_srp48_small` from ZStart 0.5 bohr ≥ the CRP's Zmin, plus the example slab of
     `potentials/H2Cu111SRP_C3v.dir` rescaled to the CRP lattice) and
     `scm_cu111_srp48_small_rescaling_off` (gates 1, 2). Vasy 7.7e-16 (ΔV = 0 there, SCM D3). Gate 1b
     allows 1e-12 of max(|V|, 1 eV): V_CRP and ΔV reach ±125 eV where V is a few tenths of an eV, and
     their absolute differences do not shrink with V. The 42620 compared points show no cancellation
     (1.39e-13 and 1.33e-13 at V_CRP 0.47, ΔV 0.59 eV); the 1e-12 rests on the whole .pes (3.94 M
     points; 2026-11-15, .pes of the 8222b8d runs, extract_pes.py at stride 1, C++ Release build):
     2.98e-13 (margin 3.4) and 3.41e-13 at V = −0.22 = 125.23 − 125.44 eV, 1.5e-13 scaled by max(|V|,
     |V_CRP|, |ΔV|, 1 eV). The terms are checked separately at identical points (`tests/scm`: ΔV
     bitwise, V_CRP within 3e-13). On that grid V − ΔV is up to 1.14e-12 of max(|V_CRP|, 1 eV) from
     the C++ V_CRP (r = 8.6 bohr, ΔV −329 eV): ΔV one ulp lower in r (the Fortran's grid is
     quad-rounded) moves by −1.42e-12 eV (unscaled: 7.5e-13, −7.50e-13 eV). Small errors fail
     (`RegressionPesScmNegative`): P17/P18 read as Å at 35244 points, θ, φ swapped once at 39499,
     1.01× the tolerance at the largest |V| there only; 0.99× at the worst point passes.
   - Gate 1c on hardwall_small, cu111_dailight_small, cu111_persson_small, crp_cu111_srp48_small and
     scm_cu111_srp48_small, two cuts each (the Fortran defaults x = y = 0, θ = 90°, φ = 0; and
     x ≈ 0.3 alat, y ≈ 0.1 blat, θ = 37°, φ = 113°; the hard wall's V depends on neither, so its two cuts
     hold the same values), 5280–9856 points per cut, each Fortran file's header (x, y, θ, φ to its 10
     decimals) checked against the cut it should be: every V within half a
     printed unit of the e15.8 value plus `potential_tolerance_of(case)` × max(|V|, 1 eV) (the printed
     half unit binds: the worst points are at 0.989–0.990 of the tolerance, i.e. the C++ agrees with all
     8 printed digits); r and Z within half a printed unit plus 1e-14 bohr (the Fortran's Z = 0 is
     4.4e-16 in the C++). Measured 2026-11-20. Negative tests (`RegressionPesCutChecks`): one Fortran
     value moved by one printed unit away from the C++ value fails at exactly that point, a NaN at exactly
     that point, and the angles handed to the PES in degrees (a porting slip) fail at most points of the
     Dai–Light cut 2. Gate 1c checks V to its 8 printed digits; the tight check is this C++-only one: a
     cut at a grid point and a DVR node (in degrees) equals the potential grid of cu111_dailight_small at
     its 5656 points with Z ≤ ZminSwitch within 1e-13 of max(|V|, 1 eV), and x and y swapped fails at
     exactly the points where the grid differs.
2. **Gas-phase rovibrational states** (same tests; Fortran `GasPhase.f90` / `FourierGridHamiltonian.f90`): the
   energy and wave function of every (v, j) level against the Fortran's full-precision dump
   `debug_asyp.dat` (`run_setup_reference.sh` sets `DoDump`): the analysed levels to 1e-12 relative
   (energy) and 1e-10 absolute (ψ), the other levels to 1e-10 relative (energy); LAPACK against
   the Fortran's Jacobi ([fortran_deviations.md](fortran_deviations.md)). Measured on hardwall_small
   (2026-11-15): 1.3e-14 and 8.7e-15. The 6-digit printed table (e.g. e(0,0) = 0.258378 eV,
   e(1,0) = 0.762842 eV) is checked by the unit test `GasPhase.ReferenceCaseSixDigitValues`.
   The other cases of gate 1 run the same comparison (`RegressionPes.Gate2RovibrationalStates`).
3. **Propagation** (`tests/test_regression_propagation.cpp`): initial wave-packet block
   (norm, |ψ(k)| values, iLow/iHigh, norm inside) and the averages table (time, norm, ⟨Z⟩, ⟨r⟩,
   rest norm) against the Fortran output (stdout; averages from averages.out), to all printed
   digits: |C++ − Fortran| ≤ half a unit in the 8th mantissa digit of E15.8 (the Fortran is not
   patched). The |ψ(k)| values (k = 0, grid edge, maximum) allow max(half a printed unit, 1e-14): at
   k = 0 or the edge they can be FFT rounding noise (user's decision 2026-11-15; looser for the edge
   value, which had a fixed 1e-15; `EXPECT_MOMENTUM_EQ`, tests/propagation_reference.hpp). Passes
   for hardwall_short (21 rows) and the full hardwall_small run (7000 steps, 141 rows; worst
   0.989645 of half a printed unit). Default case `tests/reference/hardwall_short` (same as
   hardwall_small but 1000 steps); the full 7000-step hardwall_small run only with
   `SPODVR_LONG_TESTS=1`. The hard wall has no θ/φ/x/y coupling, so gate 3 only checks that the
   angular transform and x/y FFTs act as the identity there; their correctness rests on the
   analytic unit tests (test_angular, test_propagation). Corrugated cases (60° cell, C6v):
   `cu111_dailight_short` (1000 steps, default), `cu111_dailight_short_sym` (C6v UseSymmetry,
   nPhisOptimal: nφ 12), `cu111_dailight_small` (7000 steps, long; also fine analysis, gate 5
   `Gate5DaiLight`). Status: all pass (gates 3, 4, 5). Full Dai–Light run (7000 steps): averages
   worst 0.986583 of half a printed unit, BK (1400 records) 6.41e-13 of max|c| = 0.265,
   probabilities within half a unit (worst 0.989415) except the reaction probability (see gate 4).
   The 1000-step Dai–Light cases: BK 6.25e-12 (short) and 6.288e-12 (short_sym) of the scale,
   against the tolerance 1e-10. Run times of this case: performance.md "Fortran vs C++".
   fast-math (`SPODVR_FAST_MATH`) moved the reaction
   probability to 2.0e-12 / 5.3e-12 from the Fortran (old floor 1e-12): floor raised to 1e-11;
   not the default build.
   Low-Z absorbers (upstream 5.24): `cu111_dailight_short_lowz` (default, CPU and GPU) is
   `cu111_dailight_short` plus `ZOptLowEnergy` 0.005 eV and `ZSpOptLowFactor` 4e-4 hartree, both
   from ZStart = −2 to Zana = 8. Each matters on its own (C++ runs with one removed, 2026-11-12):
   without the scattering one the averages move by 3.3e4 printed half units and the BK coefficients
   by 9.6e3 gate tolerances, without the specular one by 3.1e6 and 5.5e5 (the summed probabilities
   by only 1.2 and 2.6e3 half units: gate 4 sees the scattering one in the BK coefficients).
   Margins (CPU / GPU): averages 0.976 / 0.976, probabilities 0.989 / 0.991 of half a printed unit,
   BK 5.3e-12 / 2.3e-12 of the scale. The packet stays at Z ≳ 2 by t = 2000 (below Z = 2 then: 2e-11
   of the specular, 3e-8 of the scattering norm; C++ checkpoint, 2026-11-12), so the gate exercises
   the weak part of the profile (≤ 36 % of its peak); the strong end at Start is covered by
   `PropagationTables.LowZAbsorberProfileIsTheMirroredFormula`, the legacy eV→hartree factor of the
   strength by `Input.LowZAbsorberStrengths`. Negative test (Long suite: its outcome changes only
   with the input or the reference; a code regression fails the gate itself)
   `LowZGateFailsWithoutTheScatteringAbsorberLong`: without that absorber the C++ misses the
   Fortran by 3.3e4 half units and 9.6e3 BK tolerances.
   SCM (SystemID `CRPDistrtSurface`): `scm_cu111_srp48_small` (default, CPU and GPU; the case is
   described with its gate 1 above), 3000 steps, so that the wave reflected by the surface reaches
   Zana (before t = 3600 every BK coefficient is below 1e-4, the floor of the scale). Margins (CPU /
   GPU): averages 0.986 / 0.986, probabilities 0.988 / 0.988 of half a printed unit, BK 5.9e-13 /
   5.1e-13 of the scale 0.075. Negative test (Long suite) `ScmGateFailsWithoutTheCorrectionLong`:
   with the CRP term alone (ΔV removed, the lattice kept) the C++ misses the Fortran by 7.8e9 half
   units and 1.4e10 BK tolerances (2026-11-15).
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
   Since the reference is r1727 (2026-11-24) also every S-matrix row of the Fortran's on-the-fly
   `fortran_fort25.txt` (tests/fort25_reference.hpp; r1720 refers the initial state's phase to the specular
   analysis line, as the C++ does). The tolerances add the Fortran's printed digits and a code term
   c · max(S_scale, 1e-3), c = 1e-10, S_scale the case's largest \|S\|:

   - S: \|S_C++ − S_F\| within the printed half units of Re S and Im S (combined; a printed zero is exact) plus
     the code term;
   - P: within its printed half unit plus the code term's effect on \|S\|²;
   - the phase: within its printed half unit plus code term / \|S\| (the Fortran prints the phase of its
     unrounded S), compared where that bound is at most 1e-3 rad, i.e. \|S\| ≳ 1e3 times the code term (31624
     of the 41536 rows of cu111_dailight_short, 36534 of cu111_dailight_small, none in hardwall_short); a C++
     phase must be finite on every row;
   - each block's printed total energy.

   The code term was measured on every case with the CPU, the GPU and both fast-math builds: the largest
   difference beyond the printed digits is 5.99e-12 of max(S_scale, 1e-3), in hardwall_short with the fast-math
   CPU build (its elements are all rounding noise, \|S\| ≤ 2.4e-15, so the floor 1e-3 binds); 1.1e-12 in the
   Dai–Light cases with fast-math, at most 2.1e-13 with the default builds. c = 1e-10 is the smallest power of 10
   with a margin of 10 over 5.99e-12. Where \|S\| is large (cu111_dailight_small, scm_cu111_srp48_small) the
   printed digits bind on every row. Counts: the rows read equal the Fortran's print rule (the open channels,
   E_total > E(v, j)) applied to the C++ run and an independent count of the file's row lines; every row matches
   an open C++ state of its block's diffraction order; the number of phases compared is pinned per case. A run with
   a single diffraction order is refused: the Fortran then prints one block header per v for all rotational
   states. Negative check on the same run: the initial state's phase referred to Zana (the Fortran's convention up
   to r1719) fails exactly the rows r1720 changed (50 per Dai–Light case, 51 in hardwall_small, 21 in
   scm_cu111_srp48_small; none in hardwall_short, whose line never fixes); against the r1713 reference files the
   C++ fails exactly those rows and the Zana convention none (measured once, 2026-11-24).
   Status: all pass. Full hard-wall run (hardwall_small, 7000 steps): probabilities worst 0.989451
   of half a printed unit, BK (1400 records) 6.36e-13 of max|c| = 0.409; Fortran restart run
   worst 0.99. In hardwall_short all BK coefficients are
   rounding noise (≤ 2e-15, the packet has not reached Zana), hence the scale floor 1e-4.
5. **Fine analysis** (`tests/test_regression_fine_analysis.cpp`): the C++ analyses the Fortran's
   own BK file (hardwall_small `fortran_fort21.bin`); Fortran reference from unmodified
   onlyAnalysis runs (`tools/fortran_reference/run_fine_analysis_reference.sh`, variants full,
   truncated (tEnd 13860), small (tEnd 800), empty (tEnd 100, nen95 = 0) in
   `hardwall_small/fine_<variant>/`). The variant with NrOfEs −95 (analysis energies on fine-grid points) was
   removed with r1719, which makes −95 an input error; no Fortran-compared case has its energies on grid points
   any more, and the on-grid interpolation is covered by unit tests (`LocateOnGrid.NeighboursAndClamps`,
   `InterpolateSMatrix.GridPointsAndLinearProbability`). Compared to half a printed unit: AsympInfo numbers,
   probaf table, all of fil10 and unit 3 (spectrum), the interpolated per-energy blocks and
   fort.25 probabilities; fort.25 Re S / Im S / phase against an emulation of the reference's
   interpolation, whose ATAN and sign bugs 5.24 keeps (see fortran_deviations.md). The gate analyses
   the Fortran's BK file without a checkpoint, so every phase refers to Zana as in the Fortran
   (fortran_deviations.md, analysis D2). Status: all
   pass (worst 0.990127 of half a printed unit, in `Gate5DaiLight`).
   Independent of the Fortran: `ElasticPhase.HardWallElasticElementIsTheExactReflection` compares the
   elastic element of a hard-wall run (hardwall_small with JMax 2, i.e. 6 rotational states (JMax 0 is an input
   error since r1726), and 2 × 2 (x, y) points, 7000 steps, ≈ 20 s at 8 threads) with the exact 1D reflection amplitude (Numerov) at
   0.25–0.55 eV: phase within 1e-2 rad, \|S\| within 4e-3 (measured 6.24e-3 rad and 1.85e-3: the time-step
   error and the norm not yet across the line at tEnd); the phase referred to Zana (the Fortran's fine analysis,
   and its on-the-fly analysis up to r1719) fails at all 31 energies.
   Note: the associated Gauss-Legendre quadrature for |mJ| > 0 is dead code in the Fortran
   (`setupdvr.f90` forces `SymmetryInPhi = .FALSE.`, so only mJ = 0 is used). The rotational
   basis uses `legendre.f90 AssociatedLegendre` instead (C++ `associated_legendre`), tested
   analytically (closed forms, orthonormality, cos²θ matrix elements, transform round trip).
