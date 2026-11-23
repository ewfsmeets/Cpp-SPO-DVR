# Bugs: the highest analysed j for odd JIn, and the VIn check against VAna

**Status:** fixed upstream on 2026-11-23 in r1726 (the input checks) and r1727 (the search starts at `jin`), and
ported to the C++ (section [The C++](#the-c)). Found on 2026-11-23 while checking which input paths no reference case
covers (every Fortran-compared case has `VIn = JIn = MJIn = 0`), and verified by runs of both codes. The sections
below describe the bugs as they were in 5.19, 5.24 and upstream master up to 8c2b3fe (r1725).

**In one paragraph.** The analysis projects the outgoing wave onto the channels (v, j, mJ, diffraction order) with
v below `VAna` levels and j up to `jana`. Two of the limits are wrong:

1. **The highest analysed j for a homonuclear molecule with odd `JIn`.** When `JAna` is not given, the program
   searches the highest open j from j = 0 in steps of 2, so it only ever tests even j. The basis of an odd `JIn`
   holds only odd j, so `jana` comes out as the highest open *even* level. The highest open odd level is then left
   out of the analysis whenever it lies just above that one, and no warning is printed. With `JIn = JMax` the
   initial state itself is left out: the Fortran writes out of bounds and reports zero scattering, and the C++
   crashes.
2. **The check of `VIn` against `VAna`.** `VAna` is a number of levels (v = 0 … `VAna` − 1 are analysed), but the
   input check rejects only `VIn > VAna`. `VIn = VAna` passes, and the initial vibrational level, with every level
   above it, is then left out of the analysis.

In both cases every analysed channel is still right. What is wrong is everything that assumes the analysed channels
are complete: the scattered sum, the reaction probability 1 − sum, the averages over final states, and (bug 1)
anything indexed by the initial state.

## Contents

- [Background: which channels are analysed](#background-which-channels-are-analysed)
- [Bug 1: the highest analysed j for odd JIn](#bug-1-the-highest-analysed-j-for-odd-jin)
- [Bug 2: VIn against VAna](#bug-2-vin-against-vana)
- [Evidence](#evidence)
- [History](#history)
- [Proposed fix: Fortran](#proposed-fix-fortran)
- [The C++](#the-c)
- [Open points](#open-points)

## Background: which channels are analysed

- **The analysis basis.** The Balint-Kurti (BK) coefficients, and from them the S-matrix, are computed for the
  rotational-diffraction states of the propagation basis with j ≤ `jana`, for each of the `nviba` lowest
  vibrational levels (`RotDifBasis.f90`, `Analyzor.f90` `BalintKurtiCoeffs`; C++ `src/analysis/analysis_basis.cpp`,
  `src/analysis/balint_kurti.cpp`).
- **The limits.** `JAna` (Fortran `jana`, a highest j, inclusive) and `VAna` (`nviba`, a number of levels) are
  input keys. When they are absent (or `VAna` < 1), the program sets them to the open channels: the levels with
  E(v = 0, j) − E(0, 0) at most `EMaxWP` + (E(VIn, JIn) − E(0, 0)) + `EParallel` (`GasPhase.f90`; C++
  `src/gasphase/gas_phase.cpp`).
- **A homonuclear molecule** (`HomoNuc` true) conserves the parity of j: the propagation basis holds only the j of
  the parity of `JIn` (`angulargrid.f90` `InitAngularGrid`: `JMin` and `JMax` adjusted to that parity, `dJ` = 2; C++
  `src/grids/angular_grid.cpp:37`). From j = 1 the allowed final states are j = 1, 3, 5, …; the even j are not
  reachable and not in the basis.
- **The initial state** must be in the analysis basis. Its index `inchna` (C++ `AnalysisBasis::initial_state()`)
  labels the channel whose coefficient is taken from the specular grid, the reference of the fine analysis, and the
  initial-state tables of the output.

## Bug 1: the highest analysed j for odd JIn

### Fortran

- 5.19 `GasPhase.f90`; upstream master (8c2b3fe) line numbers in brackets.
- `SetupAsympWF`, lines 156–174 [168–186], inside the block marked "ADDED BY E. PIJPER":

  ```fortran
  DeltaJ = Setup_TheAngularGrid%dJ
  ...
    j = 0   ! Start here.
    ! totEnergy is the maximum energy available for scattering.
    totEnergy = EMaxWP/chaev + (Evinjin - E00) + EParallel
    DO
      Evj = eivaa( j, 1 )
      if ( (Evj - E00) > totEnergy ) EXIT                ! eigenValue in au.
      j = j + DeltaJ                                     ! Next j-state.
      if ( j > jmaxw ) then
        call ShowWarning( 'Rotational basis too small to describe all open channels!')
        EXIT  ! Larger j not possible.
      endif
    ENDDO
    jmaxAnalysis = j - DeltaJ
  ```

  With `DeltaJ` = 2 the search visits j = 0, 2, 4, …, whatever the parity of the basis. Line 181 [193] copies the
  result to `jana` when `JAna` is absent.
- `RotDifBasis.f90` line 572 [584] keeps the basis states with `J <= jana`; lines 661–669 [673–681] set `inchna`
  only if the initial (`jin`, `mjin`) is among them. `inchna` is a module variable without an initial value
  (`globaldata.f90` line 177).
- The warning, `GasPhase.f90` lines 267–271 [279–283], tests the level `jana + 2`:

  ```fortran
  IF( lhomo ) THEN
    IF( jana < (jmaxw - 1 ) ) THEN
      Checkenergy = eivaa( jana + 2, 1 )                 ! Check for max. J channel
      IF( ( eivaa( jin, nvin ) + ectmax/chaev ) > Checkenergy ) &
        PRINT *, '## WARNING: Not all open j channels analysed!'
  ```

  With `jana` even this is again an even level, not the next level of the basis, so the warning does not see the
  odd level that was left out.

### C++

The port reproduces it:

- `highest_open_rotational_level` (`src/gasphase/gas_phase.cpp:63`) starts at `int j = 0` and steps by
  `angular.delta_j()`.
- `AnalysisBasis` (`src/analysis/analysis_basis.cpp:97`) keeps j ≤ `jana`; `initial_` stays −1 when the initial
  state is not among them (`src/analysis/analysis_basis.cpp:104`). Nothing checks it.
- The warning (`src/gasphase/gas_phase.cpp:163`) tests `jana + delta_j`.

### What it causes

- **Who is affected:** every run with `HomoNuc` true, odd `JIn` and no `JAna`: ortho-H₂ (j = 1, 3, …) and para-D₂.
  Even `JIn` and heteronuclear molecules are not affected: there the ladder from j = 0 is the basis's own.
- **An open level is left out.** Let J_e be the highest open even level. The search returns `jana` = J_e, so the odd
  levels up to J_e − 1 are analysed. The highest open odd level is J_e − 1 or J_e + 1. In the second case, about
  every other energy window, it is allowed, open and propagated, but never analysed.
  - Correct: P of every analysed channel (each BK coefficient is a projection of its own).
  - Wrong: the scattered sum is too low by the probability of the missing level (summed over mJ and diffraction);
    the reaction probability 1 − sum is too high by the same amount; the averages over final states
    (`av erot`, `av edif`, …) are normalised by the incomplete sum.
  - From j = 1, the lost level is the Δj = 2 excitation j = 3 for `EMaxWP` + `EParallel` between 0.071 and
    0.126 eV (the levels of example 1 below: j = 3 open, j = 4 closed); at higher energies it is a higher odd level.
  - It cannot be repaired afterwards: the BK file holds only the analysed channels, so an analysis-only run of the
    same file has the same gap.
  - No warning is printed.
- **`JIn = JMax` (odd): the initial state is not analysed.** The even ladder can never reach `JMax`, so `jana` =
  `JMax` − 1 < `JIn` at any energy (only the "Rotational basis too small" warning is printed). The same happens at
  low energy when E(0, `JIn` + 1) − E(`VIn`, `JIn`) > `EMaxWP` + `EParallel` (the rotational gap of the next even
  level).
  - Fortran: `inchna` keeps its unset value (printed as 0). Once the specular analysis line is fixed,
    `BalintKurtiCoeffs` writes `balcf(inchna,iv)` (`Analyzor.f90` line 388 [400]) outside the array, every analysis
    step. The run completes with exit status 0 and reports a scattered sum of 0 and a reaction probability of 1.
  - C++: `BalintKurtiAnalyzer` writes `bk(-1, v)` (`src/analysis/balint_kurti.cpp:80`), and the run ends in a
    segmentation fault at that step.
  - The fine analysis reads the initial channel's energy `enymo(inchna,nvin)` (`AsympInfo.f90` line 220 [232]; C++
    `src/analysis/fine_analysis.cpp:158`) at index 0 / −1 (from the code; not run).
- **Workaround:** give `JAna` (e.g. `JAna = JMax`).

## Bug 2: VIn against VAna

### Fortran

- `VAna` is read into `nviba` (`InputNew.f90` line 264 [276]) and is a **number of levels**:
  - every loop over levels is `do iv = 1 , nviba` with iv = v + 1 (`Analyzor.f90` lines 319, 383, 400, 456);
  - the program's own value is the number of open levels (`GasPhase.f90` lines 122–147 [134–159]: the loop counts
    iv = v + 1 up to the first closed level, `ivmaxAnalysis = iv - 1`).
- The input check, `ValidateInput`, line 641 [713]:

  ```fortran
  IF ( getIntegerField('VAna',dummy_int) ) CALL ERROR( vin > nviba, 'VIn > VAna.' )
  IF ( getIntegerField('JAna',dummy_int) ) CALL ERROR( jin > jana, 'JIn > JAna.' )
  ```

  `jana` is a highest j, so `jin > jana` is right. `nviba` is a count, so the matching check is `vin >= nviba`:
  `VIn = VAna` passes, and v = `VIn` is then not analysed. (The output labels `VAna` "Maximum V analysed", line 1097
  [1171], which suggests a highest v.)
- A side effect of the same line: `VAna = 0` means "set by the program" (`GasPhase.f90` line 145 [157]), but the
  check runs whenever the key is present, so `VAna = 0` with `VIn` ≥ 1 is rejected.

### C++

The port reproduces the check (`src/input/validate.cpp:183`: `p.v_in > *p.v_analysis`). The analysis uses
`AnalysisBasis::n_v()` ("Fortran: nviba", a count) in the same way.

### What it causes

With `VIn = VAna` (`VAna` given):

- **Right:** the S-matrix and P of the analysed levels v < `VIn` (each is an independent projection, and the
  energies use E(`VIn`, `JIn`), which is set regardless).
- **Wrong:** the scattered sum lacks the vibrationally elastic level and every level above it, usually most of the
  scattered flux; the reaction probability 1 − sum and the averages over final states with it. The Fortran and the
  C++ print "Not all open v channels analysed!", but accept the input.
- **Out of range:**
  - C++, every on-the-fly run: the screen table "Diffraction orders of the initial state" reads `p_order` at
    v = `VIn` (`src/report/analysis_report.cpp:115`), beyond its `n_orders · n_rot · n_v` elements
    (`src/analysis/s_matrix.cpp:135`).
  - Fine analysis (`onlyAnalysis` / `alsoAnalysis`), both codes: the energy of the initial channel,
    `ev0jl = enymo(inchna,nvin)` (`AsympInfo.f90` line 220 [232]; C++ `src/analysis/fine_analysis.cpp:158`,
    `src/app/fine_analysis_run.cpp:78`), is read beyond the analysed levels. It is the reference of the
    time-to-energy transform and of the open-channel test, so the fine analysis would be wrong for every channel
    (from the code; not run).

## Evidence

Runs of the Fortran 5.24 r1713 (`-O3`, the C++ port's reference) and of upstream master 8c2b3fe (built the same way,
`GNU.mk` `O3FLAGS = -O3`), and of the C++ at commit 000b505. Each input is a reference case of `tests/reference/`
with the keys listed changed (`DVR.inp` and `input.json` alike); the executables run in a fresh directory, 4 threads.

### 1. Odd JIn: the open level j = 3 is not analysed

- Input: `hardwall_short` with `JMax` 5, `MJMax` 5, `JIn` 1, `MJIn` 0, `EMinWP` 0.06, `EMaxWP` 0.12, and
  `Potential` true (set-up only).
- Gas-phase levels (both codes, equal to the printed digits), E(0, j) − E(0, 0):

  | j | 1 | 2 | 3 | 4 | 5 |
  | --- | --- | --- | --- | --- | --- |
  | eV | 0.01424 | 0.04260 | 0.08482 | 0.14052 | 0.20923 |

- Available: 0.12 + 0.01424 = 0.13424 eV, so j = 1 and 3 are open and j = 5 is closed (the basis holds 1, 3, 5).
- Both codes: highest analysed j 2. The Fortran (r1713 and master) prints `nana = 3`: j = 1 with its three mJ.
  Correct is j ≤ 3, `nana` = 3 + 7 = 10. No warning.

### 2. JIn = JMax: the initial state is not analysed

- Input: `cu111_dailight_short` with `JMax` 3, `MJMax` 3, `JIn` 3, `MJIn` 0 (the specular analysis line fixes at
  Z = 10.88 bohr, t = 1610 a.u.).
- Fortran master: `nana = 3`, `inchna = 0`, the warning "Rotational basis too small to describe all open
  channels!"; the run ends with exit status 0 and prints, at all 41 energies, a sum of 0.0000000 and a reaction
  probability of 0.1000000D+01.
- C++: highest analysed j 2, `nana` 3; "Segmentation fault" (exit status 139) right after "analysis surface on the
  specular grid fixed at Z = 10.880000 bohr (t = 1610.00 a.u.)".

### 3. VIn = VAna: the initial vibrational level is not analysed

- Input: `cu111_dailight_short` with `VIn` 1, once with `VAna` 1 and once without `VAna` (the program then analyses
  the 3 open levels).
- A short case (2000 a.u.): only part of the reflected packet has crossed the analysis line, so the probabilities are
  small; what matters is the comparison.
- C++, `scattering_summary.tsv`:

  | E_coll (eV) | P(v = 0), `VAna` 1 | P(v = 0), 3 levels | P_scattered, `VAna` 1 | P_scattered, 3 levels |
  | --- | --- | --- | --- | --- |
  | 0.2 | 1.3617721144711988e-08 | 1.3617721144711988e-08 | 1.36e-08 | 9.96e-07 |
  | 0.4 | 7.2191437757944002e-09 | 7.2191437757944002e-09 | 7.22e-09 | 1.43e-06 |
  | 0.6 | 4.0677863893268805e-08 | 4.0677863893268805e-08 | 4.07e-08 | 1.79e-05 |

  The analysed level is bitwise the same; the sum lacks 99.8 % of the scattered flux at 0.6 eV. The table
  "Diffraction orders of the initial state (v=1, …)" is printed for a level that was not analysed.
- Fortran, `VAna` 1 (r1713 and master, identical): `nviba = 1`, "Not all open v channels analysed!", and at all
  41 energies a sum of 0.0000000 and a reaction probability of 0.1000000D+01. Without `VAna` (master): `nviba = 3`
  and a sum of 0.0000010, 0.0000014 and 0.0000179 at 0.2, 0.4 and 0.6 eV, as the C++ above.

## History

- The open-level search is in the block of `GasPhase.f90` marked "ADDED BY E. PIJPER" (lines 107–184 in 5.19). It
  reads `DeltaJ` from the angular grid, so the parity was considered for the step but not for the start.
- Both bugs are unchanged in upstream master 8c2b3fe, which includes the input-check fixes of 19–22 November 2026.
- The C++ port reproduced both. No reference case has an odd `JIn` or `VAna` given, so no gate could see them.

## Proposed fix: Fortran

**Bug 1.** Start the search at the first j of the basis, which has the parity of `jin` (`angulargrid.f90` lines
116–121):

```fortran
    IF ( thetaDoF .AND. phiDoF ) THEN
      j = Setup_TheAngularGrid%JMin   ! Start here: the lowest j of the basis (parity of jin if lhomo).
```

- E(0, `JMin`) ≤ E(0, `jin`) ≤ E(`vin`, `jin`) is always open, and E(0, j) increases with j, so the result is
  ≥ `jin`. An odd `JMax` becomes reachable, and the warning's `jana + 2` is then the next level of the basis.
- Even `jin` and `lhomo` false give `JMin` = 0: unchanged.
- In addition, stop when the initial state is not in the analysis basis (a `JAna` < `JIn` is already rejected; this
  catches every other way): set `inchna = 0` before the loop in `RotDifBasis.f90` and add after it
  `CALL ERROR( inchna < 1, 'The initial state is not in the analysis basis (jana < jin).' )`.

**Bug 2.** Compare with the number of levels, and only when `VAna` is given as a number of levels:

```fortran
    IF ( getIntegerField('VAna',dummy_int) .AND. nviba >= 1 ) &
       CALL ERROR( vin >= nviba, 'VIn >= VAna (VAna is the number of vibrational levels analysed).' )
```

**Tests:** example 1 must give `nana = 10`; example 2 must analyse j = 1 and 3 (`inchna` ≠ 0) and report a non-zero
sum; example 3 with `VAna` 1 must stop with the error, and `VAna` 2 must run.

## The upstream fix

- **r1727** (`GasPhase.f90 SetupAsympWF`): the search starts at `j = jin` instead of 0. Like the proposal above, it
  then tests only levels of the basis's parity, and its result is ≥ `jin`.
- **r1726** (`InputNew.f90 ValidateInput`): `JIn >= JMax`, `VIn >= VAna` and `JIn >= JAna` are input errors (were
  `>`). This differs from the proposal in three ways:
  - `JIn = JMax` is rejected, so example 2 is now an input error, and so is every input with JMax = 0;
  - a given `VAna = 0` ("set by the program", `GasPhase.f90`) is rejected whatever `VIn` is; leaving `VAna` out keeps
    the program's choice;
  - `JIn = JAna` is rejected, although `jana` is an inclusive highest j (`RotDifBasis.f90` keeps `j <= jana`), so a
    valid input (e.g. JIn = JAna = 0, only j = 0 analysed) is refused.
- Not changed upstream: the check that the initial state is in the analysis basis (`inchna`), and the case
  VIn = nR − 1 without `VAna`, where at most nR − 1 levels are analysed (bug 2 remains there; read from the code).

## The C++

The C++ follows upstream in all three r1726 checks (the user's decision, 2026-11-23) and ports r1727:

- `highest_open_rotational_level` (`src/gasphase/gas_phase.cpp`) starts at JIn;
- `validate_input` (`src/input/validate.cpp`): `JIn >= JMax`, `VIn >= VAna`, `JIn >= JAna`;
- [analysis D3](fortran_deviations.md#analysis-d3): VIn = nR − 1 is an input error;
- `compute_gas_phase` and `AnalysisBasis` stop when the initial state is not analysed (internal checks, unreachable
  for a valid input). This also covers the sites the C++ reached out of range in bug 1 with JIn = JMax and in bug 2:
  `balint_kurti.cpp`, `fine_analysis.cpp`, `fine_analysis_run.cpp`, `analysis_report.cpp` (the initial state's
  diffraction table and its set-up label), `fine_analysis_report.cpp` and the `initial_state_index` of analysis.json.
- Tests: `GasPhase.OddJInSearchReachesJ3` (example 1), `GasPhase.JMax7OpensJ5`, `GasPhase.OddJInWarnsForTheNextOddLevel`,
  `GasPhase.PostConditionsThrow`, `AnalysisBasis.OddJInExampleHasTenRotationalStates`,
  `AnalysisBasis.InitialStateOutsideTheBasisThrows`, the input rows `JInEqJMax`, `VInEqVAna`, `VAnaZeroGiven`,
  `JInEqJAna`, `AnalysisD3VInIsTheLastLevel`, and the designer's copy (`tools/wave_packet_designer`).
- BK and checkpoint files written before the fix by a run with odd JIn and without JAna have another number of
  analysed states; reading them now stops with exit code 4 (give the old run's JAna to read them).

## Open points

- **Reported upstream** and fixed on 2026-11-23 (r1726/r1727). Still to tell the authors: `JIn >= JAna` refuses the
  valid JIn = JAna, and VIn = nR − 1 (without VAna) leaves the initial level out of the analysis.
- **The size in production use** is estimated from the level ladder, not measured: which odd level is lost depends
  on the energy window; the loss is largest for the Δj = 2 excitation at low energy.
- **The fine-analysis consequences** (bug 1 with `JIn = JMax`, bug 2) are read from the code, not run.
