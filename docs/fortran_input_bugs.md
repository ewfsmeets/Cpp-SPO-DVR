# Fortran input checks: bugs and how the C++ handles them

The C++ program reproduces the input checks of the Fortran code (`InputNew.f90` ReadInput and
ValidateInput, and the grid checks of `Setup.f90`; see [input_schema.md](input_schema.md),
"Validation"). A few of those checks are wrong. This document describes each one:

- the Fortran code;
- what goes wrong;
- a demonstration with the unmodified Fortran executables;
- what the C++ does instead.

The short entries are in [fortran_deviations.md](fortran_deviations.md), under the same IDs.
The sections cover the flux analysis checks (flux Q1, flux Q2), the default strength that the
Fortran keeps for an absorber energy of 0 (absorber D1), the optical potentials at the low end
of the Z grid that upstream 5.24 added (low-Z D1–D4), the number of analysis energies
(analysis D1), and dtAna or dtAve just below a multiple of dtStep or dtAna (analysis D4, open in the Fortran).

**Status (2026-11-23):** upstream fixed these checks in r1714–r1719 (flux Q1 in r1714, absorber D1 in r1718,
analysis D1 in r1715–r1719, flux Q2 in r1719; low-Z D2 leaves the absorber off since r1717 instead of producing
NaN), and the C++ follows the fixed checks. The reference of the port is r1727 since 2026-11-23, so it has them. Each section starts with its status. The demonstrations below show the
Fortran before those fixes (5.19 and r1713).

## The demonstration runs

- **Executables:**
  - 5.19, the version the port started from: the archive install `v5.19`, whose `bin/spodvr.x` is
    the 5.19 distribution's `EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x` (gfortran 15,
    `-O3`, OpenMP, FFTW);
  - 5.24: the install `r1713-dce0ecc`, upstream master at commit dce0ecc (svn r1713), built with
    `O3FLAGS = -O3` in `GNU.mk` (the default flags enable fast-math); the reference of the port until
    2026-11-23. The reference since then, `r1727-3d00331`, has the fixes ([validation.md](validation.md#reference-data)).
  - The installs: [tools/fortran_reference/README.md](../tools/fortran_reference/README.md#the-reference-install).
- **Input:** `tests/reference/hardwall_small/DVR.inp` with the keys of each case appended (or, where
  a case says so, lines of it replaced).
  - Its wave packet spans EMinWP = 0.2 to EMaxWP = 0.6 eV.
  - rFluxAnalysis = 3.05 bohr is point 16 of the r grid (0.8 + 15 × 0.15) and lies below
    ROptStart = 4.55, so the checks on r pass.
- **How:** each run in a new scratch directory, 4 OpenMP threads. The runs that `timeout` stops
  have unbuffered stdout (`GFORTRAN_UNBUFFERED_PRECONNECTED=y`), so they keep all of their output.
- **Exit status:** an input error ends the Fortran program with a plain `STOP`
  (`errortrap.f90 AbortWithError`; `__MPI_STOP` without MPI). Its exit status is therefore 0, the
  same as for a successful run. Status 124 below means that `timeout` stopped the run.
- 5.19 and 5.24 print the same messages in every case below: the checks and
  `FluxAnalyzor.f90 setupFluxAnalyzor` are identical in both versions.

## flux Q1: EMaxFluxAna is checked the wrong way round

**Status:** fixed upstream in r1714 (2026-11-19): the Fortran checks EMinFluxAna > EMaxWP, EMaxFluxAna < EMinWP and
EMaxFluxAna > EMaxWP, and the wrong check is gone. The C++ has the same checks.

### The Fortran code

`InputNew.f90 ValidateInput`, lines 603–613 in 5.19 and 673–683 in 5.24 (identical):

```fortran
    if ( doFluxAnalysis ) then
      CALL ERROR( rFluxAnalysis .RLT. rmin, 'rFluxAnalysis < rmin' )
      CALL ERROR( rFluxAnalysis .RGT. (rmin+(nr-1)*delr), 'rFluxAnalysis > rmax' )
      Test = (NINT((rFluxAnalysis-rmin)/delr)*delr) .RNE. (rFluxAnalysis-rmin)
      CALL ERROR( Test, 'rFluxAnalysis not on a grid point.' )
      CALL ERROR( EMinFluxAna .RLT. EMinWP, 'EMinFluxAna < EMinWP' )
      CALL ERROR( EMaxFluxAna .RLT. EMaxWP, 'EMaxFluxAna < EMinWP' )
      CALL ERROR( EMinFluxAna .RGE. EMaxFluxAna, 'EMinFluxAna >= EMaxFluxAna' )
      CALL ERROR( deltaEFluxAna < 0, 'deltaEFluxAna < 0' )
      CALL ERROR( deltaEFluxAna .RGT. (EMaxFluxAna-EMinFluxAna), 'deltaEFluxAna > EMaxFluxAna-EMinFluxAna' )
      CALL ERROR( rFluxAnalysis .RGT. rsplmi, 'rFluxAnalysis > ROptStart' )
```

`ERROR(Test, Message)` prints `ERROR: <Message>` and stops when `Test` is true. `.RLT.`, `.RGT.`
and the other dotted operators compare with an absolute tolerance of 1e-12 (`MyOperators.f90`).

### What goes wrong

- The flux analysis computes reaction probabilities on the energy grid EMinFluxAna,
  EMinFluxAna + ΔE, … up to EMaxFluxAna (`setupFluxAnalyzor`). The two energy checks should keep
  this grid inside the energy range of the wave packet: EMinFluxAna ≥ EMinWP and
  EMaxFluxAna ≤ EMaxWP.
- The second check is the wrong way round: it stops when EMaxFluxAna < EMaxWP. Its message also
  names EMinWP, although the test compares with EMaxWP.
- So every range that ends inside the wave packet is rejected, unless it ends at EMaxWP (to
  within 1e-12 eV). Every range that ends above EMaxWP, beyond the energies the wave packet was
  built for, is accepted.

### Demonstration

**Case (a), a range inside the wave packet.** Appended to the base `DVR.inp`:

```text
doFluxAnalysis : LOGICAL : TRUE;
rFluxAnalysis : REAL : 3.05;
EMinFluxAna : REAL : 0.3;
EMaxFluxAna : REAL : 0.5;
DeltaEFluxAna : REAL : 0.01;
```

Both versions stop right after the echo of the input, with exit status 0:

```text
* Validating the input...

ERROR: EMaxFluxAna < EMinWP
```

With `Potential : LOGICAL : TRUE;` added as well (a potential-only run, nread = 6), the result is
the same: ValidateInput runs in every run mode.

**Case (b), a range that ends above the wave packet:** `EMaxFluxAna : REAL : 0.8;` instead of 0.5.
Both versions accept it:

- the potential-only run finishes with exit status 0. Its set-up prints
  `Flux analysis in r performed at ir =   16 (r =  0.30500000E+01)`;
- the propagating run sets up the flux analysis for 51 energies (0.3 to 0.8 eV in steps of 0.01)
  and propagates, until `timeout` stops it after 60 s (status 124):
  ```text
                 INFO FROM SETTING UP THE FLUX ANALYZOR
  ********************************************************************************

   Computing flux probabilty for   51 energies
   Flux analysis performed at ir =   16 (r =  0.30500000E+01)
   Memory required for flux analysis =  44.26 MB

   ----- Flux analyzor set up -----
  ```

### The C++ behaviour

- `src/input/validate.cpp` (`validate_input`) checks the range the right way round: it stops
  when EMaxFluxAna > EMaxWP. The test is tolerant (`.RGT.`), the mirror image of the
  EMinFluxAna check.
- Case (b) is then an input error (exit code 2):
  ```text
  Input error: EMaxFluxAna: EMaxFluxAna > EMaxWP (0.8 eV > 0.6 eV): the flux analysis energies must lie inside the wave packet's energy range; use EMaxFluxAna <= EMaxWP
  ```
- Case (a) passes the input checks. The flux analysis itself is not ported
  ([flux_analysis_findings.md](flux_analysis_findings.md)), so a propagating run then stops with
  exit code 3:
  ```text
  Not supported: doFluxAnalysis = true: the flux analysis is not implemented yet
  ```
  Analysis-only and potential-only runs accept the flux keys and ignore them.
- The two range checks r1714 added (EMinFluxAna > EMaxWP, EMaxFluxAna < EMinWP; tolerant) are in
  `validate_input` too.
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/FluxQ1EMax`, `FluxEMinAboveEMaxWP`,
  `FluxEMaxBelowEMinWP` and their tolerance rows, `Input.FluxQ1EnergyRangeInsideTheWavePacketIsValid`) and
  `tests/test_cli.cpp` (`FluxAnalysisGuard.*`).

## flux Q2: a zero DeltaEFluxAna passes the check

**Status:** fixed upstream in r1719 (2026-11-19): the Fortran rejects `deltaEFluxAna <= 0.0` (exact). The C++
check is exact now too.

### The Fortran code

In ValidateInput (quoted above; line 611 in 5.19, 681 in 5.24):

```fortran
      CALL ERROR( deltaEFluxAna < 0, 'deltaEFluxAna < 0' )
```

This is a plain `<`, so ΔE = 0 passes; the next check (`deltaEFluxAna .RGT.
(EMaxFluxAna-EMinFluxAna)`) passes too. `FluxAnalyzor.f90 setupFluxAnalyzor` then divides by
ΔE (lines 104–124 in 5.19, 116–136 in 5.24; identical):

```fortran
    delta_E_au = delta_E_eV * Constants_eV2Ha
    Escat_min_au = Escat_min_eV * Constants_eV2Ha
    Escat_max_au = Escat_max_eV * Constants_eV2Ha

    ! Number of energies for which flux results are computed.
    this%numScatEnergies = INT( (Escat_max_au - Escat_min_au)/delta_E_au ) + 1

    ! Allocate memory for scattering energies.
    ALLOCATE( this%ZScatEnergy(this%numScatEnergies) )
    ...
    Allocate( this%phi(nredist,this%nzCutoff,this%numScatEnergies) )
```

### What goes wrong

- EMinFluxAna < EMaxFluxAna is checked, so the quotient is a positive number divided by zero:
  +∞ (the build does not trap floating-point exceptions).
- Fortran leaves INT(+∞) undefined. gfortran converts it to −2147483648 (−2³¹), so the energy
  count becomes −2147483647. A stand-alone program with the same compiler (gfortran 15.2) and
  flags prints exactly that count, `****` for it in an `i4` field, and size 0 for the array.
- An array with a negative extent has size zero. Every flux array is therefore empty, every loop
  over the energies does nothing, and the run continues without an error or a warning. At the
  end, the flux result files contain no energies.

### Demonstration

**Case (c):** case (b) with `DeltaEFluxAna : REAL : 0.0;`. Both versions pass the input checks and
set up the flux analysis with the garbage count:

```text
 Computing flux probabilty for **** energies
 Flux analysis performed at ir =   16 (r =  0.30500000E+01)
 Memory required for flux analysis =   0.00 KB
```

They then propagate normally until `timeout` stops them after 90 s (status 124).

A short variant runs to the end: the same input with `tEnd : REAL : 200.0;` and
`nStep : INTEGER : 100;` (instead of 14000.0 and 7000). Both versions finish with exit status 0
and no error or warning about the flux analysis. Their flux output files:

- `flux_reaction.dat`, `flux_coeffs.dat` and `flux_coeffs_deriv.dat`: empty (0 bytes);
- `flux_density.dat`: only the header, with `*****` for the number of energies (`I5`):
  ```text
  *****
   0.15707963E+01  C2v
      4 0.15599150E+01
      4 0.13099150E+01
  ```

### The C++ behaviour

- `validate_input` stops when DeltaEFluxAna ≤ 0. The test is exact, as the Fortran's since r1719 (before, it
  was tolerant, like the "≤ 0" checks of the grid spacings dZ and dR).
- The error names the key (exit code 2). For case (c) with a range inside the wave packet
  (EMaxFluxAna = 0.5, so that the flux Q1 check passes first):
  ```text
  Input error: DeltaEFluxAna: DeltaEFluxAna <= 0 (0 eV): the flux analysis energy grid needs a positive spacing; use DeltaEFluxAna > 0
  ```
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/FluxQ2DeltaZero` and
  `FluxQ2DeltaNegative`, `Input.FluxDeltaEJustAboveZeroIsValid`).

## absorber D1: an absorber energy of 0 keeps the default strength

**Status:** fixed upstream in r1718 (2026-11-19): `ZOptEnergy`, `ROptEnergy`, `ZSpOptEnergy` and `ZSpOptLowEnergy`
≤ 0 are input errors (exact `<=`), in the Fortran and the C++. `ZOptLowEnergy` keeps the tolerant `.RLT. 0`, so 0
stays valid there (strength 0, no absorption, in both codes).

### The Fortran code

Each absorbing potential has an energy key and a factor key. `InputNew.f90 ReadInput` reads the
factor only when the energy is absent, and `ProcessInput` computes the strength from the energy
only when it is positive:

```fortran
    IF( ZOptEnergy > 0.0 ) vabszb = getACoeffOptPotential( ZOptEnergy/chaev, zsplma-zsplmi, totmass )
    IF( ROptEnergy > 0.0 ) vabsrb = getACoeffOptPotential( ROptEnergy/chaev, rsplma-rsplmi, rdmass )
    IF( ZSpOptEnergy > 0.0 ) vabssb = getACoeffOptPotential( ZSpOptEnergy/chaev, ssplma-ssplmi, totmass )
```

The strengths start with the defaults of `globaldata.f90 SetupGlobalDataModule`:

```fortran
           ZOptEnergy = 0.05                ! } Belong together
           vabszb = 4.380088200912914E-003  ! }
           ...
           ROptEnergy = 0.10                ! } Belong together
           vabsrb = 9.554538811058970E-003  ! }
           ...
           ZSpOptEnergy = 0.15              ! } Belong together
           vabssb = 9.120042351381136E-003  ! }
```

### What goes wrong

An energy of 0 passes the check `ZOptEnergy .RLT. 0.0`, so the input is valid. The strength is
then not computed, and the absorber keeps the default strength, which belongs to another energy
and another absorber length. One would expect an energy of 0 to switch the absorption off; instead
the Fortran absorbs with the default strength and prints it as a coefficient. The low-Z absorbers
(below) do not have this problem: their default strength is 0.

### Demonstration

**Case (d):** `ZOptEnergy : REAL : 0.0;` instead of 0.3, and `Potential : LOGICAL : TRUE;`
appended. 5.24 runs to the end (exit status 0) and prints the default strength:

```text
* Optical potential in Z absorbs with coefficient:          0.43800882E-02
```

### The C++ behaviour

- Since r1718, as the Fortran: an energy ≤ 0 of `ZOptEnergy`, `ROptEnergy`, `ZSpOptEnergy` or `ZSpOptLowEnergy` is
  an input error (exit code 2). For case (d):
  ```text
  Input error: ZOptEnergy: ZOptEnergy <= 0.
  ```
- A `ZOptLowEnergy` of 0 means no absorption: the strength is 0 and the low-end factors are exactly 1
  (`src/input/process.cpp`, `process_input`).
- Before r1718 an energy of 0 gave strength 0 for every absorber; for case (d) the set-up printed
  `CAP strengths (Z, r, Zsp) ..... 0.00000000e+00, 2.38963250e-02, 1.82989643e-02 hartree`.
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/…EnergyZero` and `…EnergyJustBelowZero`,
  `Input.AbsorberEnergiesJustAboveZeroAreValid`) and `tests/test_propagation_tables.cpp`
  (`PropagationTables.LowZEnergyZeroGivesNoAbsorption`).

## The low-Z absorbers of upstream 5.24

**Status:** since r1717 (2026-11-19) the Fortran switches an absorber on only for Start < End, so Start = End
(low-Z D2) leaves it off without a message instead of producing NaN; since r1718 `ZSpOptLowEnergy` ≤ 0 is an error
(absorber D1). Low-Z D1, D3 and D4 are unchanged. The C++ keeps its errors for D1–D3.

### The Fortran code

Upstream 5.24 adds an optical potential at the low end of each Z grid, for amplitude that reaches
the start of the grid (keys `ZOptLowStart`, `ZOptLowEnd`, `ZOptLowEnergy`, `ZOptLowFactor`, and
`ZSpOptLow…` for the specular grid). `InputNew.f90 ReadInput` reads the scattering grid's keys
like this (the specular block is the same with `ssplmilow`, `ssplmalow`):

```fortran
    IF( GetRealField( 'ZOptLowStart', zsplmilow ) .AND. GetRealField( 'ZOptLowEnd', zsplmalow ) .AND. zsplmilow <= zsplmalow ) THEN
        CALL ERROR( zsplmilow < zmin, 'ZOptLowStart was set smaller than ZStart' )
        CALL ERROR( zsplmalow > ZAnalysis, 'ZOptLowEnd was set higher than Zana' )
        IF( .NOT. GetRealField( 'ZOptLowEnergy', ZOptLowEnergy ) ) THEN
             CALL ERROR(.NOT. GetRealField('ZOptLowFactor',vabszblow),'ZOptLowFactor or ZOptLowEnergy should be present.')
             ZOptLowEnergy = -1.0
        ELSE
             CALL ERROR( ZOptLowEnergy .RLT. 0.0, 'ZOptLowEnergy < 0.' )
        END IF
    END IF
```

The defaults (`globaldata.f90`) are `zsplmilow = -1000.0` and `zsplmalow = zsplmilow - 1000.0`,
so the absorber is off unless both keys are read. `ProcessInput` computes the strength only for
a positive length:

```fortran
    IF( zsplmalow > zsplmilow .AND. ZOptLowEnergy > 0.0 ) vabszblow = getACoeffOptPotential( ZOptLowEnergy/chaev, zsplmalow-zsplmilow, totmass )
```

`Absorption.f90 SetupAbsorpPot` then sets the potential, the high-end profile mirrored:

```fortran
      IF( zsplmilow <= zsplmalow ) THEN        ! only if a range was given we add an optical potential at low Z
              abslen = zsplmalow - zsplmilow
              do iz = 1 , nz
                 if ( zgrid(iz) >= zsplmilow .AND. zgrid(iz) <= zsplmalow ) then
                   vimag(iz) = - ima * vabszblow * 0.5 * dfloat(nabpow + 1) *  &
                              ( (zsplmalow-zgrid(iz)) / abslen )**nabpow
                 endif
              enddo
      END IF
```

`PrintInput` prints the low-end blocks:

```fortran
    IF( ssplmalow > ssplmilow ) THEN
       PRINT '(A,T60)',       '* Using a quadratic optical potential on the specular grid in Z for low end:'
       PRINT '(A,T60,E15.8)', '* Optical potential in low end Z start at (bohr): ',zsplmilow
       PRINT '(A,T60,E15.8)', '* Optical potential in low end Z ends at (bohr): ',zsplmalow
       IF( ZOptEnergy > 0.0 ) THEN
        PRINT '(A,T60,E15.8)', '* Optical potential in Z low end absorbs optimally for (eV): ',ZOptLowEnergy
       ELSE
        PRINT '(A,T60,E15.8)', '* Optical potential in Z low end absorbs with coefficient: ',vabszblow
       ENDIF
```

### What goes wrong

- **low-Z D1:** with Start > End the condition `zsplmilow <= zsplmalow` is false: the absorber
  stays off, and no message says so.
- **low-Z D2** (up to r1713): with Start = End the absorber is switched on with length 0. The strength is not
  computed (the length is not positive), so it is the factor or the default 0. A grid point equal
  to Start gets `0 * (0/0)**2` = NaN in `vimag`, so its absorption factor `exp(-ima*vimag*dtlan)`
  is NaN too; the first time step spreads the NaN over the whole wave function (also where it is 0,
  since 0·NaN = NaN).
- **low-Z D3:** with only one of Start and End given, the absorber stays off, and so does an
  energy or factor without Start and End. With only End given, the off state depends on the
  compiler: the gfortran build does not evaluate `GetRealField( 'ZOptLowEnd', … )` once the first
  operand of the `.AND.` is false (case (g)). Fortran does not prescribe this; from reading the
  code, a build that evaluates it would switch on an absorber from the default Start = −1000 to
  End with strength 0, which overwrites the high-end absorber with 0 where they overlap.
- **low-Z D4:** both low-end blocks of `PrintInput` choose between energy and coefficient by the
  high-end `ZOptEnergy`, and the specular block prints the scattering grid's variables
  (`zsplmilow`, `zsplmalow`, `ZOptLowEnergy`, `vabszblow`) instead of its own.
- Start ≥ ZStart and End ≤ Zana are compared with the plain `<` and `>`, not with the 1e-12
  tolerance of `.RLT.` and `.RGT.` that the other checks use. This is not a bug, but the C++ must
  compare exactly as well, or the two programs accept different inputs.

### Demonstration

Every case appends its keys (typed, e.g. `ZOptLowStart : REAL : 4.0;`) to the base input. Cases
(e), (g) and (h) also append `Potential : LOGICAL : TRUE;` (a potential-only run). All of them end
with exit status 0 and print no error or warning.

**Case (e), low-Z D1:** `ZOptLowStart` 4.0, `ZOptLowEnd` 3.0, `ZOptLowEnergy` 0.3. The output
contains no line about a low-end absorber: neither `PrintInput`'s "for low end" block nor the
index lines of `Setup.f90`.

**Case (f), low-Z D2:** `ZSpOptLowStart` = `ZSpOptLowEnd` = −2.0 (ZStart, a grid point in both
programs), `ZSpOptLowFactor` 0.01, and a short propagation (`tEnd` 200.0, `nStep` 100). The
set-up of the absorbers prints the coefficient

```text
* Optical potential in Z at lower end of specular grid absorbs with coefficient 0.10000000E-01
```

and `averages.out` becomes NaN after the first analysis interval:

```text
      time             norm           z              r        rest norm

  0.00000000E+00  0.10000000E+01  0.11440000E+02  0.14708348E+01  0.10000000E+01

  0.10000000E+03             NaN             NaN             NaN  0.00000000E+00

  0.20000000E+03             NaN             NaN             NaN  0.00000000E+00
```

The same input with the scattering grid's keys (`ZOptLowStart` = `ZOptLowEnd` = −2.0,
`ZOptLowEnergy` 0.3) prints the coefficient 0 (`0.00000000E+00`, the energy is not used) and gives
the same NaN averages.

**Case (g), low-Z D3:** `ZOptLowStart` −2.0 with `ZOptLowEnergy` 0.3 and no `ZOptLowEnd`, and
the other way round (`ZOptLowEnd` 10.0, no `ZOptLowStart`): as in case (e), the output does not
mention a low-end absorber, not even the index lines that `Setup.f90` prints whenever
`zsplmilow <= zsplmalow`.

**Case (h), low-Z D4:** the specular grid's absorber alone: `ZSpOptLowStart` −2.0,
`ZSpOptLowEnd` 6.0, `ZSpOptLowEnergy` 0.3. The printout shows the scattering grid's defaults:

```text
* Using a quadratic optical potential on the specular grid in Z for low end:
* Optical potential in low end Z start at (bohr):          -0.10000000E+04
* Optical potential in low end Z ends at (bohr):           -0.20000000E+04
* Optical potential in Z low end absorbs optimally for (eV)-0.10000000E+01
```

while `Setup.f90` prints the right index range (points 1 to 101, Z = −2.0 to 6.0).

### The C++ behaviour

- `src/input/json_input.cpp` (`require_low_z_absorber`) makes the silent cases errors (exit code 2)
  and keeps the Fortran's checks, with the same exact and tolerant comparisons. For cases (e),
  (f) and (g):
  ```text
  Input error: ZOptLowStart: ZOptLowStart > ZOptLowEnd (4 > 3 bohr): the absorber would be empty; use ZOptLowStart < ZOptLowEnd
  Input error: ZSpOptLowStart: ZSpOptLowStart = ZSpOptLowEnd (-2 bohr): the absorber needs a positive length; use ZSpOptLowStart < ZSpOptLowEnd
  Input error: ZOptLowEnd: ZOptLowEnd is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them
  Input error: ZOptLowStart: ZOptLowStart is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them
  ```
- The set-up prints each low-Z absorber's own range and strength, and the points inside it. For
  case (h):
  ```text
  Low-Z CAP in Z (specular) ..... -2.000000 .. 6.000000 bohr, strength 8.98684162e-03 hartree (energy 0.300000 eV)
  ...
  Low-Z optical potential Z specular (indices) .. 0 .. 100
  ```
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/LowZD1…`, `LowZD2…`,
  `LowZD3…`, the exact bounds `LowZStartBelowZStart` and `LowZEndAboveZana`, and
  `Input.LowZAbsorbersAtTheirBoundsAreValid`) and `tests/test_report.cpp`
  (`Report.LowZD4PrintsEachAbsorbersOwnValues`).

## analysis D1: NrOfEs = 0 and 1 pass the check

**Status:** fixed upstream in r1715–r1719 (2026-11-19): the Fortran rejects NrOfEs < 2 (`'NrOfEs < 2.'`), which
also makes NrOfEs = −95 an input error. The C++ check is the same.

### The Fortran code

`InputNew.f90 ValidateInput` rejects only a negative number of analysis energies (line 574 in
5.19, 644 in 5.24; identical):

```fortran
    IF ( nenrgy /= Constants_NEN95 ) CALL ERROR( nenrgy < 0, 'NrOfEs < 0.' )
```

`ProcessInput` then computes the spacing of the analysis energies (lines 906–907 in 5.19, 978–979
in 5.24):

```fortran
    enrgi = EMinWP + epar
    denrgy = (EMaxWP - EMinWP) / FLOAT(nenrgy - 1)
```

and `Analyzor.f90 SetupAnalysis` the energies themselves (lines 112–114 in 5.19, 124–126 in 5.24;
NrOfEs = −95, `Constants_NEN95`, replaces `enrgi`, `denrgy` and `nenrgy` by the spectral grid first):

```fortran
      allocate ( ensctr(nenrgy), enscev(nenrgy), enscau(nenrgy) )
      do ien = 1 , nenrgy
        ensctr(ien) = enrgi + dfloat(ien-1) * denrgy    ! Total translational energy in eV.
```

### What goes wrong

- **NrOfEs = 1:** the spacing is (EMaxWP − EMinWP)/0 = +∞ (the build does not trap floating-point
  exceptions), and the one energy is `enrgi + 0·∞` = NaN. Every channel then counts as closed
  (`smatrix_calc` tests `skzsq > 0.0`, which is false for NaN), so the sum of the probabilities is
  0 and the reaction probability 1. The run ends normally, without an error or a warning.
- **NrOfEs = 0:** the arrays of energies are empty and every loop over them does nothing: the run
  propagates and prints no S-matrix, probability or reaction probability at all.

### Demonstration

**Case (i):** `NrOfEs : INTEGER : 1;` instead of 41, and a short propagation (`tEnd : REAL : 200.0;`,
`nStep : INTEGER : 100;` instead of 14000.0 and 7000). Both versions finish with exit status 0;
their analysis output (stdout from the line `PRINTING SCATTERING AND REACTION PROBABILITIES` on),
`averages.out` and fort.25 are identical. The input echo and the analysis print:

```text
* Number of energies that are on-the-fly analysed:             1
...
 total energy =             NaN
 trans energy =             NaN
 trans energy in z =             NaN
 total momentum =             NaN
...
 sum of the probabilities =  0.0000000
...
 the reaction probability =  0.1000000D+01
```

**Case (j):** the same with `NrOfEs : INTEGER : 0;`. Both versions finish with exit status 0. After
the header `PRINTING SCATTERING AND REACTION PROBABILITIES` the output goes straight on to the CPU
times, and no fort.25 is written.

### The C++ behaviour

- `src/input/validate.cpp` (`validate_input`) rejects NrOfEs < 2, as the Fortran since r1719 (before, the
  C++ rejected 0 and 1 and accepted −95). Both cases are input errors (exit code 2), e.g. for case (i):
  ```text
  Input error: NrOfEs: NrOfEs = 1: the analysis needs at least 2 energies, (EMaxWP-EMinWP)/(NrOfEs-1) apart; use NrOfEs >= 2.
  ```
- Until this check the C++ behaved as the Fortran: for case (i) its summary table showed
  `E_coll = -nan`, `P_scattered = 0`, `P_reaction = 1`, and for case (j) empty tables, with exit
  code 0.
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/AnalysisD1NrOfEsZero`,
  `AnalysisD1NrOfEsOne` and `NrOfEsNen95`, `Input.NrOfEsTwoIsValid`).

## analysis D4: dtAna or dtAve just below a multiple

**Status:** open in the Fortran (5.19, r1713 and r1727, the reference). The C++ rejects these inputs since
2026-11-24 (a deviation: [fortran_deviations.md](fortran_deviations.md#analysis-d4)).

### The Fortran code

`InputNew.f90 ValidateInput` compares the time intervals with the tolerant operators of `MyOperators.f90`
(absolute tolerance 1e-12; lines 566–580 in 5.19, 636–650 at r1727):

```fortran
    CALL ERROR( dtAnalysis .RLT. dtlan, 'dtAna < dtlan.' )
    ...
    Test = (NINT(FLOAT(nstep)*dtlan/dtAnalysis)*dtAnalysis) .RNE. (FLOAT(nstep)*dtlan)
    CALL ERROR( Test, 'Propagation time must be an integer times dtAna.' )
    CALL ERROR( dtAverages .RLT. dtAnalysis, 'dtAve < dtAna.' )
    Test = (NINT(dtAverages/dtAnalysis)*dtAnalysis) .RNE. dtAverages
    CALL ERROR( Test, 'dtAve must be an integer times dtAna.' )
```

`ProcessInput` then truncates the two ratios (lines 915–916 in 5.19, 990–991 at r1727):

```fortran
    nanspo = INT(dtAnalysis/dtlan)
    nave = INT(dtAverages/dtAnalysis)
```

### What goes wrong

A value up to 1e-12 below a multiple passes the tolerant checks and truncates to one less:

- **dtAna just below dtStep:** nanspo = 0. The Fortran stops with SIGFPE (the build traps no floating-point
  exception, so this is an integer division by zero) in `ndatnw = nstep / nanspo` (`AsympInfo.f90
  ReadOrWriteAsympInfo`, line 85 in 5.19, 97 at r1727), before it propagates. The C++ looped forever:
  `for (step = n; step <= nStep; step += n)` with n = 0 propagates one time step and appends a BK record on every
  pass, so the disk fills.
- **dtAve just below dtAna:** nave = 0. `mod(countPropSteps, nave)` (`Main.f90`, line 405 in 5.19, 418 at
  r1727) divides by zero after the first analysis step: SIGFPE in both codes.
- **dtAna just below k·dtStep (k ≥ 2):** nanspo = k − 1. The analysis stays consistent (both codes use nanspo ·
  dtStep as the interval of the BK coefficients and the S matrix), but when nStep is not a multiple of k − 1 the
  time loop (`DO propStep = nanspo, nstep, nanspo`) stops before tEnd. The C++ checkpoint is then labelled tEnd,
  and the Fortran's restart file stores no time (a restart takes tStart from its input), so a restart continues
  from the wrong time in both codes.
- **dtAve just below m·dtAna (m ≥ 2):** m − 1 analyses per average. Harmless: the averages carry their real times.

### Demonstration

`tests/reference/hardwall_small/DVR.inp` (tStart 0) with the lines of each case replaced; the C++ of 2026-11-24
before the check from `input.json` with the same values, 1 thread. Each run under `timeout 120` with a 1 GiB file
size limit. 5.19 and r1727 behave identically in every case (the same exit status, and the same averages where they print any).

| Case | Keys | Fortran 5.19 and r1727 | C++ before the check |
| --- | --- | --- | --- |
| (k) | dtStep 2.0, nStep 1, tEnd 2.0, dtAna = dtAve = 1.9999999999995 | exit 136 (SIGFPE) right after the analysis set-up ("Maximum diffraction order in analysis"), before the propagation | exit 124 (stopped by `timeout`): the printed time stays 0.00 while ⟨Z⟩ keeps falling; `bk_coefficients.bin` 5.6 MB after 120 s |
| (l) | dtStep 2.0, nStep 5, tEnd 10.0, dtAna 10.0, dtAve 9.9999999999995 | exit 136 (SIGFPE) after the averages at t = 0 | exit 136 (SIGFPE) |
| (m) | dtStep 0.8, nStep 3, tEnd 2.4, dtAna 2.4, dtAve 2.4 (2.4 / 0.8 = 2.9999999999999996) | exit 0; the last averages at t = 0.16000000E+01, not 2.4 | exit 0; "Steps (total / per analysis) 3 / 2", the last averages at t = 1.60, `checkpoint.bin` labelled t = 2.4 |

### The C++ behaviour

- `src/input/validate.cpp` (`validate_input`) adds three exact checks after the tolerant ones, with the run's own
  truncation (`steps_per_analysis`, `analyses_per_average`, `src/input/process.cpp`). Each case above is an input
  error (exit code 2):

  ```text
  Input error: dtAna: dtAna = 1.9999999999995 lies within the 1e-12 tolerance below dtStep = 2, so dtAna / dtStep truncates to 0 time steps per analysis; use dtAna >= dtStep.
  Input error: dtAve: dtAve = 9.9999999999995 lies within the 1e-12 tolerance below dtAna = 10, so dtAve / dtAna truncates to 0 analyses per average; use dtAve >= dtAna.
  Input error: dtAna: nStep = 3 is not a multiple of the 2 time steps per analysis that dtAna / dtStep truncates to; the run would stop before tEnd.
  ```

- Still allowed, as in the Fortran: dtAna just below k·dtStep when nStep is a multiple of k − 1 (the run reaches
  tEnd, with the analysis interval (k − 1)·dtStep throughout), and dtAve just below m·dtAna.
- The wave-packet designer has the same checks (`tools/wave_packet_designer/program/validate.py`).
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/AnalysisD4DtAnaTruncatesToZeroSteps`,
  `AnalysisD4DtAveTruncatesToZeroAnalyses`, `AnalysisD4NStepNotAMultipleOfTheTruncatedSteps`,
  `AnalysisD4DtAnaBelowTheBand`, `Input.AnalysisD4DtAnaEqualToDtStepIsAccepted`,
  `Input.AnalysisD4DtAveJustBelowAMultipleIsAccepted`); the designer's `tests/test_validate.py`.
