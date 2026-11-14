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
Fortran keeps for an absorber energy of 0 (absorber D1), and the optical potentials at the low end
of the Z grid that upstream 5.24 added (low-Z D1–D4).

## The demonstration runs
- **Executables:**
  - 5.19, the version the port started from: `~/work/SPO-DVR_v5.19/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x`
    (gfortran 15, `-O3`, OpenMP, FFTW);
  - 5.24: upstream master at commit dce0ecc (svn r1713), built with `O3FLAGS = -O3` in `GNU.mk`
    (the default flags enable fast-math); now the reference of the port
    ([validation.md](validation.md#reference-data)).
- **Input:** `tests/reference/hardwall_small/DVR.inp` with the keys of each case appended (or, where
  a case says so, one line of it replaced).
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
```
doFluxAnalysis : LOGICAL : TRUE;
rFluxAnalysis : REAL : 3.05;
EMinFluxAna : REAL : 0.3;
EMaxFluxAna : REAL : 0.5;
DeltaEFluxAna : REAL : 0.01;
```
Both versions stop right after the echo of the input, with exit status 0:
```
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
  ```
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
  ```
  Input error: EMaxFluxAna: EMaxFluxAna > EMaxWP (0.8 eV > 0.6 eV): the flux analysis energies must lie inside the wave packet's energy range; use EMaxFluxAna <= EMaxWP
  ```
- Case (a) passes the input checks, and a new run or a restart computes the flux analysis.
  Analysis-only and potential-only runs accept the flux keys and ignore them.
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/FluxEMax`,
  `Input.FluxEnergyRangeInsideTheWavePacketIsValid`) and `tests/test_cli.cpp`
  (`FluxAnalysisRuns.*`).

## flux Q2: a zero DeltaEFluxAna passes the check

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
```
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
  ```
  *****
   0.15707963E+01  C2v
      4 0.15599150E+01
      4 0.13099150E+01
  ```

### The C++ behaviour
- `validate_input` stops when DeltaEFluxAna ≤ 0. The test is tolerant (`.RLE.`), like the other
  "≤ 0" checks of grid spacings (dZ, dR).
- The error names the key (exit code 2). For case (c) with a range inside the wave packet
  (EMaxFluxAna = 0.5, so that the flux Q1 check passes first):
  ```
  Input error: DeltaEFluxAna: DeltaEFluxAna <= 0 (0 eV): the flux analysis energy grid needs a positive spacing; use DeltaEFluxAna > 0
  ```
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/FluxDeltaZero` and
  `FluxDeltaNegative`).

## absorber D1: an absorber energy of 0 keeps the default strength

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
```
* Optical potential in Z absorbs with coefficient:          0.43800882E-02
```

### The C++ behaviour
- An energy of 0 means no absorption: the strength is 0 and every absorption factor is exactly 1
  (`src/input/process.cpp`, `process_input`). For case (d) the set-up prints
  ```
  CAP strengths (Z, r, Zsp) ..... 0.00000000e+00, 2.38963250e-02, 1.82989643e-02 hartree
  ```
- Tests: `tests/test_propagation_tables.cpp` (`PropagationTables.AbsorberD1EnergyZeroGivesNoAbsorption`).

## The low-Z absorbers of upstream 5.24

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
- **low-Z D2:** with Start = End the absorber is switched on with length 0. The strength is not
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
```
* Optical potential in Z at lower end of specular grid absorbs with coefficient 0.10000000E-01
```
and `averages.out` becomes NaN after the first analysis interval:
```
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
```
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
  ```
  Input error: ZOptLowStart: ZOptLowStart > ZOptLowEnd (4 > 3 bohr): the absorber would be empty; use ZOptLowStart < ZOptLowEnd
  Input error: ZSpOptLowStart: ZSpOptLowStart = ZSpOptLowEnd (-2 bohr): the absorber needs a positive length; use ZSpOptLowStart < ZSpOptLowEnd
  Input error: ZOptLowEnd: ZOptLowEnd is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them
  Input error: ZOptLowStart: ZOptLowStart is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them
  ```
- The set-up prints each low-Z absorber's own range and strength, and the points inside it. For
  case (h):
  ```
  Low-Z CAP in Z (specular) ..... -2.000000 .. 6.000000 bohr, strength 8.98684162e-03 hartree (energy 0.300000 eV)
  ...
  Low-Z optical potential Z specular (indices) .. 0 .. 100
  ```
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/LowZD1…`, `LowZD2…`,
  `LowZD3…`, the exact bounds `LowZStartBelowZStart` and `LowZEndAboveZana`, and
  `Input.LowZAbsorbersAtTheirBoundsAreValid`) and `tests/test_report.cpp`
  (`Report.LowZD4PrintsEachAbsorbersOwnValues`).
