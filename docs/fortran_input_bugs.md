# Fortran input checks: bugs and how the C++ handles them

The C++ program reproduces the input checks of the Fortran code (`InputNew.f90` ReadInput and
ValidateInput, and the grid checks of `Setup.f90`; see [input_schema.md](input_schema.md),
"Validation"). A few of those checks are wrong. This document describes each one:
- the Fortran code;
- what goes wrong;
- a demonstration with the unmodified Fortran executables;
- what the C++ does instead.

The short entries are in [fortran_deviations.md](fortran_deviations.md), under the same IDs.
Later sections will cover more input checks: the optical potentials at the low end of the Z grid
(added in upstream 5.24), and the default absorber strength that the Fortran applies silently.

## The demonstration runs
- **Executables:**
  - 5.19, the version the port started from: `~/work/SPO-DVR_v5.19/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x`
    (gfortran 15, `-O3`, OpenMP, FFTW);
  - 5.24: upstream master at commit dce0ecc (svn r1713), built with `O3FLAGS = -O3` in `GNU.mk`
    (the default flags enable fast-math); now the reference of the port
    ([validation.md](validation.md#reference-data)).
- **Input:** `tests/reference/hardwall_small/DVR.inp` with the flux keys of each case appended.
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
- Case (a) passes the input checks. The flux analysis itself is not ported yet, so a propagating
  run then stops with exit code 3:
  ```
  Not supported: doFluxAnalysis = true: the flux analysis is not implemented yet
  ```
  Analysis-only and potential-only runs accept the flux keys and ignore them.
- Tests: `tests/test_input.cpp` (`Checks/InputViolation.IsRejected/FluxEMax`,
  `Input.FluxEnergyRangeInsideTheWavePacketIsValid`) and `tests/test_cli.cpp`
  (`FluxAnalysisGuard.*`).

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
