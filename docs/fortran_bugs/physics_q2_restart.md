# physics Q2: a restart is not checked against the run it continues

A bug of the Fortran SPO-DVR code (upstream r1727, git `3d00331`, and earlier versions as stated below), found in an
audit of the active build for errors that make a computed quantity wrong. The candidate was read by an auditor, checked
by an independent reviewer who saw the code but not the auditor's conclusion, and run where possible. Line numbers refer to
`sources.dir/SOURCES.dir/` at `3d00331`. The short entry is in [fortran_deviations.md](../fortran_deviations.md) under
the same ID; the other bugs of the audit are listed in [README.md](README.md).

**Status:** open in the Fortran (5.19 and r1727); not yet reported upstream (2026-11-25). The C++ refuses each of
these restarts.

## The demonstration runs

- **Executable:** the reference install `r1727-3d00331` (upstream git `3d00331` built with `O3FLAGS = -O3`;
  [tools/fortran_reference/README.md](../../tools/fortran_reference/README.md#the-reference-install)).
- **Inputs:** the `DVR.inp` of a reference case in `tests/reference/` with the keys of each case changed or appended
  (and its data files, as `tools/fortran_reference/reference_config.sh` `stage_inputs` stages them).
- **How:** each run in a new scratch directory through `reference_config.sh` (`run_clean`: a clean environment, 4
  OpenMP threads). The C++ runs used the program at git `48e5bdd`, 4 threads.

## The Fortran code

A restart (`NewCalc : FALSE`, nread = 1) reads the wave function from fort.8 and the BK coefficients of the earlier
run from fort.22. fort.8 holds no time (`WaveFunction.f90 SaveWaveFunc`, lines 1258–1264); fort.22 holds the number of records and the
analysis interval `deltat` of the earlier run. `Analyzor.f90 GetSMatrixElemsOfOldRun` adds the earlier records to the
S-matrix time integral (lines 447–470):

```fortran
      read( ibkin_file ) ndatold     ! rec1
      read( ibkin_file ) deltatold   ! rec2
      ...
      dtmul = nanspo * dtlan
      do istep = 1 , ndatold
        read (ibkin_file) balcf
        time = istep * deltatold
        do ien = 1 , nenrgy
          phaz = exp( ima * enscau(ien) * time )
          ...
              smtrx(is,iv,ien) = smtrx(is,iv,ien) + balcf(is,iv) * phaz * dtmul
```

The new records are labelled from the input's tStart (`Main.f90` line 402):

```fortran
           the_time = tottym + dfloat(propStep - nstep) * dtlan
```

ValidateInput checks only that a restart has tStart ≠ 0 and that tEnd = tStart + nStep·dtStep (`InputNew.f90` lines
621–626). The basis labels of fort.22 (records 4–9) are read into dummies (`RotDifBasis.f90 ReadFromFile`), and each
old record is read into `balcf` of the size the new input gives. `AsympInfo.f90 ReadOrWriteAsympInfo` copies the old
records into the new BK file with the new `deltat` in its record 2 (lines 137–164).

## What goes wrong

- **dtAna changed:** the old records keep their times (`istep*deltatold`) but are weighted with the new interval
  (`dtmul`), so the earlier part of the time integral is scaled by dtAna_new / dtAna_old. The new BK file then labels
  every record with the new interval, so a later restart or fine analysis puts the old records at the wrong times.
- **tStart not the end of the earlier run:** nothing compares tStart with ndatold·deltatold. The new records get a
  phase e^{iE(tStart − T_old)} relative to the old ones, and the two parts of the S matrix interfere wrongly.
- **Analysis basis changed** (a different EMaxWP, JAna, VAna or DifOrdAna changes the number of analysed states):
  with fewer states the read takes the leading part of each old record, so coefficients land in the wrong states (in
  column-major order); with more it fails ("I/O past end of record"). From the code and gfortran's record semantics;
  not run.
- No message in any case.

## Demonstration

`hardwall_small` run to t = 7000 (tEnd 7000.0, nStep 3500), then three restarts from its fort.9 and fort.21 (copied
to fort.8 and fort.22) to the end: (ok) `NewCalc` FALSE, tStart 7000.0, nStep 3500, the other keys unchanged; (a) the
same with dtAna 20.0 instead of 10.0; (b) tStart 7100.0 and tEnd 14100.0. Every run exits with status 0 and prints no
error or warning beyond the two of every run of this case ("Rotational basis too small…", "Diffractional basis too
small…"). The hard wall reflects everything, so the sum of the probabilities should be 1 at every energy; the unsplit
run (`tests/reference/hardwall_small/fortran_stdout.txt`) gives 0.994–1.006, and (ok) reproduces it exactly
(and `fortran_restart_stdout.txt`, which the reference tooling makes the same way):

| E_tot (eV) | unsplit | restart (ok) | (a) dtAna 20 | (b) tStart 7100 |
| --- | --- | --- | --- | --- |
| 0.45838 | 0.9939394 | 0.9939394 | 1.5724319 | 0.3356800 |
| 0.55838 | 1.0023228 | 1.0023228 | 1.5444369 | 1.1635500 |
| 0.65838 | 1.0011205 | 1.0011205 | 3.2709078 | 0.3155920 |
| 0.75838 | 0.9994264 | 0.9994264 | 4.5814986 | 1.8221434 |
| 0.85838 | 0.9972464 | 0.9972464 | 5.2341844 | 2.5102554 |

Largest difference from the unsplit run over the 41 energies: 0 (ok), 4.24 (a), 1.51 (b).

## The C++ behaviour

- A restart from Fortran-format files checks fort.22 against the input: record 1 (the numbers of states), record 2
  (dtAna, exactly), the basis labels and the size of every record (`src/io/fortran_files.cpp` lines 152–185), and
  tStart against the time of the last BK record (`src/app/run.cpp` lines 173–175). Native checkpoints carry the time
  steps and the analysis basis among their restart parameters (`src/io/checkpoint.cpp`).
- Cases (a) and (b) with the C++ (`RestartFormat` "fortran", `RestartDir` the Fortran run's directory) stop with exit
  code 4:

  ```text
  Checkpoint error: the analysis time step (dtAna) differs from the one in '…/fort.21'
  Checkpoint error: tStart = 7.1000000000000000e+03 but the checkpoint in '…' is at t = 7.0000000000000000e+03
  ```

## Proposed fix

In the restart branch, compare fort.22 with the input before using it: `deltatold` with `nanspo*dtlan`, `ndatold*deltatold`
with tStart, records 1 and 4–9 with the current basis; stop with a message on a mismatch.
