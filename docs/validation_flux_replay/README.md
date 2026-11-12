# Flux analysis: replay of one Fortran analysis step

One-off checks of the C++ flux analysis (`src/analysis/flux_analysis.*`; Fortran `FluxAnalyzor.f90`)
against the reference Fortran (5.24 r1713 built with `-O3`, [validation.md](../validation.md#reference-data)),
made on 2026-11-12 before the flux analysis was wired into the program. They test the arithmetic of
the flux analysis to the last bit, on the Fortran's own wave function, so that a later comparison of
whole runs only has to account for the difference in ψ. Nothing here runs in ctest and no data is
committed: [run_checks.sh](run_checks.sh) rebuilds and reruns all three checks in a scratch directory.

## Results (2026-11-12)

| Check | What is compared | Result |
|---|---|---|
| Replay | hardwall_small, one analysis step, from the Fortran's ψ: c_n·ψ_c and c_n·∂ψ_c/∂r for 3 energies × 113 planes × 240 (rotational function, x, y) points, against the Fortran's `flux_coeffs.dat` and `flux_coeffs_deriv.dat` | 162 720 of 162 720 complex values equal, all bit for bit (32 544 of them nonzero); the C++ writer's two files are byte-identical to the Fortran's |
| 1 | b(k) at 401 momenta from below to above the wave packet's range; the 56 weights D(ir) of the r-derivative at ir_c = 9 (nR 56, dR 0.15 bohr); the number of flux energies of 4 ranges | 401 of 401 and 56 of 56 bit for bit; the Fortran prints 4, 6, 3 and 22 energies, the counts the unit test `FluxEnergies.CountDropsTheEndPointAsTheFortran` requires of `flux_energy_count` |
| 2 | five analysis steps accumulated on a synthetic ψ (2 rotational functions, 1 × 2 (x, y) points, nR 8, 2 planes, 3 energies): φ and φ′; P(E) and P(x, y; E) | φ and φ′ 48 of 48 bit for bit; P(E) 3 of 3 and P(x, y; E) 6 of 6 within half a printed unit of the Fortran's E15.8 files |

Check 2 also showed, once, with a modified C++ build (not part of run_checks.sh) that the complex
division of c_n matters: with A_n = k·conj(d)/|d|² (d = 2πM·b(k)) instead of the division, 3 of the
48 values are bitwise equal.

## The replay

1. **The Fortran part** ([run_fortran.sh](run_fortran.sh)). The reference executable runs
   hardwall_small with the flux keys (rFluxAnalysis 2.0 bohr, grid point 9 of the r grid; EMinFluxAna
   0.2 = EMinWP, EMaxFluxAna 0.6 = EMaxWP, DeltaEFluxAna 0.15 eV: the energies 0.2, 0.35 and 0.5 eV)
   in two parts, split by a restart:
   - part1: t = 0 … 7000 (3500 steps of 2);
   - part2: t = 7000 … 7010, one analysis step (dtAna 10), restarted from part1's fort.9 and fort.21,
     with part1's flux coefficient files zero-filled ([zero_payload.py](zero_payload.py) keeps the
     gfortran record markers, so the layout is exact). After the step the files hold 0 + c_n·ψ_c and
     0 + c_n·∂ψ_c/∂r;
   - part2_noflux, the control: part2 without the flux keys. Its fort.21 and fort.9 are
     byte-identical to part2's: the flux analysis does not change the propagation.
2. **The wave function and the time.** In each analysis step Main.f90 calls `Propagate`,
   `BalintKurtiCoeffs` and `storeFluxCoeffs(the_time)`, and it saves ψ (`SaveWaveFunc`) after the time
   loop: part2's output fort.9 is the ψ that `storeFluxCoeffs` used. `the_time = tottym +
   dfloat(propStep − nstep)·dtlan` is tEnd = 7010 for the one step, the time the C++ time loop gives
   that step.
3. **The C++ part** ([replay_flux.cpp](replay_flux.cpp)): the part-2 input (hardwall_small's
   `input.json` with the same keys), ψ from part2's fort.9 (`read_fortran_wave_function`), and eviroi
   from record 2 of its fort.21. The Fortran's eviroi, not the C++ gas-phase energy, so that the
   replay isolates the flux arithmetic. Then `make_flux_energies` (edifin = 0 here),
   `flux_derivative_weights`, ψ_c and `flux_r_derivative` per (rotational function, x, y, Z), one step
   of `FluxAccumulator`, and the Fortran's files read with `read_fortran_flux_coefficients`.
4. **The comparison** is by value (real and imaginary parts with ==), as for every real × complex
   product of the port: gfortran can flip the sign of an exact zero there
   ([fortran_deviations.md](../fortran_deviations.md), "signed zeros"). Here the bits agree too,
   because adding the zero-filled start value turns −0 into +0 on both sides.

**The zeros.** 80 % of the values are exact zeros: the hard wall couples the initial state (j = 0,
mJ = 0) only to the rotational functions with mJ = 0 (j = 0, 2 and 4: 3 of the 15), and ψ stays exactly
0 in the others. The nonzero values, 16 272 per file, are the 3 energies × 113 planes × 3 rotational
functions × 16 (x, y) points. Check 2 covers a ψ without zeros.

**Control.** With the time of the step's start (7000) instead of its end, only the 65 088 zeros of
φ agree, and the largest difference is 0.23 of max |φ|: the replay sees the time of the step.

**Not covered here:** the C++ propagation (ψ is the Fortran's) and the C++ backends' ψ_c and ∂ψ_c,
which their own tests compare with `flux_r_derivative`.

## Checks 1 and 2

Small Fortran drivers call the reference routines directly and print their results as hexadecimal bit
patterns; a C++ program reads them and compares.
- **Check 1** ([flux_f1_check.f90](flux_f1_check.f90), [flux_f1_check.cpp](flux_f1_check.cpp)):
  `getBkCoeffWavePacket` for the hardwall_small packet (EMinWP 0.2, EMaxWP 0.6 eV, H₂, Z0 11.44
  bohr) against `wave_packet_momentum_amplitude`; `getWaveFuncDerivAtZr` applied to unit vectors
  against `flux_derivative_weights`; `setupFluxAnalyzor` prints its energy counts (0.2–0.6 eV in steps
  of 0.1, 0.08 and 0.15, and 0.40–0.95 in steps of 0.025), to be read against the unit test.
- **Check 2** ([flux_f2_check.f90](flux_f2_check.f90), [flux_f2_check.cpp](flux_f2_check.cpp)):
  `storeFluxCoeffs` on a synthetic ψ with nonzero edifin and eviroi, then `readWriteFluxCoeffs` and
  `computeFluxProbability`; the C++ `FluxAccumulator` (2 slots, so the flush is exercised) and
  `flux_results` on the same ψ and times.

## Rerun

```
docs/validation_flux_replay/run_checks.sh build /tmp/flux_checks   # a new directory; about 3 minutes
```

It needs the reference install ([tools/fortran_reference/README.md](../../tools/fortran_reference/README.md)),
gfortran and FFTW to link the Fortran drivers, and a built `build` directory (the default, strict
build). It prints each check's result and exits with 0 when all pass. The pieces on their own:

| File | Use |
|---|---|
| [run_checks.sh](run_checks.sh) `BUILD_DIR SCRATCH_DIR` | everything below |
| [run_fortran.sh](run_fortran.sh) `SCRATCH_DIR` | the Fortran part of the replay (part1, part2, part2_noflux) |
| [zero_payload.py](zero_payload.py) `IN OUT` | a gfortran file with its record payloads zeroed |
| [build_helper.py](build_helper.py) `BUILD_DIR HELPER.cpp OUT_DIR` | compiles and links one C++ program against a build |
| [replay_flux.cpp](replay_flux.cpp) `input.json RUN_DIR [OUT_DIR]` | the C++ part of the replay |
| flux_f1_check.*, flux_f2_check.* | checks 1 and 2 |
