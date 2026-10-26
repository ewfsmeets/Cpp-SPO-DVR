# Cpp-SPO-DVR

Modern C++ rewrite of the Fortran 6D quantum wavepacket code in `~/work/SPO-DVR`
(diatomic molecule–surface scattering: Z, r, x, y, θ, φ; split-operator propagation;
FFT grids in Z/r/x/y, Gauss-Legendre/FFT DVR in θ/φ; reduced "specular" (Z,r) grid
for the incoming channel).

## Goals (in priority order)
1. Reproduce Fortran results to within numerical limits (target: ~1e-10 relative on
   wavefunction norms/averages and S-matrix elements for the hard wall test system).
2. Clean, structured, readable, easy to modify, with an extensive unit test suite.

## Fortran reference
- Sources: `~/work/SPO-DVR/procfiles.dir/SOURCES.dir/` (active build: GFORTRAN.mk, OpenMP, FFTW,
  `-fdefault-real-8` ⇒ every Fortran `REAL`/`COMPLEX` is double).
- Active input reader: `InputNew.f90` (`__USE_NEW_INPUT_FORMAT__` defined). `Input.f90` is dead.
- Ignore: MPI code (`MPI_*.f90`, `#if defined(USE_MPI)` blocks), `*_hack.f90`,
  all PES modules except `HardWallPES.f90`, CRPlibrary, PES4Cu100.
- **`fftw_2_scs.f90` is NOT dead**: it is the FFTW wrapper behind every `ffts.f90` call in the active
  build (FFTW_ESTIMATE plans). It emulates the SCS sign convention: r2c with Sign = +1 returns
  `Scale·conj(FFTW r2c)`; c2r with Sign = −1 conjugates the input before FFTW c2r; output ×Scale.
  `Inverse…FFT` uses scale `1/(Scale·N)`. Check the complex (`zzfftm`) wrapper before porting the
  propagator FFTs.
- Key files ↔ C++ responsibility:
  - `Main.f90` → program driver / run modes
  - `InputNew.f90` (ReadInput/ValidateInput/ProcessInput) → `input/` (JSON + validation)
  - `Setup.f90` (MiscellaneousSetup) → grid + derived-parameter construction
  - `setupdvr.f90`, `angulargrid.f90`, `legendre.f90`, `transform.f90` → angular grid + FBR⇄DVR
  - `SplitOperator.f90` → propagator; `KineticOp.f90`, `GridEnergies.f90` → kinetic operators
  - `Potential.f90` (setup + PerformEXPPotentialOp) → potential operator; `HardWallPES.f90` → hard wall
  - `FourierGridHamiltonian.f90` + `jacobi.f90` → `gasphase/` FGH (C++ diagonalises with LAPACK)
  - `TwoGrids.f90` → specular⇄scattering projection; `Absorption.f90` → optical potentials
  - `MomentumShifting.f90`, `WaveFunction.f90` (CoherentStateFunc, save/read), `GasPhase.f90` (FGH)
  - `Analyzor.f90`, `AsympInfo.f90`, `RotDifBasis.f90`, `FluxAnalyzor.f90`, `Averages.f90` → analysis
- Reference executable (built with gfortran, OpenMP, FFTW, -O3):
  `~/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x`. It reads `DVR.inp` from the
  current working directory and writes output files there — always run it in a scratch directory,
  never inside the Fortran source tree. Rebuild with `cd ~/work/SPO-DVR/procfiles.dir && ./DoMake`.
- Fortran `DVR.inp` syntax: `Name : value;` (`#` comment lines, strings in `"…"`, logicals TRUE/FALSE).
  Only fields pre-registered with `Set…Field` in `InputNew.f90 ReadInput` may omit the type; all
  others need `Name : INTEGER|REAL|LOGICAL|STRING : value;` (else "ERROR CODE RETURNED: 4").
  Valid FFT sizes (`ffts.f90 GetNearestFFTSize`): even, 2^h·3^i·5^j·7^k·11^m with i ≤ 2, j,k,m ≤ 1.
- When porting a routine, cite the Fortran origin in a short comment (`// Fortran: SplitOperator.f90 Propagate`).

## Algorithm essentials (do not change semantics without a reference comparison)
- Propagate one analysis step = `[T(dt/2) V(dt) T(dt/2)]^nanspo` with adjacent half steps merged:
  `T/2 V (T V)^(nanspo-1) T/2`; on the final T/2 no potential is applied.
- One T(dt') block: (1) specular grid ops (if nproj>0) → (2) FFT kinetic in Z then r on full scattering
  WF → (3) per (iz,ir): momentum shift, FFT x,y, ×exp(-iT_xy dt'), iFFT, unshift; then potential op.
- Potential op per (iz,ir): ×exp(-iV_optZ dt)·exp(-iV_optR dt); if iz≤nzpot && nproj>0 add specular
  part to channel `ichn`; iz≤nzpot: rot(dt/2)·[FBR→DVR · exp(-iV dt) · DVR→FBR]·rot(dt/2); else
  rot(dt/2)·exp(-i vasy(r) dt)·rot(dt/2); then project specular part back out (TwoGrids Add/Sub).
- Stored potential is premultiplied by −dt (`vdvr = −V·dt`), only for Z < zasymp (`nzpot` points).
- FFT naming in Fortran is inverted: "Backward" = coordinate→momentum. Scale factors are explicit:
  z,r: `d/√(2π)`; x,y: `√(dx/nx)`. Verify sign convention against `ffts.f90` before porting.
- Run modes (Fortran `nread`): 0 new, 1 restart, 2 analysis-only, 6 potential-only.
  `nproj`: 1 = merge specular into scattering grid once norm beyond scat grid < tolsp; 2 = two grids
  throughout. `nproj = 0` is rejected by the input reader.

## Architecture rules
- **Potentials**: abstract `Potential` interface (V(x,y,Z,r,θ,φ) in eV, asymptotic V(r), system
  parameters: alat, blat, zasymp, skew angle, point group). Concrete impls only. Currently only
  `HardWallPotential`. Factory selects by `SystemID` string.
- **FFT**: abstract FFT interface (plan creation for strided/batched 1D complex transforms, in-place
  execute, explicit scale). Concrete impl only `FftwFft` (FFTW3, double). No other library code.
- **Grids**: grid objects (uniform 1D grid, momentum grid, specular/scattering Z grids, angular grid,
  rotational FBR basis) expose points, spacing, size, k-values, kinetic energies. No free-floating
  grid arrays.
- **Arrays**: multi-dimensional data = flat `std::vector` + explicit index/offset helpers. Keep the
  Fortran memory order (first index fastest): ScatWave `(irot, iy+ny*ix, ir+nr*iz)`, SpecWave `(ir, iz)`,
  potential `(iφ+nφ*iθ, ixy, ir+nr*iz)`. This keeps FFT strides and summation order identical.
- **Constants**: one simple `Constants` object (header-only, `constexpr`): π, 2π, √π, √(2π), deg→rad,
  imaginary unit, `nen95` sentinel (−95). Values copied verbatim from `constants.f90`.
  `Constants_EPS` is computed at run time in Fortran (2⁻⁵⁴): use `fortran_machine_epsilon()`.
- **UnitConversion**: one simple object holding the conversion factors plus straightforward named
  procedures (`ev_to_hartree(x)`, `hartree_to_ev(x)`, `amu_to_me(x)`, `bohr_to_angstrom(x)`,
  `angstrom_to_bohr(x)`, `hartree_to_wavenumber(x)`, `deg_to_rad(x)`, …). Two named instances:
  - `legacy` — the factors the Fortran actually uses almost everywhere: `chaev = 27.2114` (Ha→eV),
    `camau = 1822.89` (amu→mₑ), `cboan = 0.529177` (bohr→Å), `chacm = 219475.0` (Ha→cm⁻¹).
  - `precise` — `Constants_eV2Ha = 0.036749308286`, `u2aum = 1822.88853006`,
    `Bohr2Ang = 0.529177249`, `Ha2cmInv = 219475.0`.
  The Fortran mixes them: `FluxAnalyzor.f90` uses `precise` eV↔Ha, everything else uses `legacy`.
  Reproduce this exactly (pass the right instance explicitly); never hard-code a factor elsewhere.
- **Real comparisons**: Fortran `.RGT./.REQ./…` (`MyOperators.f90`) use an *absolute* tolerance of
  1e-12. Provide equivalent helpers and use them in input validation so accept/reject decisions match.
- **Pointers**: the only pointer types allowed are `std::unique_ptr` and `std::shared_ptr`. No raw
  pointers (owning or non-owning), no `new`/`delete`. Non-owning access uses references or
  `std::span`. The only exceptions are the C-API boundaries inside the FFTW implementation
  (`src/fft/fftw_fft.cpp`: `fftw_plan`, `vector.data()` passed to FFTW calls; plans wrapped in
  `std::unique_ptr` with a custom deleter `fftw_destroy_plan`) and the LAPACK backend
  (`src/linalg/lapack_eigensolver.cpp`). Keep them local to those files.
- **Eigensolvers**: abstract `SymmetricEigensolver` interface; only concrete impl `LapackEigensolver`
  (dsyevd). Other backends (e.g. a port of `jacobi.f90`) plug in without touching callers.
- **No globals.** Configuration is an immutable struct produced by input processing; pass explicitly.
- Parallelism: OpenMP only (no MPI). Results must not depend on thread count.
- Precision: `double` / `std::complex<double>` everywhere.

## Input / CLI
- JSON input file; keys mirror the Fortran `InputNew.f90` field names where sensible (document any
  renames in the input schema doc/test). Unknown keys are an error.
- Reproduce ALL consistency checks from: `InputNew.f90` ReadInput + ValidateInput, `Setup.f90`
  (ZAna/rFlux on grid, inside absorber), `setupdvr.f90` (J parity), `WaveFunction.f90`
  CoherentStateFunc (initial norm outside region < tolsp), `transform.f90` (nφ ≥ 2·MJMax+1 — silent
  RETURN in Fortran, must be a hard error here). Each check gets a unit test (valid + violating case).
- Invalid input → clear error message naming the parameter, non-zero exit. Never silently continue.
- CLI: `spodvr <input.json> [--threads N] [--output-dir DIR] [--no-color] [--help] [--version]`. Keep it minimal.
- Program behavior: read input, run to completion (like the Fortran code).

## Output
- Output does NOT need to match the Fortran text verbatim. Make it modern and readable: clear section
  headers, aligned tables, units in column headers, sensible number formatting, optional ANSI color
  (auto-disabled when stdout is not a TTY or with `--no-color`).
- Numerical content must still be comparable to the Fortran output: every quantity used as a
  regression gate must be printed (or written to a data file) with full precision (`%.16e`).
- Keep presentation separate from computation (a small reporting/formatting layer); numerical
  modules return data, they do not print.

## Acceptance gates (per milestone)
Reference case: `tests/reference/hardwall_small/DVR.inp` (Fortran) — Fortran stdout stored next to it.
**All gate values below are specific to the parameters in that `DVR.inp`** (grids, masses, JMax, …).
If `DVR.inp` changes, the Fortran reference must be regenerated and the gate values updated; the
C++ regression tests must use the equivalent JSON input (`tests/reference/hardwall_small/input.json`)
derived from exactly the same parameters.
1. **Hard wall PES — gas-phase potential**: V(r) at Z = zasymp = 7.0, θ = π/2 (Fortran
   "INFO FROM SETTING UP THE ROVIBRATIONAL ASYMPTOTIC POTENTIAL" table, r = 0.8 … 9.05, 56 points)
   reproduced to ~1e-14 relative (full 6D potential ~1e-13, see the quad-precision quirk). Sanity: V(r = 1.4) = 0 exactly.
2. **Gas-phase rovibrational states** (`GasPhase.f90` / `FourierGridHamiltonian.f90`): energies
   e(v,j) for v = 0,1 and j = 0…4 (Fortran "INFO FROM SETTING UP THE ASYMPTOTIC ROVIBRATIONAL
   STATES"), e.g. e(0,0) = 0.258378 eV, e(1,0) = 0.762842 eV; also ⟨r⟩, psi(1), psi(nr), norms.
   The printed Fortran values have only 6 digits — dump full precision from Fortran when a tighter
   gate is needed.
3. Later gates (to be defined): initial wavepacket norms, averages vs time, S-matrix / probabilities.
   Open item: the associated Gauss-Legendre quadrature for |mJ| > 0 (`gauss_associated_legendre`)
   is only checked for exactness and against the earlier C++ port, not against Fortran. The
   rotational FBR basis needs it, so that milestone needs a Fortran reference for it.

## Known Fortran bugs/quirks (fixed in C++; record every deviation here)
- `nZSpecular` read depends on `nproj` before `ProjOp` is read (always uses default 2).
- Error message `'EMaxFluxAna < EMinWP'` for the check `EMaxFluxAna < EMaxWP`.
- `InitGenComplexTransform` silently returns on invalid nφ.
- `tolsp` (`TolSp`) has no default in `globaldata.f90` (uninitialised if absent). C++: required key.
- `nZSpecular < nZ` is not checked by Fortran (out-of-bounds zgrid access). C++: hard error.
- Kept on purpose (Fortran semantics reproduced): `MJMax <= 0` ⇒ the basis/grid uses MJMax = JMax
  (`setupdvr.f90`/`transform.f90` "MjMax < 1"), while ValidateInput checks |MJIn| ≤ |MJMax| with the
  raw value. `VAna < 1` and absent `JAna` are replaced by the program-determined open-channel limits
  (`GasPhase.f90`).
- Not ported: `.pes` file caching (`ReadPesFromFile`/`WritePesToFile`); the potential is always
  computed. `UseRCutOff` is silently forced off by Fortran `ProcessInput`; C++ rejects
  `UseRCutOff: true` / `RCutOff`. Also rejected (NotSupportedError): `SystemMix` (+`SystemID2`,
  `SystemMixCoeff`), `UseSymmetry`/`nPhisOptimal` (C2v symmetry is disabled in Fortran anyway),
  `write2DPES` (+ coordinates), `PESMixingAlfa`, `CRPlibDataSet`.
- Giving both `ZOptEnergy` and `ZOptFactor` (same for R/ZSp) is an error (Fortran silently ignores
  the factor).
- The check `ZSpOptStart < ZAsymptotic` is unreachable (implied by the Z0 checks); kept for parity.
- `-fdefault-real-8` also promotes `DOUBLE PRECISION`, `DFLOAT`, `DBLE` and `…d0` literals to
  REAL(16) (quad precision). Expressions such as `rmin + dfloat(i-1)*delr` (Setup.f90 grids) are
  therefore evaluated in quad and rounded once to double, i.e. correctly rounded (verified: the
  Fortran r grid differs from plain-double `rmin + i*dr` by 1 ulp at 13 of 56 points). Accepted by
  the user as within numerical limits: C++ keeps plain double grids. Effect on gate 1: vasy table
  agrees to ≤1e-14, the full 6D potential to ~5e-14 relative (gate tolerance 1e-13).
  The same applies to every momentum quantity in Setup.f90 (lines ~129–131, 454–525): delkz,
  delksp, delkr, the kinetic energies ekinz/ekinr/ensp used in exp(−iT·dt), and zk/rkg/zksp. These
  are ~1 ulp off in C++ too (≈1e-12 phase after 7000 steps): first suspect if a propagation gate
  drifts. When porting any routine, look for `DFLOAT`, `DBLE`, `…d0` literals and
  `DOUBLE PRECISION`; all are quad precision in this build.
- The CoherentStateFunc check (initial norm outside region < tolsp) is not done yet; it belongs to
  the wave-packet milestone.

## Build & test
```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```
- C++20, CMake ≥ 3.24, Ninja. GCC 13 on this machine.
- Dependencies: FFTW3 (system, `libfftw3-dev`), LAPACK (system, `liblapack-dev`), nlohmann/json and
  GoogleTest via FetchContent.
- Compiled with `-ffp-contract=off` (no FMA contraction) to stay comparable with the Fortran build.
- Fortran reference data: `tools/fortran_reference/run_setup_reference.sh` (potential-only run with
  DoDump; writes `fortran_setup_stdout.txt`, `debug_asyp.dat`, `pes_reference.txt` into
  `tests/reference/hardwall_small/`). Regression tests skip when these files are missing.
- Warnings on (`-Wall -Wextra -Wpedantic`); treat as errors in CI.
- `tools/check_bitwise_refactor.sh [base-ref]`: bitwise before/after check for refactors (see
  Conventions).
- Unit tests live in `tests/`, one test file per module; regression tests compare against reference
  output from the Fortran code (stored under `tests/reference/`, small hard wall cases only).

## Readability (as important as correctness)
- Write code for a human reader first: a physicist who knows the Fortran code should be able to
  follow a C++ function without reading its Fortran origin side by side.
- Comments: concise but descriptive. Say *what* a block computes and *why* (physics meaning,
  Fortran quirk being reproduced), not how the syntax works. One short line is usually enough;
  avoid long narrative comments and avoid restating the code.
- Short functions with one job; name intermediate quantities instead of nesting long expressions
  (keep the Fortran evaluation order where it matters for bitwise results, and say so).
- Prefer plain loops and simple data structures over clever templates, variant tables,
  member-pointer tricks or dense lambdas. Avoid deep nesting; return early.
- Keep index arithmetic in named helpers (`Index3D`, grid accessors), not inline in kernels.
- No cast noise: `-Wno-sign-conversion` is set, so index containers with `int` directly and let
  int→double convert implicitly. But declare floating-point variables as `double`, never `auto`,
  when initialised from an `int` (`auto ord = order;` silently turns `amj / ord` into integer
  division).
- Each `.cpp` file starts with a 1–3 line comment stating its responsibility and Fortran origin;
  headers document each type/function they declare (with its Fortran origin).

## Conventions
- Atomic units internally; convert eV/Å/amu only at input/output boundaries, always through
  `UnitConversion` (see Architecture rules).
- Names: types `PascalCase`, functions/variables `snake_case`, members `trailing_underscore_`.
- Prefer descriptive names over Fortran abbreviations (`dtlan` → `time_step`), but note the Fortran
  name in the declaration comment when porting.
- Do not "optimize" summation order or fuse loops in numerical kernels before a regression test
  exists for them.
- Refactors that are meant to keep the numbers unchanged must leave them *bitwise* unchanged.
  Passing tests is not enough (e.g. a wrong Newton starting guess still converges). Run
  `tools/check_bitwise_refactor.sh [base-ref]`, which compares the potential-only output files
  of the base commit and the working tree (1 and 4 threads). Code paths not exercised by that run
  need a direct old-vs-new bitwise comparison.
