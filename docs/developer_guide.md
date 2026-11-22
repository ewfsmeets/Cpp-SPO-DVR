# Developer guide

How the C++ code is organised, the rules it follows, and what to touch (and which test notices a
forgotten step) for the common changes. The method itself, as the code implements it:
[method.md](method.md). Building, formatting and the test suite: [CONTRIBUTING.md](CONTRIBUTING.md).
The GPU backend has its own guide: [src/cuda/README.md](../src/cuda/README.md); the CRP library port:
[crp_library.md](crp_library.md).

## Module map

One line per directory of `src/`, with the Fortran files it replaces. Every `.cpp` file starts with a
comment naming its responsibility and, where there is one, its Fortran origin.

| Directory | What it holds | Fortran |
|---|---|---|
| `app/` | the program driver (`run.cpp`: input, grids, PES, gas phase, propagation, output), the command line, the backend factory (the only `#ifdef SPODVR_HAVE_CUDA`), the fine-analysis run | `Main.f90` |
| `input/` | the JSON reader and the checks made while reading, `validate_input`, the derived parameters (`build_config`), the plot-key checks | `InputNew.f90` (ReadInput, ValidateInput, ProcessInput), `Setup.f90` (scalar part) |
| `core/` | constants, the two unit sets, tolerant real comparisons, errors, aligned vectors, the loop schedule | `constants.f90`, `MyOperators.f90` |
| `grids/` | the Z, r, x, y grids with their momenta, the (θ, φ) DVR grid, Gauss–Legendre, associated Legendre functions, the rotational basis | `Setup.f90`, `setupdvr.f90`, `angulargrid.f90`, `legendre.f90` |
| `potential/` | the `Potential` interface and the PESs (hard wall, Dai–Light, Persson, the CRP adapter, the SCM correction in `scm/`), the factory, the potential on the grid, the surface-symmetry self-check, the 2D cut | `H2SurfaceSystem.f90`, `Potential.f90`, the PES modules, `DistortedSurfacePES.dir/` |
| `crplib/` | the CRP library port (standard library only, own constants and units) | `CRPlibrary.dir/` |
| `fortran_text/` | the list-directed reader shared by `crplib/` and the SCM code | gfortran `READ(unit,*)` |
| `symmetry/` | surface operations, point groups from generators, the lattice, the grid symmetry (representatives and orbits) | `C*vSurfaceSymmetry.f90` (generalised) |
| `gasphase/` | the Fourier grid Hamiltonian of the vibration and the rovibrational states | `FourierGridHamiltonian.f90`, `GasPhase.f90` |
| `propagation/` | the operator stages (kinetic, potential, specular, angular transform), their tables, the split-operator schedule, the CPU backend, the time loop, averages, the analysis surface, the initial packet, the memory estimate | `SplitOperator.f90`, `KineticOp.f90`, `Potential.f90`, `Absorption.f90`, `TwoGrids.f90`, `transform.f90`, `WaveFunction.f90`, `Averages.f90` |
| `cuda/` | the GPU backend (optional) | — |
| `analysis/` | the analysis basis, Balint-Kurti coefficients, the S matrix, the fine analysis, the norm budget (what each absorber removes per analysis window), the energy-resolved sums | `RotDifBasis.f90`, `Analyzor.f90`, `AsympInfo.f90` |
| `frames/` | the reductions of ψ for `dtFrames` and `FrameExtras` | — |
| `io/` | native checkpoint and BK files, the Fortran restart files, gfortran unformatted records, series directories | `WaveFunction.f90`, `AsympInfo.f90` |
| `fft/` | the FFT interface and its FFTW3-API implementation (oneMKL or FFTW), valid FFT sizes | `ffts.f90`, `fftw_2_scs.f90` |
| `linalg/` | the symmetric eigensolver interface and its LAPACK implementation | `jacobi.f90` (replaced by LAPACK) |
| `report/` | everything printed and every data file (numerical modules return data and print nothing) | the Fortran's output statements |

## The propagation in one page

- `WavePacketPropagation` (`propagation/propagation.cpp`) runs the time loop: per analysis step one
  call of the backend, then the analysis (BK coefficients, S-matrix sums), the averages (every dtAve)
  and the plot data. Observers (`PropagationObserver`) receive the analysis steps; `app/run.cpp`
  writes the checkpoint after the loop.
- `PropagationBackend` (`propagation/propagation_backend.hpp`) is the seam between the time loop and
  the stages: one virtual call per analysis step. `CpuBackend` holds ψ in host memory; `CudaBackend`
  in GPU memory.
- One analysis step is the list `analysis_step_schedule` (`propagation/split_operator_schedule.hpp`):
  T/2 V, (T V) n − 1 times, T/2. Both backends loop over the same list; what one T block and one V
  block do: [method.md](method.md#one-analysis-step).
- The stages: `ScatteringKineticOperator` (T_Z + T_r, T_xy), `PotentialOperator` (absorbers, rotation,
  the coupling to the specular grid, exp(−iV dt) in the DVR), `SpecularOperator` (the incoming channel
  on its own grid), `AngularTransform` (FBR ⇄ DVR). Each owns plain tables precomputed by
  `make_propagation_tables` and is applied to the whole wave function, OpenMP over fixed work units
  inside. The GPU backend has a twin of each, reading the same tables.
- Memory layouts (first index fastest, as in the Fortran): the scattering wave `(irot, iy + ny·ix,
  ir + nr·iz)`, the specular wave `(ir, iz)`, the potential `(iφ + nφ·iθ, ixy, ir + nr·iz)`. The GPU
  keeps ψ as `(izr, irot, ixy)`.
- The stored potential is −V·dt, and only for Z below z_asymptotic (`nz_potential` points).
- The CPU's angular transform runs over the (x, y) columns of a (Z, r) point in passes of
  `AngularTransform::columns_per_pass` columns (a constant byte limit, so the bits never depend on the
  machine), and T_xy runs on the block itself when a group holds one point, there is no momentum
  shift and the block is a multiple of 64 bytes (`ScatteringKineticOperator::transforms_in_place`).
- Results never depend on the thread count: every point is independent, every reduction has a fixed
  order, and work units have a fixed shape.

## Coding rules

The same rules as for the existing code; reviews check them.

- **Pointers:** only `std::unique_ptr` and `std::shared_ptr`. No raw pointers, owning or not, and no
  `new`/`delete`; non-owning access uses references, `std::reference_wrapper` or `std::span`. The
  exceptions are C-API boundaries, each kept local to its file: the FFTW calls in
  `fft/fftw_fft.cpp`, LAPACK in `linalg/lapack_eigensolver.cpp`, the allocator of `AlignedVector`,
  the `char*` calls of the binary file helpers (`io/binary_file.hpp`), `std::from_chars` in the
  list-directed reader, and device pointers in `src/cuda/` (owned by `DeviceBuffer`; kernels and
  launchers take raw device pointers).
- **No globals.** The configuration is an immutable `Config` built by `build_config` and passed
  explicitly.
- **Arrays:** multi-dimensional data are a flat `std::vector` (or `AlignedVector` for large arrays)
  with the index arithmetic in named helpers. Large arrays are not zeroed on allocation: every element
  is written, in parallel with the loops' own schedule (first touch), before it is read.
- **Constants and units:** one `Constants` object; the two Fortran unit sets `units::legacy` and
  `units::precise`, always passed explicitly, never a factor written into the code. `src/crplib/`
  keeps the CRP library's own units.
- **Readability:** code for a physicist who knows the Fortran. Short functions with one job; plain
  loops and simple data structures, no clever templates or member-pointer tricks; comments say what a
  block computes and why, not how the syntax works. A floating-point variable initialised from an
  `int` is declared `double`, never `auto` (`auto ord = order;` turns `amj / ord` into integer
  division).
- **Origins:** a ported routine names its Fortran origin in a short comment
  (`// Fortran: SplitOperator.f90 Propagate`); every `.cpp` file starts with its responsibility.
- **Numbers:** `double` and `std::complex<double>` everywhere; no FMA, no fast-math in production
  (`-ffp-contract=off -fno-fast-math`). A refactor meant to keep the numbers must keep them bitwise
  (`tools/check_bitwise_refactor.sh`; CONTRIBUTING, "Numerically neutral changes and rounding").
- **NaN fails every comparison:** reduce deviations with `nan_aware_max` and write checks as
  `EXPECT_LE(worst, tol)` ([validation.md](validation.md)).
- **Fortran semantics:** reproduce them unless a deviation is decided, and record every deviation in
  [fortran_deviations.md](fortran_deviations.md) with a test. The tolerant comparisons of
  `MyOperators.f90` are `fortran_compare::rgt` and friends (absolute 1e-12); use them only where the
  Fortran uses `.RGT.` and the like, plain operators where it uses `>`.

## Checklists

Each step names the test or check that notices when it is forgotten.

### A new potential energy surface

1. A class implementing `Potential` (`potential/potential.hpp`): `value()` in eV with the coordinate
   conventions stated there (skewed x, y; θ from the surface normal, φ from the first lattice vector),
   `system_parameters()`, `name()`. Copy an existing PES as a template.
2. Its `SystemID` in `make_potential` and `potential_is_implemented` (`potential/potential_factory.cpp`)
   — if it is missing there, any input with it fails with "unknown system ID".
3. Its data files, if any, in `potentials/` (tracked: the root `.gitignore` lets `potentials/` in although
   it ignores `*.dat` elsewhere).
4. A Fortran reference case in `tests/reference/` (`tools/fortran_reference/`, with a `MANIFEST.tsv`;
   for gate 1c also a `pes_cuts.json`) — `ReferenceManifest.FilesMatch` fails on a case without its
   manifest. Then list the case in the `INSTANTIATE_TEST_SUITE_P` lists of
   `tests/test_regression_pes.cpp` (`RegressionPes`, and `RegressionPesCut` for the cuts) and in
   `tools/fortran_reference/regenerate_all.sh`: nothing notices a case missing there, and a missing
   reference file makes a regression test skip (green) unless the case is in
   `required_reference_cases`. Name it in [validation.md](validation.md) (`Docs.Consistency` checks the
   gate cases). A PES without a Fortran counterpart has no gate 1: say what replaces it.
5. With `UseSymmetry` the PES self-check (`SymmetryTolerance`) runs at set-up; test the symmetry like
   the existing PESs (`tests/test_potential_symmetry.cpp`).
6. The wave-packet designer's system parameters per SystemID
   (`tools/wave_packet_designer/program/potential_systems.py`): nothing notices it missing.
7. Docs: `SystemID` in [input_schema.md](input_schema.md), the row in [port_status.md](port_status.md),
   the README table — nothing checks the SystemID list; check it by hand.

### A new input key

1. Add the field, with its Fortran default (`globaldata.f90`), to `InputParameters` (`input/config.hpp`);
   derived quantities go into `build_config` (`input/process.cpp`).
2. Read it in `input/json_input.cpp` (optional keys keep existing inputs valid) and check it where the
   Fortran does (while reading, or in `input/validate.cpp`) — a violating row with the check's message in the
   `Checks/InputViolation` table of `tests/test_input.cpp`; an unknown key is an error, so a key
   that is never read fails every input that uses it.
3. Document it in [input_schema.md](input_schema.md) — `Docs.Consistency` fails on an undocumented key.
4. If a restart must not change it: add it to `restart_parameters` (`io/checkpoint.cpp`) —
   `Restart.MismatchingInputIsAnError` shows the pattern.
5. Add it to the wave-packet designer's key table (`tools/wave_packet_designer/program/json_input.py`)
   — `WavePacketDesigner.InputReading` fails otherwise — and port its checks there if the designer
   uses it (the `Shadows:` lines name the C++ files).
6. If it switches an addition to the time loop: free when off — `TimeLoopCalls.*` count the backend
   calls (CONTRIBUTING, "Additions to the time loop").

### A new large array

1. Allocate it as `AlignedVector`, first written in parallel with the schedule of the loops that use
   it.
2. Add it to the memory estimate (`propagation/memory_estimate.cpp`, and `estimate_gpu_memory` in
   `cuda/cuda_backend.cpp` for the GPU) — `MemoryEstimate.MatchesAllocations` and
   `CudaBackend.MemoryEstimateMatchesTheAllocations` compare the estimates with the allocations; run
   the tests of the `build-poison` variant (`SPODVR_POISON_UNINITIALISED`): a read before the first
   write shows up as NaN.
3. Follow in the designer's copy (`tools/wave_packet_designer/program/memory_estimate.py`) and
   regenerate `expected_cases.json` — `MemoryEstimate.MatchesTheWavePacketDesignerCases` (and its
   `CudaBackend` twin) fail until it does.

### A new output file

1. Write it in `report/` (numerical code returns the data).
2. Document it in [output_files.md](output_files.md) — `Docs.Consistency` checks the file list.
3. If the plot tool reads it: its reader in `tools/plots/readers/` and the fixture runs
   (`PlotRuns.*`).

### Porting a Fortran routine

1. Cite the routine in a short comment; keep the Fortran's order of operations where it decides the
   bits.
2. Check its quad-precision literals and expressions (`-fdefault-real-8` promotes `…d0` and `DBLE`)
   and its `**` with integer exponents (`fortran_integer_power`).
3. Compare with a 5.24 reference built with `-O3` (`tools/fortran_reference/`), never the upstream
   default flags.
4. Record any deviation in [fortran_deviations.md](fortran_deviations.md) with a test; a new gate
   or a looser tolerance needs a negative test that makes it fail.
