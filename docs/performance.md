# Performance: build options, libraries, reproducibility

## Build options (CMake)

| Option | Default | Meaning |
|---|---|---|
| `SPODVR_USE_MKL` | `ON` | Use Intel oneMKL if found (`MKLConfig.cmake`; `MKLROOT` or `/opt/intel/oneapi/mkl/latest`). oneMKL then provides BLAS/LAPACK (the system libraries are not linked). Always the **sequential** oneMKL: the program calls it from its own OpenMP threads. |
| `SPODVR_FFT_LIBRARY` | `AUTO` | Implementation of the FFTW3 API: `AUTO` (oneMKL's FFTW3 interface if oneMKL is found, else FFTW), `MKL`, `FFTW`. With `FFTW` and oneMKL both present, FFTW is linked first so that the `fftw_*` symbols resolve to it. |
| `SPODVR_ARCH` | `x86-64-v3` | Value for `-march`. `x86-64-v3` (AVX2) runs on any recent x86-64 CPU; `x86-64-v4` (AVX-512) or `native` for a specific machine. Results of our own code do not depend on it (no FMA contraction). |
| `SPODVR_PGO` | `OFF` | Profile-guided optimisation, see below. |
| `SPODVR_FAST_MATH` | `OFF` | Experimental: `-ffast-math` instead of the flags below (FMA contraction, reassociation, vectorised math functions). Changes results at ~1e-12 on the gate grids and by 7.8e-12 on the thesis grid; passes the gates with the reaction-probability floor of 1e-11 except, with GCC 15, the hard-wall V<sub>asy</sub> table (1.6e-14 against 1e-14), see "Fortran vs C++" below. Native restart files are not interchangeable with the default build (README "Build variants"). Speed experiments only. |
| `SPODVR_POISON_UNINITIALISED` | `OFF` | Debug: fill uninitialised `AlignedVector` elements with NaN (any read before write then shows up as NaN). |
| `CMAKE_BUILD_TYPE=Release` | | `-O3` plus link-time optimisation. |

Floating-point flags (always): `-ffp-contract=off -fno-fast-math` (results comparable with the Fortran
code and independent of the CPU), plus `-fcx-fortran-rules -fno-math-errno -fno-trapping-math`, which
change no result but let GCC vectorise (complex multiplication without the C99 NaN-recovery call is
exactly the gfortran semantics). Never use `-ffast-math`/`-Ofast` or FMA contraction in production
builds: they change the results (measured below; `SPODVR_FAST_MATH` exists only for experiments).

### Profile-guided optimisation (about 15% faster)
```
cmake -S . -B build-pgo -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_PGO=GENERATE -DSPODVR_PGO_DIR=$PWD/pgo
cmake --build build-pgo
./build-pgo/spodvr <representative input.json> --threads <N>   # a short run (tens of steps) is enough
cmake -S . -B build-pgo -DSPODVR_PGO=USE
cmake --build build-pgo
```
Results are bitwise identical to the normal build.

### Running a oneMKL build
Inside the build tree the oneMKL library path is in the executable's RPATH. An installed or copied
executable needs the oneMKL libraries on the library path, e.g. `source /opt/intel/oneapi/setvars.sh`
or `LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib`. The program prints the libraries it uses in its
header line.

## Reproducibility
- Results never depend on the number of OpenMP threads (fixed work units, fixed summation order,
  identical per-thread FFT plans).
- **oneMKL code path (CNR)**: unless the environment variable `MKL_CBWR` is set, the program selects
  `MKL_CBWR_AUTO`: oneMKL uses the best code path for the CPU (e.g. AVX2 here, AVX-512 on newer
  machines) with fixed blocking and deterministic reductions. Results are then identical from run to
  run on the same machine, but may differ in the last bits between machines with different vector
  extensions. For identical results across machines set `MKL_CBWR=AVX2` (or `COMPATIBLE`), at some
  cost in speed on AVX-512 machines. The active branch is printed in the program header.
- Builds with different FFT libraries (oneMKL vs FFTW) or LAPACK libraries differ at the level of
  1e-14 relative; all regression gates pass with every combination.

## Memory and NUMA
- Large arrays (`AlignedVector`: wave function, potential grid, exp(−iV dt) table) are 64-byte aligned
  and are *not* zeroed on allocation; they are first written in parallel, per (Z, r) group, with the
  schedule of the loops that use them. On multi-socket (NUMA) machines this places every memory page
  on the socket of the thread that works on it ("first touch"). The first writes and the compute
  loops use the same (Z, r) groups and schedule (see "Threads and loop scheduling").
- exp(−iV dt) = (cos A, sin A) with A = −V dt is precomputed once if that complex table needs at most
  4 GiB (`PotentialOperator::default_table_limit`); otherwise cos/sin are evaluated in every step, as in
  the Fortran code (production grids: the potential alone is tens of GB, so the on-the-fly path is
  used and the memory is that of the real potential). Both paths give bitwise identical results.
- Work arrays per thread are cache-sized: (r, Z) units of ≤ 4 adjacent (irot, iy) values,
  (x, y) groups of r points ≤ 1 MiB; the DVR work array of one (Z, r) point (as in the Fortran code).
- The program prints a "MEMORY ESTIMATE" block before it computes the potential
  (`propagation/memory_estimate.*`): the size of each large array, the total while the propagation
  is set up (potential and exp(−iV dt) table exist together) and during the time loop, and the
  physical memory, with a warning if the run does not fit. Small tables and FFT plans are not
  counted; the medium grid estimate (5.0 GiB) matches the measured peak resident size (4.99 GiB).

## Threads and loop scheduling
- **Single NUMA node** (one socket, e.g. the Ryzen 5800X): dynamic schedule, and only as many threads
  as physical cores (`--threads 8` on the 8-core/16-thread 5800X). Work per (Z, r) point is not
  uniform (see "Fortran vs C++" below), and dynamic scheduling evens it out; hyperthreads add
  contention rather than speed. (An earlier measurement with the static schedule found 16 threads
  ~20% faster than 8, probably because the extra threads filled the idle time of the static split.)
- **Multiple NUMA nodes** (multi-socket): static schedule, and all available threads, bound to cores
  (`OMP_PROC_BIND=close OMP_PLACES=cores`, or `spread` to use all sockets' memory bandwidth). Static
  keeps every thread on the data it first touched (see "Memory and NUMA"); with a dynamic schedule a
  thread works on pages of the other socket.
- Implementation: `--schedule dynamic|static` (`src/core/loop_schedule.*`), set once at start-up
  with `omp_set_schedule` (`OMP_SCHEDULE` is ignored); the default is `static` with several NUMA
  nodes (`/sys/devices/system/node`), else `dynamic`; the header line prints the choice. Every loop
  that computes on, or first writes, the per-(Z, r) arrays runs over the same groups of
  `xy_group_points` consecutive r points with `schedule(runtime)`. These are the fused T_xy + V loop,
  `apply_xy` and `PotentialOperator::apply`, and the first writes of the wave function, the potential
  grid and the exp(−iV dt) table. So with `static,1` thread t owns groups t, t+T, ... in every one of
  them (first touch matches). Only chunk 1 is offered: `dynamic,1` (one NUMA node) and `static,1`
  (round robin, balanced to within one group; several NUMA nodes). The results are bitwise identical
  with either.
- Not covered: `apply_zr` (T_Z + T_r) works on units that span all (Z, r) points, so on a NUMA machine
  it reads memory of every socket whatever the schedule.
- To be investigated: the group size, larger chunks (the Fortran uses (nZ − nzpot)·nR/threads
  points), `static,1` vs `dynamic,1` on one node, and a default thread count of the physical cores.
  Plan: `~/work/loop_scheduling_optimization.md`.

## Benchmarking checklist
- Use a realistic grid: `tools/benchmark/run_benchmark.sh [build-dir] [threads] [steps]` (the
  "medium" grid below; ~5 GB of memory). The tiny test grids give misleading rankings (on them
  FFTW beats oneMKL; at realistic sizes it is the other way round or equal).
- **Rebuild every build directory you compare** (`cmake --build build-xxx`); a stale binary once led
  to a wrong conclusion here.
- Repeat each timing at least twice: on WSL2 the same run varies by 20% and occasionally by 2x.
- Fortran vs C++ (whole program): run them one after the other on an idle machine with the same
  `OMP_NUM_THREADS`, under `/usr/bin/time -f "%e %U"`. CPU time / wall time = threads busy on
  average: CPU time compares the work, the ratio shows idle threads (libgomp spin-waiting counts as
  CPU time in both codes, so the ratio is approximate). On WSL2 the program's own "Propagation time"
  (steady_clock) came out ~3.6% *larger* than `/usr/bin/time`'s wall time of the whole run
  (792 s vs 764 s): different clocks; compare like with like only.
- Compare the printed "Propagation time", not the wall time (which includes set-up and the
  potential grid).
- `spodvr_operator_timing <input.json>` (built in `build/tools/`) shows where the time goes.
- Check the results of every change: `tools/check_bitwise_refactor.sh` (must be bitwise identical)
  or `tools/compare_propagation.sh` (rounding-level changes, ≤ 1e-13), then the gates
  (`ctest`, and `SPODVR_LONG_TESTS=1 ctest -R Long` for the full 7000-step run).

## Measurements (AMD Ryzen 7 5800X, 16 threads, WSL2; propagation time only)
| Case | Baseline (faithful port) | Optimised, oneMKL FFT | Optimised, FFTW | Optimised + PGO (oneMKL) |
|---|---|---|---|---|
| hardwall_short grid (nx = ny = 4, JMax = 4; 100 steps) | 2.8 s | 2.4 s | 2.0 s | |
| medium (nx = ny = 16, JMax = MJMax = 8, 45 rot. states, 1.8 GB wave function; 20 steps) | 32.7 s | 19.0–20.0 s | 19.9 s | 16.9 s |

For the tiny test grids FFTW is faster (oneMKL has a large per-call overhead for 4×4 and length-10
transforms); at realistic sizes oneMKL is as fast or faster. Per time step on the medium case (16
threads): T_Z + T_r ≈ 0.4 s, T_xy ≈ 0.17 s, potential (incl. angular transforms) ≈ 0.33 s.

## Fortran vs C++ (full Dai–Light run, 8 threads, 2026-11-05)
Case `tests/reference/cu111_dailight_small`: 7000 steps on a 176 × 56 × 8 × 8 grid (nZ × nR × nX × nY)
with 15 rotational states (JMax = MJMax = 4; nθ = 3, nφ = 10), plus the on-the-fly analysis. AMD Ryzen 7
5800X (8 cores, Zen 3, no AVX-512), WSL2, `OMP_NUM_THREADS=8`, the three runs one after the other on an
otherwise idle machine; whole program, `/usr/bin/time`. Measured with GCC 13 / gfortran 13, before the
switch to GCC 15 (2026-11-06); not re-measured since.

| Binary | Wall time | CPU time (user) | Threads busy (CPU / wall) |
|---|---|---|---|
| Fortran (gfortran -O3, FFTW, no `-march`) | 680 s | 5255 s | 7.7 |
| C++, default flags (x86-64-v3, no FMA, LTO, oneMKL) | 764 s | 5040 s | 6.6 |
| C++, `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native`, PGO (trained on 200 steps of the same case) | 675 s | 4489 s | 6.7 |

- The C++ needs less CPU time than the Fortran (−4% with the default flags), but keeps the 8 threads
  less busy (6.6 of 8 against 7.7): it loses on wall time through load imbalance, not through slower
  kernels. At the Fortran's parallel efficiency the default build would take about 650 s.
- Main cause (from the code, not yet measured per thread): the fused T_xy + V loop
  (`apply_xy_and_potential`) runs over (Z, r) groups with `schedule(static)`. Here one group is a
  whole Z row, and only rows iz < nzpot (125 of 176; medium grid 113 of 176) do the angular transform
  and exp(−iV dt); the others are one scalar multiply. A static contiguous split gives threads 6–7
  almost only cheap rows. The Fortran (`SplitOperator.f90` KineticAndPotentialOp) uses
  `SCHEDULE(DYNAMIC, chunkSize)` for exactly this reason. Minor: the specular grid operator is serial.
- Since then this loop uses `schedule(dynamic, 1)` (bitwise identical results). On the grids of
  `inputs/thesis_fig5_3_*` (JMax 28, 435 rotational functions, 8 threads, 40 steps of 5 a.u.) the
  propagation time went from 176.9 to 154.5 s (DZ, −13%) and from 230.3 to 208.6 s (ER, −9%). The full
  Dai–Light benchmark above has not been re-run yet.
- fast-math + native + PGO: −12% wall, −11% CPU time against the default build. `-march=native` adds
  little on Zen 3 (no AVX-512: the same AVX2 as x86-64-v3); the gain is fast-math and PGO.
- fast-math changes results by ~1e-12 (norm 2e-12 relative, probabilities 1.5e-12 absolute against the
  default build; the limit for rounding-level optimisations is 1e-13). Gates: all regression gates of
  the regular suite and the long hard wall gate pass. The long Dai–Light gate failed on 2 of 1435
  probabilities with the then floor of 1e-12, both reaction probabilities (1 − sum): deviations
  2.0e-12 (E = 0.42 eV; default build 0.36e-12) and 5.3e-12 (E = 0.44 eV, printed half unit 5e-12;
  default build 4.1e-12). The floor was then raised to 1e-11 (reasoning in
  `tests/analysis_reference.hpp`), with which the fast-math build passes. Three unit tests with exact
  double comparisons fail by 1 ulp. Not the default build.
- Re-measured with GCC 15 (2026-11-07/08, [validation report](validation_thesis_fig5_3_dz/README.md)):
  - Both fast-math builds (`build-fast`, `build-cuda-fast`) fail gate 1's hard-wall V<sub>asy</sub> table by
    1.6e-14 against 1e-14; the default build passes. All other regular gates pass.
  - On the thesis grid (`inputs/thesis_fig5_3_dz`, 7000 steps, 435 rotational functions) the effect
    is larger: `build-cuda-fast` against `build-cuda` 7.8e-12 (norm-weighted relative S-matrix
    difference), 7.8× the figure above. Most of it is a uniform amplitude offset of 6.4e-12. Against
    the Fortran the default GPU build is 3.5e-13 off, the fast-math builds 7.5e-12 (CPU) and 7.9e-12 (GPU).
  - A default build refuses a fast-math build's native checkpoint and BK files (and vice versa): the
    stored effective `ZOptFactor` differs by 1 ulp and the reader compares exactly (README "Build variants").
- Earlier measurement with 16 threads (hyperthreads, includes test overhead): Fortran 637 s, C++ 816 s.

## CRP library (SystemID `CRPlibrary`, 2026-11-09)
Setting up the potential of a CRP data set costs microseconds per point, and it dominates the
set-up. Measurements on the synthetic data sets of the library tests: the gate cases
`tests/reference/crp_*_synthetic_small` with nX = nY = 12 and JMax = MJMax = 8 (nθ = 5, nφ = 18;
C4v with `nPhisOptimal` nφ = 20). That is 80 × 32 × 144 × 90 = 3.3·10⁷ PES evaluations, in
potential-only runs with 8 threads on an otherwise idle machine (GCC 15 / gfortran 15). The Fortran
time is the wall time of the whole run; the C++ time is the "6D potential computed in" line.

| Data set | Fortran | C++ | C++ with `UseSymmetry` + `nPhisOptimal` (points computed) | C++, one thread, per evaluation |
|---|---|---|---|---|
| C6v | 22.5 s | 17.4 s | 2.2 s (19 of 144 (x, y)) | 2.9 µs |
| C3v | 44.7 s | 26.7 s | 6.2 s (31 of 144) | 4.6 µs |
| C4v | 13.4 s | 12.0 s | 2.6 s (28 of 144) | 1.9 µs |

- Per evaluation: `build/tools/spodvr_crp_timing <data set>` (the QD wrapper at 98496 points
  spread over the cell, Z = −1..8, r = 0.6..4.1 and the orientations; best of 3).
- `UseSymmetry` (the Fortran switches it off for C3v and C4v) evaluates the PES only at the
  symmetry-unique (x, y) points; the potential differs from the full one by rounding only (self-check
  ≤ 6e-13 eV; [crp_library.md](crp_library.md)). Recommended for production runs.
- Profile (perf, one thread, `-O2` without LTO): the MURII cluster sum of the atom-surface PES takes
  45 % (C3v: 256 atoms per atom, most beyond the cut-off of the corrugation function), 11 % (C6v)
  and 13 % (C4v); the grid lookup of the molecular cuts (`locate_interval`, one per cut although
  the cuts share their grid) 9–15 %; libm (cos/sin of the site functions, acos/atan2) 12–22 %; the
  rest is spread over the interpolation and the derivatives (which the QD wrapper discards).

## What was optimised (and how it relates to the Fortran code)
| Change | Results |
|---|---|
| Compiler flags above, LTO, `-march` | bitwise identical |
| T_xy and V applied per group of (Z, r) points in one pass (as Fortran `KineticAndPotentialOp`) | bitwise identical |
| exp(−iV dt) precomputed (≤ 4 GiB) | bitwise identical |
| First-touch parallel initialisation | bitwise identical |
| T_Z + T_r as one 2D (r, Z) FFT on gathered contiguous units (Fortran: 1D FFTs in Z, then r) | ≤ 1e-14 |
| T_xy as one 2D (x, y) FFT per group of points, gathered with a single batch dimension | ≤ 1e-14 |
| Angular transform: all (x, y) columns of a point with one batched φ FFT | ≤ 1e-14 |
| Merged scalar factors (optical × rotation × gas phase) | ≤ 1e-14 |
| oneMKL BLAS/LAPACK and FFTs | ≤ 1e-14 |

Not done:
- FFTW wisdom / `FFTW_MEASURE`: oneMKL (the default) ignores planner flags, and measured plans would
  make results depend on the planning run.
- BLAS GEMM for the θ sums of the angular transform: the matrices are tiny for the test grids; the
  sums are written in a batched, GEMM-shaped form over the (x, y) columns, so a oneMKL/cuBLAS batched
  GEMM can replace them when profiling of production grids shows a need.

## GPU backend
The CUDA backend (`--device cuda`) is described, with its measurements, in [cuda.md](cuda.md).
