# Performance: build options, libraries, reproducibility

## Build options (CMake)

| Option | Default | Meaning |
|---|---|---|
| `SPODVR_USE_MKL` | `ON` | Use Intel oneMKL if found (`MKLConfig.cmake`; `MKLROOT` or `/opt/intel/oneapi/mkl/latest`). oneMKL then provides BLAS/LAPACK (the system libraries are not linked). Always the **sequential** oneMKL: the program calls it from its own OpenMP threads. |
| `SPODVR_FFT_LIBRARY` | `AUTO` | Implementation of the FFTW3 API: `AUTO` (oneMKL's FFTW3 interface if oneMKL is found, else FFTW), `MKL`, `FFTW`. With `FFTW` and oneMKL both present, FFTW is linked first so that the `fftw_*` symbols resolve to it. |
| `SPODVR_ARCH` | `x86-64-v3` | Value for `-march`. `x86-64-v3` (AVX2) runs on any recent x86-64 CPU; `x86-64-v4` (AVX-512) or `native` for a specific machine. Results of our own code do not depend on it (no FMA contraction). |
| `SPODVR_PGO` | `OFF` | Profile-guided optimisation, see below. |
| `SPODVR_POISON_UNINITIALISED` | `OFF` | Debug: fill uninitialised `AlignedVector` elements with NaN (any read before write then shows up as NaN). |
| `CMAKE_BUILD_TYPE=Release` | | `-O3` plus link-time optimisation. |

Floating-point flags (always): `-ffp-contract=off -fno-fast-math` (results comparable with the Fortran
code and independent of the CPU), plus `-fcx-fortran-rules -fno-math-errno -fno-trapping-math`, which
change no result but let GCC vectorise (complex multiplication without the C99 NaN-recovery call is
exactly the gfortran semantics). Never use `-ffast-math`/`-Ofast` or FMA contraction: they change the
results.

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
  and are *not* zeroed on allocation; they are first written in parallel with the static OpenMP
  schedule of the loops that use them. On multi-socket (NUMA) machines this places every memory page
  on the socket of the thread that works on it ("first touch").
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

## Threads
- Use all hardware threads: on the 8-core/16-thread Ryzen 5800X, 16 threads are ~20% faster than 8.
- On multi-socket machines bind threads to cores, e.g. `OMP_PROC_BIND=close OMP_PLACES=cores` (or
  `spread` to use all sockets' memory bandwidth); with first touch this keeps data local. On a single
  socket pinning made no difference.

## Benchmarking checklist
- Use a realistic grid: `tools/benchmark/run_benchmark.sh [build-dir] [threads] [steps]` (the
  "medium" grid below; ~5 GB of memory). The tiny test grids give misleading rankings (on them
  FFTW beats oneMKL; at realistic sizes it is the other way round or equal).
- **Rebuild every build directory you compare** (`cmake --build build-xxx`); a stale binary once led
  to a wrong conclusion here.
- Repeat each timing at least twice: on WSL2 the same run varies by 20% and occasionally by 2x.
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

## Notes for a future CUDA implementation
- The stage classes (`ScatteringKineticOperator`, `PotentialOperator`, `SpecularOperator`,
  `AngularTransform`) work on whole flat arrays; a GPU backend replaces a stage (or the fused
  T_xy + V stage) as a whole.
- All FFT layouts have one batch dimension (`FftLayout`, `FftLayout2D`) and map onto
  `cufftPlanMany`; oneMKL's FFTW3 interface has the same restriction, which is why the (x, y)
  groups are gathered first.
- Tables are plain arrays; reductions use a fixed order (on a GPU this needs a deterministic
  reduction scheme to keep thread-count independence).
