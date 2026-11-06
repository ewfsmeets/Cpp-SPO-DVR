# GPU backend (CUDA)

`spodvr --device cuda` propagates the wave packet on an NVIDIA GPU. The wave packet stays in GPU
memory for the whole run. After every analysis step only a few small pieces are copied back to the
host:
- the specular wave;
- the nr (Z, r) blocks at Zana (the input of the BK analysis);
- every `dtAve`, one density value per (Z, r) point.

Everything else is unchanged: set-up, potential, gas phase, analysis, S-matrix, output files,
checkpoints and restarts (also from or into Fortran files, and between CPU and GPU runs).

## Building
```bash
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_USE_CUDA=ON
cmake --build build-cuda
OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R Cuda --output-on-failure
```
| Option | Default | Meaning |
|---|---|---|
| `SPODVR_USE_CUDA` | `OFF` | build the backend. Without it nothing CUDA-related is compiled, and `--device cuda` stops with a clear error |
| `CMAKE_CUDA_COMPILER` | nvcc from `PATH`, else `/usr/local/cuda/bin/nvcc` | the CUDA compiler |
| `CMAKE_CUDA_ARCHITECTURES` | `native` (the GPU of the build machine) | set e.g. `89-real` (RTX 40xx), `80-real` (A100), `90-real` (H100), or a list, for a build that runs elsewhere |
| `CMAKE_CUDA_HOST_COMPILER` | `$CUDAHOSTCXX`, else the C++ compiler of the build | host compiler of nvcc. Each CUDA release accepts GCC only up to a limit (`crt/host_config.h`) |
| `SPODVR_CUDA_SYNC_DEBUG` | `OFF` | wait for and check every kernel right after its launch (errors are reported where they happen) |

- **Toolchain:** one compiler for host code and kernels. Tested with CUDA 13.3 (NVIDIA's apt
  repository, see the README) and GCC 15 (needs CUDA ≥ 13.3). Developed with CUDA 12.0 (Ubuntu 24.04 packages, g++-12
  as nvcc's host compiler, C++17 kernels); untested with 12.0 since the kernels are C++20.
- **Tests:** the GPU tests in `tests/cuda/` are only built with the option, and they skip on a
  machine without a GPU.
- **Host code:** only kernels and their launch functions are in `.cu` files (C++20, nvcc). All host
  code of the backend is in `.cpp` files compiled by the normal C++ compiler: they keep nvcc away
  from the project headers and the host code inside LTO.
- **Compile flags:** `spodvr_warnings` and `spodvr_fp` apply to C++ only. The kernels get
  `-fmad=false` (no FMA contraction, as on the CPU). With `SPODVR_FAST_MATH=ON` they get
  `--use_fast_math --extra-device-vectorization` instead (see "Fast build" under Performance).
  Both builds add `-lineinfo` (source lines for Nsight; the generated SASS is identical).
- **cuFFT callback:** `src/cuda/fft_callbacks.cu` is not part of the library. CMake builds it in
  three steps:
  - nvcc compiles it to LTO-IR (`-gencode arch=compute_75,code=lto_75 -dc -fatbin`; compute_75 IR
    runs on every GPU from sm_75 on);
  - `cmake/embed_binary.cmake` turns that into a byte array (`build-*/generated/`);
  - `src/cuda/fft_callbacks.cpp` hands the array to cuFFT, which links it into its FFT kernel when
    the plan is made (`cufftXtSetJITCallback`).

  At run time cuFFT needs `libnvJitLink` and `libnvrtc` of its own CUDA release on the library
  path; both are part of the toolkit. Without them the GPU run stops with an error that names
  them.

  The callback must:
  - be at global scope and not `extern "C"`; otherwise cuFFT fails with
    `CUFFT_NVJITLINK_FAILURE`;
  - multiply with `__dmul_rn`/`__dadd_rn`, because cuFFT sets the link options and nvJitLink
    contracts to FMA by default.

## Running
```bash
./build-cuda/spodvr input.json --device cuda --threads 8 --output-dir run
```
- **Before the potential is computed:** the program prints the GPU and the CUDA library versions,
  plus a GPU MEMORY ESTIMATE. It stops if the free GPU memory is too small, so a run never fails
  after a long set-up.
- **Threads:** `--threads` still applies to the host work (potential, gas phase, BK analysis).
- **Run modes without a propagation:** `--device cuda` has no effect (potential-only, analysis-only).
- **Default device:** `--device cpu`. Results never depend silently on the hardware.

## Design
`PropagationBackend` (`src/propagation/propagation_backend.hpp`) is the boundary. The time loop
(`WavePacketPropagation::run`) makes one virtual call per analysis step, then reads what the
analysis needs:
- `CpuBackend` wraps the OpenMP operator stages;
- `CudaBackend` (`src/cuda/cuda_backend.cpp`) holds the device stages;
- `src/app/backend_factory.cpp` chooses one of them, and is the only file that knows whether the
  CUDA backend was built.

### Device layout
On the GPU, ψ is stored as `(izr, irot, ixy)` with izr = ir + nr·iz fastest (`DeviceLayout`). Host
files and the CPU code keep `(irot, ixy, izr)`. With izr fastest, every stage is either one batched
cuFFT with one batch dimension or a kernel whose neighbouring threads read neighbouring memory, and
none of the CPU's gather/scatter code is needed:

| Stage | GPU | CPU twin |
|---|---|---|
| T = T_Z + T_r + T_xy, one stage (methodology exception, below) | momentum shift (`multiply_by_blocked_factors`); one 2D (Z, r) cuFFT over all n_rot·nxy contiguous blocks; one strided 2D (x, y) cuFFT, batches at distance 1; `multiply_by_kinetic_factors` (both factor sets in one pass); both inverses; unshift | two stages: `ScatteringKineticOperator::apply_zr` (gathers ≤ 4 values per unit), then `apply_xy_group` (gathers 1 MiB groups, fused with V) |
| V, Z ≥ zasymp | `apply_gas_phase_step` (elementwise) | `PotentialOperator::apply_block` |
| V, Z < zasymp | per slab of (Z, r) points, see below | `apply_block`, `apply_interaction`, `AngularTransform` |
| specular grid | 2D cuFFT, multiply, inverse, `apply_specular_gas_phase_and_absorbers` | `SpecularOperator` |

Kernels of the interaction region, one slab at a time. The slab is the largest divisor of the
number of (Z, r) points whose DVR scratch fits in 1 GiB:
1. `absorb_and_transform_to_dvr`: one thread per (Z, r, x, y) column. Absorption, rot(dt/2) and the
   specular wave added to the initial channel, then the Legendre sums FBR → DVR.
2. cuFFT over φ: length nφ, batch columns·nθ at distance 1. The scratch layout is (col, θ, φ),
   without the CPU's nφ+1 padding.
3. exp(iA) of the column's (x, y) representative at the mapped φ (surface symmetry). It is applied
   by a cuFFT load callback of the inverse φ FFT (`src/cuda/fft_callbacks.cu`), so the DVR scratch
   is not read and written once more. The kernel `multiply_by_exp_potential` does the same when
   there is no φ FFT (nφ = 1), and serves as the reference in the tests.
4. Inverse φ cuFFT, then `transform_to_fbr_and_rotate`: Legendre sums DVR → FBR, then rot(dt/2).
5. `project_onto_specular`: one thread per (Z, r) point; the sum over (x, y) runs in the CPU order.

The tables are the CPU's. `potential_step_factors`, `specular_gas_phase_factors` and the
`AngularTransform` matrices are shared, and exp(iA) is computed on the host with the CPU's
`exp_i`, then re-laid out with izr fastest.

### One method, two executions
The CPU and the GPU backend implement the same method:
- the same operator sequence `T/2 V (T V)^(n-1) T/2`;
- the same stages within T and V (table above);
- the same tables;
- the same order of operations per element.

They differ only in how the stages are executed: memory layout, gather/scatter, the CPU's
fusion of T_xy and V per group, and the GPU's slabs. This keeps every GPU stage traceable to its
CPU twin, and it is what the per-stage CPU-vs-GPU tests in `tests/cuda/` rely on.

**Rule for changes to the GPU backend:**
- **Allowed in the GPU backend alone:** a change to *how* a stage runs. Examples: layout, fusing
  adjacent elementwise work into a kernel or a cuFFT callback, streams and overlap, launch shapes.
- **Never in the GPU backend alone:** a change to *what* is computed: the operator
  decomposition, the order of the factors per element, the tables, or the order of a reduction
  (the inside of an FFT excepted). Such a change goes into both backends in the same commit, or
  into neither.
- A deliberate exception is listed below, with its reason and the test that compares it with
  the CPU.

**Methodology exceptions:**
1. **Merged kinetic stages** (`DeviceKineticOperator::apply`). The GPU applies T_Z + T_r and T_xy
   as one stage:
   - shift → (Z, r) and (x, y) cuFFTs → one multiplication with both factors → inverses → unshift;
   - the CPU applies two stages, T_Z + T_r, then T_xy (with V).

   Why it is equivalent: T_Z + T_r and T_xy commute, and so do the momentum shift (a factor per
   ixy) and the (Z, r) transform. Per element the multiplication keeps the CPU's order: the
   (Z, r) factor with its scales first, then the (x, y) factor with its scales. Only the FFT
   order differs.

   Why the GPU does it: one pass over ψ less, −29 ms per DZ time step (−4.4%). The CPU keeps
   two stages, because its cache-blocked T_Z + T_r and its T_xy fused with V per group would be
   lost by the merge.

   Test: `CudaBackend.MergedKineticStageMatchesTheCpuStages` (≤ 1e-13, full and half step, with
   and without the momentum shift, hard wall and Dai-Light grids). Against the GPU before the
   merge:
   - averages ≤ 1e-14 relative;
   - BK coefficients ≤ 6e-13 of max |c|;
   - all GPU gates pass, the long ones included.

## Numerics
- **No FMA:** `-fmad=false` everywhere. Complex products are formed as GCC does with
  `-fcx-fortran-rules` (`Complex` in `src/cuda/device_complex.hpp`), and every kernel keeps the CPU
  operation order. For the same input, the elementwise kernels, the Legendre sums and the specular
  projection are bitwise equal to the CPU stages.
- **FFT scale factors:** cuFFT does not scale. The scale factors are applied in the multiplication
  kernel between the transforms, as two separate multiplications. Folding their product into the
  factor table rounds it once more, and that bias repeats in every step: the norm drifted from the
  CPU's by 2.6e-13 relative in 1000 steps (hardwall_short). Applied separately, the drift is halved.
- **What remains is the FFT library:** cuFFT rounds differently from oneMKL/FFTW. After 1000 steps
  of hardwall_short, the CPU and GPU averages differ by ≤ 1.0e-13 relative (norm), against
  1.2e-13 between the CPU with oneMKL and the CPU with FFTW. This is the same level the CPU build
  already accepts between FFT libraries (docs/validation.md: within ~1e-13, all gates pass).
- **Densities** use CUDA's `hypot`, which can differ from glibc's by an ulp.
- **Determinism:** no atomics, fixed launch shapes, everything on the default stream, so two GPU
  runs are bitwise identical. That holds for the same GPU model and library versions only, so no
  GPU reference files are stored.
- **Restart:** a GPU run split by a restart (native or Fortran files) is bitwise identical to the
  unsplit GPU run.

## Memory
GPU memory estimate for the DZ benchmark (below):

| Array | Size |
|---|---|
| ψ | 4.8 GiB |
| exp(−iV dt) table (19 symmetry representatives of 144 (x, y) points) | 1.0 GiB |
| DVR scratch (slab of 500 points) | 1.0 GiB |
| cuFFT work area (limit; the DZ grid needs < 1 MiB with cuFFT 12.3 / CUDA 13.3) | ≤ 0.5 GiB |
| BK plane, copy buffer | small |

The CUDA context needs a few hundred MiB on top (`CUDA_MODULE_LOADING=LAZY` is set at start-up).
- **Host memory:** the host allocates ψ only to read a restart or to write the checkpoint at the end.
  The host MEMORY ESTIMATE block still describes a CPU run.
- **Production grids** (> 64 GB) do not fit on one GPU of this class. The options are surface
  symmetry (C6v: the potential ~12× smaller), a data-centre GPU (80–141 GB), or a multi-GPU version
  (not implemented).

## Performance
Measure with `build-cuda/tools/spodvr_cuda_operator_timing <input.json> [repetitions]`. It times:
- each stage and each kernel of the potential stage;
- a whole analysis step as the program runs it. That step has one kinetic stage more than n
  time steps (`T/2 V (T V)^(n-1) T/2`), so its time per step is what the full program reaches
  (DZ: 625 ms against 588 ms for the sum of the stages).

Profiling:
- **Nsight Systems 2026.1** (`/usr/local/cuda/bin/nsys`) works on this WSL2 machine, for CUDA
  API and kernel traces:
  ```bash
  nsys profile --trace=cuda --sample=none --cpuctxsw=none -o trace ./build-cuda/spodvr in.json --device cuda
  nsys stats --report cuda_gpu_kern_sum trace.nsys-rep
  ```
  A DZ trace showed the GPU idle for less than 1 ms per time step, so the host analysis never
  holds up the GPU.
- **ncu** needs access to the GPU performance counters, which this WSL2 machine does not grant.

DZ benchmark (`inputs/thesis_fig5_3_dz`: 435 rotational states, 128 × 40 × 12 × 12, ψ 4.8 GiB,
19 C6v representatives), RTX 4070 Ti against 8 threads of a Ryzen 7 5800X, per full time step
(10 repetitions; single runs on this WSL2 GPU vary by up to ~30%, it also drives the desktop;
measured before the kinetic stages were merged, see "Measured optimisations" below):

| Stage | GPU | CPU (`spodvr_operator_timing`) |
|---|---|---|
| T_Z + T_r (GPU: 2 cuFFTs 126 ms + 1 multiply pass) | 150 ms | 1044 ms |
| T_xy (GPU: 2 cuFFTs 132 ms + 1 multiply pass) | 158 ms | 549 ms |
| Potential | 301 ms | 1952 ms |
| One time step | **609 ms** | 3567 ms (fused T_xy + V) |

GPU potential stage by kernel, summed over the 8 slabs:

| Kernel | Time |
|---|---|
| Legendre sums FBR → DVR | 77 ms |
| φ FFTs | 2 × 48 ms |
| exp(iA) | 59 ms |
| Legendre sums DVR → FBR | 68 ms |
| gas phase + projection | 6 ms |

Whole program, 80 steps: 52.1 s of propagation time (CPU: 305.0 s; the averages of the two runs
agree to 1.0e-14 relative, the BK coefficients to 1.8e-15 absolute) (0.65 s per step) against ≈ 3.6 s per step on
the CPU: **5.5×**. The 10 000 steps of the benchmark would take ≈ 1.8 h instead of ≈ 10 h.

What limits the speed on this card:
- **Memory passes:** one pass over ψ (read + write) takes 25.5 ms (400 GB/s). A 2D (Z, r) cuFFT
  takes 2 passes, a strided 2D (x, y) cuFFT 2.7 (a contiguous one would take 2), and the
  multiplication kernels run at full bandwidth.
- **FP64 instructions:** the angular transform needs ≈ 17 G FP64 operations per direction per step.
  GeForce cards execute FP64 at 1/64 of the FP32 rate, and without FMA (`-fmad=false`, required for
  CPU-comparable results) every multiply and add is a separate instruction: ≈ 0.3 T per second here.
  So the two Legendre kernels are compute-bound at 65–77 ms.
- **Tried without gain:**
  - Keeping the per-thread coefficient arrays in registers (templated sizes, unrolled loops). It
    changed nothing and was reverted.
  - Starting the Legendre sums from their first term instead of 0 + … (2026-11-06; CPU and GPU,
    byte-identical results). It removes ≈ 6% of the FP64 additions of FBR → DVR, yet the kernel
    times did not change (70.6 → 70.7 ms and 66.3 → 67.7 ms), so it was reverted. The count of
    FP64 additions is evidently not what limits these kernels, although FMA (fast build) does
    speed them up.

Measuring a change: compare with **interleaved A/B runs** (base and new build alternately,
`spodvr_cuda_operator_timing`, "one analysis step … per time step"). Runs hours apart differ by
≈ 5% on this machine: the GPU also drives the desktop, and clocks vary. Interleaved runs repeat
within 0.5 ms. The isolated stage lines of the timing tool depend on what ran before them in
the same process; use them as a breakdown, not for comparisons between builds.

Measured optimisations (DZ, strict build, interleaved A/B, per time step):

| Change | Before | After |
|---|---|---|
| Merged kinetic stages (2026-11-06; one multiplication pass instead of two: 26.4 ms instead of 2 × 24.8 ms) | 656.6 ms | 627.5 ms (−4.4%) |
| exp(iA) as cuFFT load callback of the inverse φ FFT (2026-11-06; kernel 58.5 ms gone, φ FFT 48 → 57 ms; this A/B drifted, so pairs with the base run after the new one: −45 and −44 ms) | 616.9 / 598.7 ms | 571.6 / 554.9 ms (≈ −7%) |

Next steps, to be judged on the production GPU (on an A100/H100 the kernels become memory-bound,
and fewer passes pay off directly):
- **More cuFFT callbacks:** the kinetic multiplication (`multiply_by_kinetic_factors`) as the
  store callback of the forward (x, y) FFT, or as the load callback of the inverse one. That saves
  one more pass over ψ.
- **Even/odd split of the θ sums** (heteronuclear molecules only): the θ grid is symmetric about
  π/2 and P_J^m(−x) = ±P_J^m(x), so the split halves the FP64 work of the Legendre sums. For
  homonuclear molecules the grid already covers only [0, π/2] and every mJ block has a single
  parity, so there is nothing to gain there. It changes the summation order: CPU and GPU together.
- **A fused angular kernel** with an in-kernel φ transform: removes the DVR scratch and its four
  passes, but needs more FP64 work for the DFT.

### Fast build (`build-cuda-fast`: FMA in the kernels)
```bash
cmake -S . -B build-cuda-fast -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_USE_CUDA=ON \
      -DSPODVR_FAST_MATH=ON -DSPODVR_ARCH=native
```
`SPODVR_FAST_MATH` also switches the kernels to `--use_fast_math --extra-device-vectorization`.
For our FP64 kernels only the FMA contraction (`-fmad=true`) matters; the rest of
`--use_fast_math` (flush-to-zero, approximate division, square root and intrinsics) acts on FP32
only. The host code gets `-ffast-math -march=native`. No PGO: in a GPU run the host code is a small
part of the time. (With CUDA 12.0 the PGO training failed because GCC 13 objects were mixed with
nvcc's g++-12 host objects; with one compiler for both this cause is gone, but it was not retried.)

Measured on the RTX 4070 Ti (2026-11-06):

| | `build-cuda` | `build-cuda-fast` |
|---|---|---|
| DZ, Legendre sums FBR → DVR / DVR → FBR | 68 / 64 ms | 38 / 50 ms |
| DZ, one time step (timing tool, 10 repetitions) | 571 ms | 536 ms (−6%) |
| DZ, 200 steps, propagation time | 124.9 s | 115.0 s (−8%) |
| cu111_dailight_small (JMax 4), 7000 steps, propagation time | 109.3 s | 108.6 s |

The gain is in the Legendre sums only (FFTs and multiplication passes are memory-bound); it is
larger for larger JMax and smaller on A100/H100, where FP64 is not the limit.

Results, fast against strict build:
- All GPU tests pass, including the Fortran gates 3 and 4 (also the long ones) and the CPU-vs-GPU
  analysis-step tests (≤ 1e-13).
- Full cu111_dailight_small (7000 steps): norm 2.0e-12 relative, rest norm 3.9e-12, ⟨Z⟩ 2.2e-14,
  ⟨r⟩ 4.3e-15; S-matrix elements 1.0e-12 absolute, probabilities 1.6e-12 absolute. The same level
  as the CPU fast-math build (docs/performance.md).
- DZ, 200 steps: averages ≤ 2.7e-14 relative, S ≤ 3e-14 absolute.

So the fast build is within the regression tolerances, but no longer comparable to the CPU build
to the last bits; use the strict build when comparing with the CPU or the Fortran code.

## Files
| File | Content |
|---|---|
| `src/cuda/cuda_backend.{hpp,cpp}` | `CudaBackend`, the factory, device info, GPU memory estimate |
| `src/cuda/device_kinetic_operator.{hpp,cpp}` | the kinetic stage T (merged) and the specular grid (plans, tables, launch order) |
| `src/cuda/device_potential_operator.{hpp,cpp}` | the potential stage (tables, exp(iA), slab loop) |
| `src/cuda/*_kernels.cu`, `kernels.hpp` | the kernels and their launch functions |
| `src/cuda/cuda_fft.{hpp,cpp}` | cuFFT plans (64-bit, slabbed batches, shared work area, optional LTO load callback) |
| `src/cuda/fft_callbacks.{hpp,cu,cpp}`, `cmake/embed_binary.cmake` | the exp(iA) load callback (LTO-IR, embedded at build time) |
| `src/cuda/device_buffer.{hpp,cpp}`, `device_complex.hpp`, `cuda_error.{hpp,cpp}`, `launch.cuh` | device memory, complex arithmetic, layout, errors, launch helpers |
| `tests/cuda/` | cuFFT vs FFTW, CPU vs GPU per analysis step, determinism, memory estimate, restart, gates 3/4 |
| `tools/benchmark/cuda_operator_timing.cpp` | per-stage and per-kernel timing |
