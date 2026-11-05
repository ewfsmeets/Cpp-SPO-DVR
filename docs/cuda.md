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
| `CMAKE_CUDA_ARCHITECTURES` | `89-real` | GPU generation (RTX 40xx). Use e.g. `80-real` (A100) or `90-real` (H100), or a list |
| `CMAKE_CUDA_HOST_COMPILER` | `/usr/bin/g++-12` if present | host compiler of nvcc. CUDA 12.0 rejects GCC 13 (`crt/host_config.h`); CUDA ≥ 12.4 accepts it |
| `SPODVR_CUDA_SYNC_DEBUG` | `OFF` | wait for and check every kernel right after its launch (errors are reported where they happen) |

- **Toolchain:** CUDA 12.0 from the Ubuntu packages (nvcc, cuFFT 11.0).
- **Tests:** the GPU tests in `tests/cuda/` are only built with the option, and they skip on a
  machine without a GPU.
- **Host code:** only kernels and their launch functions are in `.cu` files (C++17, nvcc). All host
  code of the backend is in `.cpp` files compiled by the normal C++ compiler. nvcc 12.0 cannot
  compile some C++20 constructs of the project headers (parenthesised aggregate initialisation).
- **Compile flags:** `spodvr_warnings` and `spodvr_fp` apply to C++ only. The kernels get
  `-fmad=false` (no FMA contraction, as on the CPU).

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
| T_Z + T_r | one 2D (Z, r) cuFFT over all n_rot·nxy contiguous blocks, `multiply_by_periodic_factors`, inverse | `ScatteringKineticOperator::apply_zr` (gathers ≤ 4 values per unit) |
| T_xy | one strided 2D (x, y) cuFFT, batches at distance 1; `multiply_by_blocked_factors` (also the momentum shift) | `apply_xy_group` (gathers 1 MiB groups) |
| V, Z ≥ zasymp | `apply_gas_phase_step` (elementwise) | `PotentialOperator::apply_block` |
| V, Z < zasymp | per slab of (Z, r) points, see below | `apply_block`, `apply_interaction`, `AngularTransform` |
| specular grid | 2D cuFFT, multiply, inverse, `apply_specular_gas_phase_and_absorbers` | `SpecularOperator` |

Kernels of the interaction region, one slab at a time. The slab is the largest divisor of the
number of (Z, r) points whose DVR scratch fits in 1 GiB:
1. `absorb_and_transform_to_dvr`: one thread per (Z, r, x, y) column. Absorption, rot(dt/2) and the
   specular wave added to the initial channel, then the Legendre sums FBR → DVR.
2. cuFFT over φ: length nφ, batch columns·nθ at distance 1. The scratch layout is (col, θ, φ),
   without the CPU's nφ+1 padding.
3. `multiply_by_exp_potential`: exp(iA) of the column's (x, y) representative at the mapped φ
   (surface symmetry).
4. Inverse φ cuFFT, then `transform_to_fbr_and_rotate`: Legendre sums DVR → FBR, then rot(dt/2).
5. `project_onto_specular`: one thread per (Z, r) point; the sum over (x, y) runs in the CPU order.

The tables are the CPU's. `potential_step_factors`, `specular_gas_phase_factors` and the
`AngularTransform` matrices are shared, and exp(iA) is computed on the host with the CPU's
`exp_i`, then re-laid out with izr fastest.

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
  already accepts between FFT libraries (CLAUDE.md: within ~1e-13, all gates pass).
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
| cuFFT work area (limit; the DZ grid needs 34 MiB) | ≤ 0.5 GiB |
| BK plane, copy buffer | small |

The CUDA context needs a few hundred MiB on top (`CUDA_MODULE_LOADING=LAZY` is set at start-up).
- **Host memory:** the host allocates ψ only to read a restart or to write the checkpoint at the end.
  The host MEMORY ESTIMATE block still describes a CPU run.
- **Production grids** (> 64 GB) do not fit on one GPU of this class. The options are surface
  symmetry (C6v: the potential ~12× smaller), a data-centre GPU (80–141 GB), or a multi-GPU version
  (not implemented).

## Performance
Measure with `build-cuda/tools/spodvr_cuda_operator_timing <input.json> [repetitions]`, which
times each stage and each kernel of the potential stage. ncu needs access to the GPU performance
counters (not granted on this WSL2 machine), and the nsys of the Ubuntu packages cannot import its
own traces.

DZ benchmark (`inputs/thesis_fig5_3_dz`: 435 rotational states, 128 × 40 × 12 × 12, ψ 4.8 GiB),
RTX 4070 Ti, per full time step:

| Stage | Time |
|---|---|
| T_Z + T_r (of which cuFFT 196 ms) | 245 ms |
| T_xy (of which cuFFT 137 ms) | 161 ms |
| Potential: Legendre FBR → DVR 75 ms, φ FFTs 2 × 48 ms, exp(iA) 60 ms, DVR → FBR 70 ms, rest 7 ms | 309 ms |
| **Total** | **716 ms** |

The whole program ran 80 steps in 54 s of propagation time, against ≈ 4.4 s per step on 8 CPU
threads: **≈ 6.5×**.

This GeForce card computes FP64 at 1/64 of its FP32 rate (≈ 0.63 TFLOP/s), and the FFTs and
Legendre sums of one step need ≈ 130 GFLOP. A step is therefore mostly FP64-bound here: a pure
memory pass over ψ takes 25.5 ms (400 GB/s). A data-centre GPU (A100/H100: full FP64 rate,
2–3 TB/s) changes the balance completely. Optimisations should be judged on such a GPU.

## Files
| File | Content |
|---|---|
| `src/cuda/cuda_backend.{hpp,cpp}` | `CudaBackend`, the factory, device info, GPU memory estimate |
| `src/cuda/device_kinetic_operator.{hpp,cpp}` | T_Z + T_r, T_xy and the specular grid (plans, tables, launch order) |
| `src/cuda/device_potential_operator.{hpp,cpp}` | the potential stage (tables, exp(iA), slab loop) |
| `src/cuda/*_kernels.cu`, `kernels.hpp` | the kernels and their launch functions |
| `src/cuda/cuda_fft.{hpp,cpp}` | cuFFT plans (64-bit, slabbed batches, shared work area) |
| `src/cuda/device_buffer.{hpp,cpp}`, `device_complex.hpp`, `cuda_error.{hpp,cpp}`, `launch.cuh` | device memory, complex arithmetic, layout, errors, launch helpers |
| `tests/cuda/` | cuFFT vs FFTW, CPU vs GPU per analysis step, determinism, memory estimate, restart, gates 3/4 |
| `tools/benchmark/cuda_operator_timing.cpp` | per-stage and per-kernel timing |
