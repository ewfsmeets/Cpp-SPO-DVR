```text
>>=============================================================================================<<
 _____                  _____  _   _ ______   ___     _____ ______  _____   ______  _   _ ______
/  __ \   _      _     /  __ \| | | ||  _  \ / _ \   /  ___|| ___ \|  _  |  |  _  \| | | || ___ \
| /  \/ _| |_  _| |_   | /  \/| | | || | | |/ /_\ \  \ `--. | |_/ /| | | |  | | | || | | || |_/ /
| |    |_   _||_   _|  | |    | | | || | | ||  _  |   `--. \|  __/ | | | |  | | | || | | ||    /
| \__/\  |_|    |_|    | \__/\| |_| || |/ / | | | |  /\__/ /| |    \ \_/ /  | |/ / \ \_/ /| |\ \
 \____/                 \____/ \___/ |___/  \_| |_/  \____/ \_|     \___/   |___/   \___/ \_| \_|
>>=============================================================================================<<
```

# Cpp-SPO-DVR

A modern C++20 rewrite of the Fortran SPO-DVR code: six-dimensional quantum wave-packet dynamics
of a diatomic molecule scattering from a surface (coordinates Z, r, x, y, θ, φ). It uses
split-operator propagation, FFT grids in Z, r, x and y, and a Gauss-Legendre/FFT DVR in θ and φ.

## About this project

This project is a full reimplementation of the original SPO-DVR Fortran 90 code, written by
Dr. Mark F. Somers during his PhD at Leiden University:

- M. F. Somers, *Six-Dimensional Quantum Dynamics of Dissociative Chemisorption of Hydrogen on
  Metal Surfaces*, PhD thesis, Leiden University (2004),
  [repository](https://scholarlypublications.universiteitleiden.nl/handle/1887/3250648).
- G. J. Kroes and M. F. Somers, *Six-dimensional dynamics of dissociative chemisorption of H₂ on
  metal surfaces*, J. Theor. Comput. Chem. 4, 493 (2005),
  [DOI](https://www.worldscientific.com/doi/abs/10.1142/S0219633605001647).

I am very grateful that I could use the original Fortran code during my own PhD. Making it faster
and porting it to CUDA has been a years-long endeavour: faster calculations open the door to more
avenues of approach with quantum dynamics calculations.

Over the years SPO-DVR has been my go-to code whenever I wanted to apply newly learned programming
skills to highly complex, performance-critical code. I tinkered with and refactored small parts of
it by hand, but it took Claude to pull all the ideas together into a working program and to
benchmark its accuracy. Without Claude this would have been too big for a hobby project.

## Status

- Implemented and checked against the Fortran code: set-up, the hard wall test potential, the
  gas-phase rovibrational states, the wave-packet propagation, the on-the-fly asymptotic
  analysis (S-matrix and state-to-state probabilities) and the fine analysis of the BK time series
  (onlyAnalysis / alsoAnalysis; two Fortran interpolation bugs fixed, see
  docs/fortran_fine_analysis_bugs.md).
- Potentials: the hard wall test potential and the two analytic H₂/Cu(111) LEPS potentials of the
  Fortran code (Dai–Light and Persson), checked against the Fortran potential on the grid; the
  Dai–Light propagation and analysis are checked against the Fortran code as well.
- CRP potentials from a data set (SystemID `CRPlibrary`, input key `CRPlibDataSet`): the Fortran
  CRP library of M. Wijzenbroek (C3v, C4v and C6v surfaces), ported as the self-contained library
  `src/crplib/` (V and analytic gradients). The library is checked against the Fortran library on
  synthetic data sets (bitwise) and on the real H₂/Cu(111) SRP48 data set in `potentials/`, and the
  program against the Fortran potential on the grid for both.
- Surface symmetry (`UseSymmetry`): the potential is computed and stored only for symmetry-unique
  (x, y) points, for the point groups C2v, C4v, C6v, C3v and C3, after a self-check that the
  potential really has that symmetry.
- Implemented: checkpoint/restart (native files and Fortran-compatible files).
- Optional GPU backend (CUDA, `--device cuda`): the wave packet stays in GPU memory; passes the same
  Fortran regression gates as the CPU code, and agrees with it to ≈1e-13 (see
  [docs/cuda.md](docs/cuda.md)).
- Not yet implemented: the flux analysis. A propagating run with `doFluxAnalysis: true` is refused
  with a "not implemented yet" error ([docs/input_schema.md](docs/input_schema.md)).
- Not yet implemented: the other Fortran potentials (their data files are in `potentials/`).
- Fortran version: the port follows SPO-DVR 5.19 (2023), and every comparison uses it. The changes
  of the later upstream versions (up to 5.24, 2025) are not ported yet; see
  [docs/fortran_deviations.md](docs/fortran_deviations.md).

## Dependencies

- **Build tools**: a C++20 compiler with OpenMP (tested with GCC 15, the Ubuntu 26.04 default), CMake ≥ 3.24,
  Ninja, pkg-config (CMake looks up FFTW3 through it).
- **Math libraries**, one of:
  - **Intel oneMKL** (recommended; tested with 2026.2, found through its `MKLConfig.cmake`):
    BLAS/LAPACK and the FFTs (through its FFTW3 interface);
  - FFTW3 + LAPACK (fallback when oneMKL is not installed).
- nlohmann/json and GoogleTest are downloaded by CMake (internet access needed on the first
  configure).
- **Optional, for the GPU backend**: an NVIDIA GPU with its driver, and the CUDA toolkit (nvcc,
  cudart, cuFFT, and at run time nvJitLink and NVRTC for the cuFFT callback). Recommended: the newest
  CUDA 13 the GPU driver supports (tested with 13.3,
  the first release that accepts GCC 15; nvcc uses the same compiler as the rest of the build). FP64 speed
  matters: data-centre GPUs (A100, H100) run FP64 at 1/2 of the FP32 rate, GeForce cards at 1/64.

On Debian/Ubuntu (also under WSL2):

**1. Build tools** (always needed):
```bash
sudo apt install build-essential cmake ninja-build pkg-config
```
Any C++20 GCC works; the project is tested with GCC 15 (Ubuntu 26.04's default). For another
compiler, configure a new build directory with e.g. `CXX=g++-14 cmake ...`; nvcc then uses it too.

**2a. Intel oneMKL** (recommended), from Intel's apt repository:
```bash
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg]" \
     "https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-mkl-devel
```
CMake finds oneMKL in `$MKLROOT` or in `/opt/intel/oneapi/mkl/latest` (no `setvars.sh` needed for
building). The program always uses the sequential, LP64 oneMKL.

**2b. Without oneMKL**, FFTW3 and LAPACK instead:
```bash
sudo apt install libfftw3-dev liblapack-dev
```

**3. CUDA** (optional, for `--device cuda`), CUDA 13 from NVIDIA's apt repository (repository
`wsl-ubuntu` under WSL2, `ubuntu2604` / `ubuntu2404` on a native Ubuntu 26.04 / 24.04):
```bash
nvidia-smi                                  # the GPU must be visible; note "CUDA Version" (newest CUDA the driver runs)
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install cuda-toolkit-13-3          # newest 13.x not above the driver's CUDA version; installs into /usr/local/cuda
```
- Do not install a toolkit newer than the driver supports: it compiles, but the kernels then fail
  at run time.
- **WSL2:** the GPU driver is the Windows NVIDIA driver (update it in Windows). Never install a
  Linux driver inside WSL; the `wsl-ubuntu` repository contains none.
- `/usr/local/cuda/bin` need not be on `PATH`: CMake finds nvcc there. It also sets the GPU
  architecture from `nvidia-smi` (override with `-DCMAKE_CUDA_ARCHITECTURES`, e.g.
  `"80-real;90-real"` for A100 and H100).
- Ubuntu's own packages (`sudo apt install nvidia-cuda-toolkit`) work as well, but are older
  (Ubuntu 24.04: CUDA 12.0 with g++-12 as nvcc's host compiler; untested since the kernels became C++20).
- nvcc compiles with the C++ compiler of the build. If the CUDA release rejects that GCC
  ("unsupported GNU version"; GCC 15 needs CUDA ≥ 13.3), set `CUDAHOSTCXX` or
  `-DCMAKE_CUDA_HOST_COMPILER` to a GCC it accepts.

## Building
```bash
git clone https://github.com/ewfsmeets/Cpp-SPO-DVR.git
cd Cpp-SPO-DVR
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```
The configure step reports which libraries it uses, e.g.
`SPO-DVR: BLAS/LAPACK from oneMKL 2026.2 (sequential), FFT library MKL`.
How to run the tests: see [Testing](#testing).

With the GPU backend:
```bash
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_USE_CUDA=ON
cmake --build build-cuda
```
The kernels are compiled for the GPU of the build machine (`native`); for another GPU set
`-DCMAKE_CUDA_ARCHITECTURES` (`89-real` RTX 40xx, `80-real` A100, `90-real` H100). The configure step prints the nvcc version, host compiler and architectures.

Useful CMake options:

| Option | Default | Meaning |
|---|---|---|
| `-DSPODVR_USE_MKL=OFF` | `ON` | ignore oneMKL; use the system LAPACK and FFTW |
| `-DSPODVR_FFT_LIBRARY=FFTW` | `AUTO` | use the real FFTW for the FFTs even when oneMKL is present (`AUTO` picks oneMKL when found) |
| `-DSPODVR_ARCH=native` | `x86-64-v3` | `-march` value (`x86-64-v3` = AVX2; `x86-64-v4` = AVX-512) |
| `-DSPODVR_PGO=GENERATE/USE` | `OFF` | profile-guided optimisation, two builds and a training run (see docs/performance.md) |
| `-DSPODVR_FAST_MATH=ON` | `OFF` | experimental: `-ffast-math` (with CUDA also `--use_fast_math`, i.e. FMA, in the GPU kernels); changes results at the 1e-12 to 1e-11 level (see Benchmark) |
| `-DSPODVR_USE_CUDA=ON` | `OFF` | build the GPU backend (`--device cuda`, NVIDIA GPU, CUDA toolkit; see [docs/cuda.md](docs/cuda.md)) |
| `-DSPODVR_BUILD_TESTS=OFF` | `ON` | skip the tests |

Performance notes, reproducibility across machines and the PGO workflow are in
[docs/performance.md](docs/performance.md).

### Build variants and `tools/build_all.sh`
The development builds live side by side in git-ignored directories.
`tools/build_all.sh [variant ...]` configures and builds them (all by default; the CUDA ones only
when nvcc is found). It does not use PGO; for that, see docs/performance.md:

| Directory | Options | Purpose |
|---|---|---|
| `build` | (defaults) | the reference build: oneMKL, no FMA, LTO |
| `build-fftw` | `SPODVR_FFT_LIBRARY=FFTW` | oneMKL LAPACK with the real FFTW |
| `build-nomkl` | `SPODVR_USE_MKL=OFF` | system LAPACK + FFTW: the fallback that must keep working |
| `build-poison` | `SPODVR_POISON_UNINITIALISED=ON` | NaN-fills new large arrays: finds reads before the first write |
| `build-fast` | `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fast CPU build (results differ at the 1e-12 to 1e-11 level; add PGO by hand for the fastest) |
| `build-cuda` | `SPODVR_USE_CUDA=ON` | GPU backend, no FMA (results comparable with the CPU build) |
| `build-cuda-fast` | `SPODVR_USE_CUDA=ON`, `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fastest GPU build (FMA in the kernels; no PGO) |

Native restart files (`checkpoint.bin`, `bk_coefficients.bin`) only work between builds with the
same `SPODVR_FAST_MATH` setting. A restart or `onlyAnalysis` run checks the stored set-up parameters
exactly, and the fast-math builds compute one of them (the effective `ZOptFactor`) 1 ulp
differently, so a default build refuses a fast-math build's native files and vice versa. Files in
the Fortran format (`CheckpointFormat` `"fortran"` or `"both"`, read with `RestartFormat`
`"fortran"`) carry no parameter list and can be used across builds.

### Formatting and linting
clang-format runs as a git pre-commit hook; clang-tidy is run on demand. Both versions are pinned
in `.pre-commit-config.yaml` and installed by pre-commit itself (no system packages needed):
```bash
pipx install pre-commit                                      # once per machine
pre-commit install                                           # once per clone
pre-commit run --all-files                                   # format everything
pre-commit run clang-tidy --hook-stage manual --all-files    # lint (needs build/compile_commands.json)
```

## Testing
```bash
OMP_NUM_THREADS=8 ctest --test-dir build --output-on-failure                # everything (~5 min)
OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R Cuda --output-on-failure   # GPU tests (~1 min)
```
- **Run ctest serially** (no `-j`): many tests already use all OpenMP threads. Set
  `OMP_NUM_THREADS` to the number of physical cores (not hyperthreads).
- **What is tested:**
  - unit tests, one program per module (`tests/test_<module>.cpp`);
  - regression tests against the original Fortran code (`tests/test_regression_*.cpp`, test names
    `Regression…`/`Gate…`): potential, gas-phase states, propagation (norms and averages),
    on-the-fly analysis (probabilities, BK coefficients) and fine analysis. The Fortran reference
    output is stored in `tests/reference/`;
  - GPU tests (`tests/cuda/`, names `Cuda…`): only in builds with `-DSPODVR_USE_CUDA=ON`. They skip on
    a machine without a GPU.
- **Selecting tests** (ctest `-R`/`-E` take a regular expression on the test name;
  `ctest --test-dir build -N` lists them all):
  ```bash
  ctest --test-dir build -E Regression           # unit tests only (~1.5 min)
  ctest --test-dir build -R Gate3                # one acceptance gate
  ./build/tests/test_input --gtest_filter='Input.*'   # one test program directly (GoogleTest)
  ```
- **Long tests:** the full 7000-step regression runs (hard wall and Dai–Light; the latter alone takes
  about 13 minutes on 8 cores) are skipped unless `SPODVR_LONG_TESTS` is set:
  ```bash
  SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build -R Long --output-on-failure
  SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R 'Cuda.*Long' --output-on-failure
  ```
- **Skipped regression tests:** a regression test skips (it does not fail) when its Fortran
  reference file is missing. The scripts in `tools/fortran_reference/` regenerate the files; they
  need a Fortran reference install ([tools/fortran_reference/README.md](tools/fortran_reference/README.md)).
  `ctest -R ReferenceManifest` fails when a file listed in a case's `MANIFEST.tsv` is missing or
  edited, or a case has no `MANIFEST.tsv` ([MANIFEST.tsv](tools/fortran_reference/README.md#manifesttsv)).
- **Other build variants** ([Build variants](#build-variants-and-toolsbuild_allsh)): run the same
  ctest command in their build directory. Use `build-poison` after touching the uninitialised-memory
  code. The fast-math builds (`build-fast`, `build-cuda-fast`) are expected to fail a few unit tests
  by one ulp, and with GCC 15 also the hard-wall V<sub>asy</sub> table of gate 1 (1.6e-14 against 1e-14);
  their other gates pass.
- **Docs consistency:** `python3 -I -B tools/check_docs.py [--spodvr build/spodvr]` (under 1 s)
  checks the docs against the code: input keys, output files, links and anchors, decision IDs,
  local paths, command-line options (with `--spodvr`), gates and environment variables; it also
  lists numbers repeated across docs (report-only). ctest runs it as `Docs.Consistency`
  (`ctest --test-dir build -R Docs`). Accepted findings, each with a reason, are in
  `tools/check_docs_allowlist.tsv`; a row that no longer suppresses anything fails the check.
- **Numerically neutral changes** (refactors) must leave the results bitwise unchanged:
  `tools/check_bitwise_refactor.sh [base-ref]` (default `HEAD`). It builds both trees in a scratch
  directory; to configure without network access, set `FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON` to
  the nlohmann_json sources. For changes that may alter rounding:
  `tools/compare_propagation.sh [base-ref] [tolerance]`. The final wave packet of two builds, or of
  the CPU and the GPU: `spodvr_dump_wave` and `spodvr_compare_wave` ([docs/performance.md](docs/performance.md)).
- **Production-grid validation:** beyond the gates, the full thesis calculation
  (`inputs/thesis_fig5_3_dz`, 7000 steps) was run with the Fortran, the C++ CPU and the C++ GPU
  code and compared state by state: [docs/validation_thesis_fig5_3_dz](docs/validation_thesis_fig5_3_dz/README.md).
  With the default (strict) GPU build, all 111 reaction probabilities match the Fortran's 7 printed digits.

## Running
The program reads one JSON input file:
```bash
./build/spodvr inputs/hardwall_small/input.json --threads 8 --output-dir runs/hardwall_small
```

Command-line options:
```
spodvr <input.json> [--threads N] [--schedule dynamic|static] [--device cpu|cuda] [--output-dir DIR] [--no-color]
       [--help] [--version]
```
- `--device cpu|cuda`: propagate on the CPU (default) or on an NVIDIA GPU, where the wave packet
  stays in GPU memory (needs a build with `-DSPODVR_USE_CUDA=ON`; see [docs/cuda.md](docs/cuda.md)).
  The results agree with the CPU's to rounding (≈1e-13, the level of a different FFT library).
- `--threads N`: number of OpenMP threads (default: all). The results do not depend on it. Use the
  number of physical cores: hyperthreads do not speed up these memory- and FFT-bound loops.
- `--schedule dynamic|static`: how the (Z, r) points are shared out over the threads. `dynamic`
  balances the uneven work best (one NUMA node, i.e. one CPU socket); `static` keeps every thread on
  the memory it first wrote (several NUMA nodes). Default: chosen from the number of NUMA nodes. The
  results do not depend on it.
- `--output-dir DIR`: where the data files are written (default: the current directory).

Output:
- **Screen:** an input summary, the grids, a memory estimate, the surface symmetry (with
  `UseSymmetry`), the gas-phase states, the initial wave packet, the analysis basis, a table of averages during the propagation (on a terminal with a live progress line:
  step, time per step, estimated remaining time), and compact probability tables (one row per
  energy) at the end. Colour and the progress line are off when the output is not a terminal (a
  pipe, `tee`, a file) or with `--no-color`.
- **Data files** (in the output directory, full precision, tab-separated with a commented header
  that explains every column; they load with `pandas.read_csv(f, sep="\t", comment="#")`):

| File | Content |
|---|---|
| `scattering_summary.tsv` | per energy: total scattering and reaction probability, P(v), average final energies |
| `s_matrix.tsv` | per energy and open final state (v, j, mJ, n, m): S-matrix element, probability, phase |
| `probabilities_rovib.tsv` | per energy and (v, j, mJ), summed over diffraction |
| `probabilities_diffraction.tsv` | per energy, (v, j, mJ) and diffraction order |
| `analysis.json` | run metadata: basis sizes, initial state, energy grid, file list |
| `averages.dat` | time, norm, ⟨Z⟩, ⟨r⟩, rest norm during the propagation |
| `checkpoint.bin`, `bk_coefficients.bin` | final wave packet and Balint-Kurti coefficients, for a restart |
| `asymptotic_potential.dat`, `gas_phase_states.dat` | set-up data |

With `alsoAnalysis` (or an analysis-only run, `onlyAnalysis`) the fine-analysis results are written
to `fine_analysis/` in the output directory. The file formats are described in
[docs/output_files.md](docs/output_files.md).

### Continuing a run (restart)
Set `"NewCalc": false`, `"tStart"` to the end time of the earlier run, and `"RestartDir"` to its
output directory, and use a *different* `--output-dir`. All other parameters must stay the same
(the program checks this). The S-matrix of the restarted run covers the whole propagation; a
split run gives bitwise the same results as one long run. Runs can also be continued from the
Fortran code's files (`"RestartFormat": "fortran"`: fort.8/fort.9 and fort.22/fort.21) and
written for it (`"CheckpointFormat": "fortran"` or `"both"`); see
[docs/input_schema.md](docs/input_schema.md).

A run of the provided example (7000 time steps) takes about 3.5 minutes with 8 threads on a desktop
machine (AMD Ryzen 7 5800X).

### Running a oneMKL build from another location
Inside the build tree the oneMKL library path is built into the executable. If you copy or install
the executable elsewhere, make the oneMKL libraries findable first:
```bash
source /opt/intel/oneapi/setvars.sh        # or: export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib
```
Optional: `export MKL_CBWR=AVX2` gives bitwise identical results across machines with different
vector extensions (see docs/performance.md).

### Inputs
`inputs/` holds ready-to-run example inputs:

| Directory | Content |
|---|---|
| `inputs/hardwall_small/` | The hard wall reference case (H₂ on a hard wall, Cu(211)-like grid, 0.2–0.6 eV). `input.json` is the C++ input. `DVR.inp` is the same calculation for the Fortran reference code. |
| `inputs/thesis_fig5_3_dz/`, `inputs/thesis_fig5_3_er/` | H₂ + Cu(111) on the Dai–Zhang and Persson (ER) PES, the two curves of Fig. 5.3 of M. F. Somers' PhD thesis, at the grid sizes of the thesis (except a finer Z grid, see their READMEs; JMax 28, ψ ≈ 6 GB; hours per run). Each directory has a README with the parameter choices. The DZ case is the GPU benchmark. |

- The JSON keys mirror the Fortran `DVR.inp` field names. All keys, defaults, units and checks
  are described in [docs/input_schema.md](docs/input_schema.md).
- How to choose the Z grids, the wave packet, the analysis lines and the absorbers for a new
  input: [inputs/README.md](inputs/README.md).
- Unknown keys and inconsistent parameters are reported as errors.
- The test copies of these inputs (with the Fortran reference output) live in
  `tests/reference/`, together with the corrugated test cases (`cu111_dailight_*`,
  `cu111_persson_small`, and `crp_*_synthetic_small`, which must be run from the repository root
  because their `CRPlibDataSet` path is relative to it), whose `input.json` files can be run in
  the same way.
- The hard wall and test grids are deliberately tiny so that they run quickly. Production
  calculations use much larger grids (more than 64 GB of memory).

## Benchmark
Full Dai–Light run (`tests/reference/cu111_dailight_small`): 7000 time steps on a 176 × 56 × 8 × 8
(Z, r, x, y) grid with 15 rotational states, plus the on-the-fly analysis. AMD Ryzen 7 5800X
(8 cores), 8 threads, the three runs one after another:

| Binary | Wall time | CPU time | Threads busy on average |
|---|---|---|---|
| Fortran | 680 s | 5255 s | 7.7 |
| C++, current flags | 764 s | 5040 s | 6.6 |
| C++, fast-math + native + PGO | 675 s | 4489 s | 6.7 |

The C++ code needs less CPU time than the Fortran code, but keeps the threads less busy, so it is
slower in wall time; this is still being investigated. The fast-math build passes the regression
gates except the hard-wall V<sub>asy</sub> table (by 1.6e-14, with GCC 15), and its results differ from the
default build: by ~1e-12 on the gate grids, by 7.8e-12 on the thesis grid (GPU builds)
([validation report](docs/validation_thesis_fig5_3_dz/README.md#strict-build-gpu-run-added-2026-11-08)).
It is therefore not the default.
Details in [docs/performance.md](docs/performance.md).

GPU (`--device cuda`), DZ case of `inputs/thesis_fig5_3_dz` as of 3c1d91b (128 × 40 × 12 × 12 grid, 435 rotational
states, 4.8 GiB wave function), NVIDIA RTX 4070 Ti (a GeForce card: FP64 at 1/64 rate) against
8 threads of the Ryzen 7 5800X, per time step:

| Build | Time per step | Speed-up |
|---|---|---|
| CPU, default flags | 3.6 s | 1× |
| GPU, `build-cuda` (no FMA; agrees with the CPU to ≈1e-14 here) | 0.61 s | 5.9× |
| GPU, `build-cuda-fast` (FMA; S matrix of the thesis run 7.8e-12 from `build-cuda`) | 0.54–0.57 s | 6.3–6.7× |

The 10 000 steps of the benchmark take about 1.5–1.8 h on the GPU instead of about 10 h. Details
and the limits of this card in [docs/cuda.md](docs/cuda.md).

## Repository layout
| Path | Content |
|---|---|
| `src/` | the program: `input/`, `grids/`, `potential/`, `symmetry/`, `gasphase/`, `propagation/`, `analysis/`, `io/`, `fft/`, `linalg/`, `core/`, `report/`, `app/`, `cuda/` (optional GPU backend) |
| `tests/` | unit tests (one file per module) and regression tests against Fortran reference output (`tests/reference/`) |
| `tools/` | Fortran reference scripts, benchmark and comparison tools, the docs consistency checker |
| `docs/` | input schema, output files, validation against the Fortran (`validation.md`, `fortran_deviations.md`), performance notes, GPU backend |
| `inputs/` | example inputs |
| `potentials/` | data files of the Fortran potentials, one directory per system (verbatim copy of `potential.dir` of the upstream Fortran repository) |

## License
This program is free software, licensed under the GNU General Public License, version 3. See
[LICENSE](LICENSE).
