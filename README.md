```text
>>=======================================================================================<<
||                                                                                       ||
||   _____                          _____ ______  _____           ______  _   _ ______   ||
||  /  __ \   _      _             /  ___|| ___ \|  _  |          |  _  \| | | || ___ \  ||
||  | /  \/ _| |_  _| |_   ______  \ `--. | |_/ /| | | |  ______  | | | || | | || |_/ /  ||
||  | |    |_   _||_   _| |______|  `--. \|  __/ | | | | |______| | | | || | | ||    /   ||
||  | \__/\  |_|    |_|            /\__/ /| |    \ \_/ /          | |/ / \ \_/ /| |\ \   ||
||   \____/                        \____/ \_|     \___/           |___/   \___/ \_| \_|  ||
||                                                                                       ||
>>=======================================================================================<<
```

# Cpp-SPO-DVR

A modern C++20 rewrite of the Fortran SPO-DVR code: six-dimensional quantum wave-packet dynamics
of a diatomic molecule scattering from a surface (coordinates Z, r, x, y, θ, φ). It uses
split-operator propagation, FFT grids in Z, r, x and y, and a Gauss-Legendre/FFT DVR in θ and φ.

Status:
- Implemented and checked against the Fortran code: set-up, the hard wall test potential, the
  gas-phase rovibrational states, the wave-packet propagation, the on-the-fly asymptotic
  analysis (S-matrix and state-to-state probabilities) and the fine analysis of the BK time series
  (onlyAnalysis / alsoAnalysis; two Fortran interpolation bugs fixed, see
  docs/fortran_fine_analysis_bugs.md).
- Potentials: the hard wall test potential and the two analytic H₂/Cu(111) LEPS potentials of the
  Fortran code (Dai–Light and Persson), checked against the Fortran potential on the grid; the
  Dai–Light propagation and analysis are checked against the Fortran code as well.
- Surface symmetry (`UseSymmetry`): the potential is computed and stored only for symmetry-unique
  (x, y) points, for the point groups C2v, C4v, C6v, C3v and C3, after a self-check that the
  potential really has that symmetry.
- Implemented: checkpoint/restart (native files and Fortran-compatible files).
- Not yet implemented: the flux analysis and the Fortran potentials that need data files not
  included in the Fortran code.

## Dependencies

- **Build tools**: GCC ≥ 13 (C++20, OpenMP), CMake ≥ 3.24, Ninja, pkg-config.
- **Math libraries**, one of:
  - **Intel oneMKL** 2024 or later (recommended): BLAS/LAPACK and the FFTs (through its FFTW3
    interface);
  - FFTW3 + LAPACK (fallback when oneMKL is not installed).
- nlohmann/json and GoogleTest are downloaded by CMake (internet access needed on the first
  configure).

On Debian/Ubuntu (also under WSL2):

**1. Build tools** (always needed):
```bash
sudo apt install build-essential cmake ninja-build pkg-config
```

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

## Building
```bash
git clone https://github.com/ewfsmeets/Cpp-SPO-DVR.git
cd Cpp-SPO-DVR
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure     # unit and regression tests (a few minutes)
```
The configure step reports which libraries it uses, e.g.
`SPO-DVR: BLAS/LAPACK from oneMKL 2026.2 (sequential), FFT library MKL`.

Useful CMake options:

| Option | Default | Meaning |
|---|---|---|
| `-DSPODVR_USE_MKL=OFF` | `ON` | ignore oneMKL; use the system LAPACK and FFTW |
| `-DSPODVR_FFT_LIBRARY=FFTW` | `AUTO` | use the real FFTW for the FFTs even when oneMKL is present (`AUTO` picks oneMKL when found) |
| `-DSPODVR_ARCH=native` | `x86-64-v3` | `-march` value (`x86-64-v3` = AVX2; `x86-64-v4` = AVX-512) |
| `-DSPODVR_PGO=GENERATE/USE` | `OFF` | profile-guided optimisation, about 15% faster (see docs/performance.md) |
| `-DSPODVR_FAST_MATH=ON` | `OFF` | experimental: `-ffast-math`; about 10% faster, changes results at the 1e-12 level (see Benchmark) |
| `-DSPODVR_BUILD_TESTS=OFF` | `ON` | skip the tests |

Performance notes, reproducibility across machines and the PGO workflow are in
[docs/performance.md](docs/performance.md).

### Formatting and linting
clang-format runs as a git pre-commit hook; clang-tidy is run on demand. Both versions are pinned
in `.pre-commit-config.yaml` and installed by pre-commit itself (no system packages needed):
```bash
pipx install pre-commit                                      # once per machine
pre-commit install                                           # once per clone
pre-commit run --all-files                                   # format everything
pre-commit run clang-tidy --hook-stage manual --all-files    # lint (needs build/compile_commands.json)
```

## Running
The program reads one JSON input file:
```bash
./build/spodvr inputs/hardwall_small/input.json --threads 8 --output-dir runs/hardwall_small
```

Command-line options:
```
spodvr <input.json> [--threads N] [--schedule dynamic|static] [--output-dir DIR] [--no-color] [--help] [--version]
```
- `--threads N`: number of OpenMP threads (default: all). The results do not depend on it. Use the
  number of physical cores: hyperthreads do not speed up these memory- and FFT-bound loops.
- `--schedule dynamic|static`: how the (Z, r) points are shared out over the threads. `dynamic`
  balances the uneven work best (one NUMA node, i.e. one CPU socket); `static` keeps every thread on
  the memory it first wrote (several NUMA nodes). Default: chosen from the number of NUMA nodes. The
  results do not depend on it.
- `--output-dir DIR`: where the data files are written (default: the current directory).

Output:
- **Screen:** an input summary, the grids, a memory estimate, the surface symmetry (with
  `UseSymmetry`), the gas-phase states, the initial wave packet, the analysis basis, a table of averages during the propagation, and compact probability tables
  (one row per energy) at the end.
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

- The JSON keys mirror the Fortran `DVR.inp` field names. All keys, defaults, units and checks
  are described in [docs/input_schema.md](docs/input_schema.md).
- Unknown keys and inconsistent parameters are reported as errors.
- The test copies of these inputs (with the Fortran reference output) live in
  `tests/reference/`, together with the corrugated test cases (`cu111_dailight_*`,
  `cu111_persson_small`), whose `input.json` files can be run in the same way.
- The example grids are deliberately tiny so that they run quickly. Production calculations use
  much larger grids (more than 64 GB of memory).

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
slower in wall time; this is still being investigated. The fast-math build passes all regression
gates, but its results differ from the default build at the 1e-12 level, so it is not the default.
Details in [docs/performance.md](docs/performance.md).

## Repository layout
| Path | Content |
|---|---|
| `src/` | the program: `input/`, `grids/`, `potential/`, `symmetry/`, `gasphase/`, `propagation/`, `analysis/`, `io/`, `fft/`, `linalg/`, `core/`, `report/`, `app/` |
| `tests/` | unit tests (one file per module) and regression tests against Fortran reference output (`tests/reference/`) |
| `tools/` | Fortran reference scripts, benchmark and comparison tools |
| `docs/` | input schema, output files, performance notes |
| `inputs/` | example inputs |

## License
This program is free software, licensed under the GNU General Public License, version 3. See
[LICENSE](LICENSE).
