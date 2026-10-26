# Cpp-SPO-DVR

A modern C++20 rewrite of the Fortran SPO-DVR code: six-dimensional quantum wave-packet dynamics
of a diatomic molecule scattering from a surface (coordinates Z, r, x, y, θ, φ). It uses
split-operator propagation, FFT grids in Z, r, x and y, and a Gauss-Legendre/FFT DVR in θ and φ.

Status:
- Implemented: set-up, the hard wall test potential, the gas-phase rovibrational states and the
  wave-packet propagation. All of them are checked against the Fortran code.
- In progress: the asymptotic analysis (S-matrix, state-to-state probabilities) and
  checkpoint/restart.

## Dependencies

| Dependency | Version | Notes |
|---|---|---|
| C++ compiler | GCC ≥ 13 (C++20, OpenMP) | |
| CMake | ≥ 3.24 | |
| Ninja | any | the recommended generator |
| pkg-config | any | |
| **Intel oneMKL** | 2024 or later | **recommended**: BLAS/LAPACK and the FFTs (through its FFTW3 interface) |
| FFTW3 + LAPACK | any | the fallback when oneMKL is not installed |
| nlohmann/json, GoogleTest | — | downloaded automatically by CMake (internet access needed on the first configure) |

### Intel oneMKL (recommended)
On Debian/Ubuntu (also under WSL2), install oneMKL from Intel's apt repository:
```bash
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-mkl-devel
```
- CMake finds oneMKL in `$MKLROOT` or in `/opt/intel/oneapi/mkl/latest`; no `setvars.sh` is
  needed for building.
- The program always uses the sequential, LP64 oneMKL.

### Without oneMKL
```bash
sudo apt install build-essential cmake ninja-build pkg-config libfftw3-dev liblapack-dev
```
The build tools themselves (`build-essential cmake ninja-build pkg-config`) are also needed with
oneMKL.

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
| `-DSPODVR_BUILD_TESTS=OFF` | `ON` | skip the tests |

Performance notes, reproducibility across machines and the PGO workflow are in
[docs/performance.md](docs/performance.md).

## Running
The program reads one JSON input file:
```bash
./build/spodvr inputs/hardwall_small/input.json --threads 16 --output-dir runs/hardwall_small
```

Command-line options:
```
spodvr <input.json> [--threads N] [--output-dir DIR] [--no-color] [--help] [--version]
```
- `--threads N`: number of OpenMP threads (default: all). The results do not depend on it.
- `--output-dir DIR`: where the data files are written (default: the current directory).

Output:
- **Screen:** an input summary, the grids, the asymptotic potential, the gas-phase states, the
  initial wave packet, and a table of averages during the propagation.
- **Data files** (in the output directory, full precision):
  - `asymptotic_potential.dat`;
  - `gas_phase_states.dat`;
  - `averages.dat` (time, norm, ⟨Z⟩, ⟨r⟩, rest norm).

A run of the provided example (7000 time steps) takes about 3 minutes with 16 threads on a desktop machine (AMD Ryzen 7 5800X).

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
  `tests/reference/`.
- The example grids are deliberately tiny so that they run quickly. Production calculations use
  much larger grids (more than 64 GB of memory).

## Repository layout
| Path | Content |
|---|---|
| `src/` | the program: `input/`, `grids/`, `potential/`, `gasphase/`, `propagation/`, `analysis/`, `fft/`, `linalg/`, `report/`, `app/` |
| `tests/` | unit tests (one file per module) and regression tests against Fortran reference output (`tests/reference/`) |
| `tools/` | Fortran reference scripts, benchmark and comparison tools |
| `docs/` | input schema, performance notes |
| `inputs/` | example inputs |
