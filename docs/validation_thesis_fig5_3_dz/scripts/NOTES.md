# Reference notes

Copied from the task prompt (checked while preparing it). Additions found during the work are
under "Additions" at the end.

- **Logs:** `script` writes CRLF. In the C++ logs the averages rows and the "analysis surface …
  fixed" line sit inside `\r` progress redraws: split on `\r` first, then strip
  `\x1b\[[0-9;?]*[A-Za-z]`, then drop progress-bar fragments (`▕`). Take averages from
  `averages.dat` / `averages.out`. The Fortran stdout is block-buffered (complete only at exit).
- **Joins:** C++ prints 0-based indices, Fortran 1-based (counts are equal). Join on quantum
  numbers and energy index, never on row order. The Fortran per-energy "total energy" is
  absolute, like C++ `E_total`.
- **Fortran numbers:** `D`→`E`, 3-digit exponents drop the letter, `****` = overflow; fields are
  separated by at least one blank.
- **On-the-fly diffraction-order probabilities** are on Fortran stdout; unit 10/fil10 is written
  only by the fine analysis, which did not run (the input echo "file 10" is a legacy label).
- **fort.25** (writeProbs): per energy a header, then per (v, j, mJ, order) rows `n m ReS ImS
  |S|² phase` (d13.6/f10.6), v printed as iv−1. Rows exist only for open rovibrational channels,
  including their closed diffraction channels (S = 0), which `s_matrix.tsv` omits. Parser:
  `tests/fine_analysis_reference.hpp` `read_fortran_fort25`. The on-the-fly fort.25 has never been
  compared in a gate (gate 5 uses the fine-analysis one).
- **Full-precision Fortran data:** `.pes` (V in hartree before the −dt factor, symmetry-reduced;
  record layout as in `tools/fortran_reference/extract_pes.py`, whose loader needs ~8 GB here:
  memory-map instead; symmetryIndex is a header record, then one record per (Z, r) block);
  `fort.21` (BK coefficients of every analysis step; reader `src/io/fortran_files.cpp`
  `read_fortran_bk_file` (checks sizes, dt, labels; time = istep·dtAna); same element order as
  `bk_coefficients.bin`); `fort.9` (final ψ, one record per (Z, r) block, same (irot, ixy, izr)
  order as `checkpoint.bin`; record 3 = final analysis-surface state); the V_asy table on stdout.
- **C++ 6D potential** is not written (`"Potential": true` writes only the gas-phase files):
  compute it with a helper following `tests/test_regression_pes.cpp` (build_config →
  make_spatial_grids → make_angular_grid → make_potential → compute_potential_grid; seconds).
  V and exp(−iV dt) of the GPU run are computed on the host (`app/run.cpp`, `cuda/cuda_backend.cpp`).
- **Helper build recipe:** `lib*.a` hold slim LTO objects (machine code is generated at the
  helper's link, so it shares source and flags with the run, not machine code). Compile with g++
  15 using the flags of `build-fast/compile_commands.json` (incl. the `-D` defines and the oneMKL
  `include/fftw` path; nlohmann headers in `build-fast/_deps/`), link with the line in
  `build-fast/build.ninja` (`-O3 -DNDEBUG -flto=auto -fno-fat-lto-objects`, oneMKL sequential,
  libgomp). Do not `add_subdirectory` the repo (FetchContent clones from GitHub). Same for
  `build-cuda-fast`.
- **Replay recipe** (precedent: the restart path re-accumulates through
  `add_earlier_analysis_step`, `app/run.cpp`; `tests/test_restart.cpp` checks bitwise
  `s_matrix.tsv`): `build_config(read_input_file)`, potential/angular/spatial grids,
  `compute_asymptotic_potential`, `compute_gas_phase` (FftwFft, LapackEigensolver),
  `RotationalBasis`, `AnalysisBasis(…, parallel_kinetic_energies(config))`,
  `SMatrixAccumulator(config, basis, gas)`, then `add(bk, time)` per record from `read_bk_file` /
  `read_fortran_bk_file`, `results()`, `write_scattering_files`. The specular replacement after
  the analysis line is fixed is already in the stored coefficients (both codes).
- **`onlyAnalysis` commands** (from the repo root; no jq installed):
  `python3 -I -c 'import json,sys;j=json.load(open(sys.argv[1]));j.update(onlyAnalysis=True,RestartDir=sys.argv[2],RestartFormat=sys.argv[3]);json.dump(j,open(sys.argv[4],"w"),indent=1)' inputs/thesis_fig5_3_dz/input.json runs/thesis_fig5_3_dz_fortran fortran <scratch>/oa_fortran.json`,
  then `OMP_NUM_THREADS=8 build-fast/spodvr <scratch>/oa_fortran.json --threads 8 --output-dir runs/thesis_fig5_3_dz_oa_fortran`.
  C++ runs: their run directory with `native`. `fortran` reads fort.22, else fort.21; `native`
  reads `bk_coefficients.bin`; neither needs ψ. Output: `fine_analysis/*.tsv` (interpolated to
  the analysis energies; `spectrum/` has no `s_matrix.tsv` with writeSMatrix false). < 1 min each.
- **Analysis line (DynamZSpAna):** the C++ log gives Z and t of the fixing, the Fortran only Z.
  Pin the Fortran step through `fort.21`: from then on the initial-channel coefficient is the
  specular projection (Analyzor.f90 BalintKurtiCoeffs; docs/fortran_deviations.md). The margin of
  that decision is not in any output.
- **Timings:** the `/usr/bin/time` block in each log; C++ "Propagation time"; Fortran "CPU TIMES"
  block (wall seconds, omp_get_wtime) and `fort.99`. GPU memory: only the C++ estimate (not
  measured).
- **No counterpart** in the other program (say so, don't invent one): e.g. Fortran nmoder,
  nstatw, dφ; C++ θ weights, symmetry self-check, memory estimate. FGH diagonalisation: Fortran
  Jacobi, C++ LAPACK.
- **`.gitignore`** ignores `*.dat`, `fort.*`, `fil10*`, `*.o`, `*.a` and `runs/` (not `*.log`,
  not `__pycache__`): rename copies in the docs folder (e.g. `averages_cpu.txt`). pre-commit is
  installed: trailing-whitespace, end-of-file-fixer, check-json and check-yaml run on every
  commit (clang-format only on C++).
- **Thesis PDF:** 217 pages, vector figures. The Read tool shows pages (`pages` parameter; find
  Fig. 5.3 via the list of figures). Exact curve extraction: pymupdf `get_drawings` (D1).

## Additions (found during the work)

### Inventory: per item, file and precision in each program
F = Fortran (cleaned stdout `logs/fortran.txt` unless noted), C = C++ (CPU and GPU identical
layout; cleaned logs `logs/cpu.txt`, `logs/gpu.txt`; files in the run directories).

| Item | Fortran (precision) | C++ (precision) |
|---|---|---|
| 1 input | `DVR.inp` (text); "THE INPUT READ" echo (E15.8, I, L); "SETUP INTERNALLY" (E15.8) | `input.json` (text); INPUT/GRIDS sections (%.6f or %.16e) |
| 2 grids | GRID PARAMETERS (I, Δφ 17 digits); MISC SETUP (F 6 decimals); `.pes` header records: zgrid/rgrid/xgrid/ygrid end points, θ and φ nodes (double) | GRIDS section (%.6f; dk %.16e; θ table %.16e); `asymptotic_potential.dat` r (%.16e) |
| 3 MISC SETUP | F 6 decimals, integers 1-based | GRIDS/INPUT (%.6f, %.16e), indices 0-based |
| 4 potential | `.pes` (V in hartree before −dt, symmetry-reduced, double) | not written: helper `potential_check` (double) |
| 5 V_asy | stdout table (list-directed, 17 significant digits) | `asymptotic_potential.dat` (%.16e) and log table |
| 6 rovib states | stdout (e, erel F 6 dec.; psi D 7 digits; norm, ⟨r⟩ F 7 dec.); `fort.21` record 2 eviroi (double, hartree) | log table (%.16e, psi %.8e); `gas_phase_states.dat` (%.16e hartree) |
| 7 diffraction basis | stdout (I, warning text) | log ASYMPTOTIC ANALYSIS (I, warning text); labels in `bk_coefficients.bin` |
| 8 rotation basis | stdout (I); `fort.21` labels PntToRot, J, mJ (I) | log; `bk_coefficients.bin` labels J, mJ (I) |
| 9 rot-dif basis | stdout (I); `fort.21` labels n, m, enymo (double hartree) | log; `bk_coefficients.bin` labels n, m, order, channel energies (double hartree) |
| 10 optical potentials | stdout coefficients (E15.8); indices in MISC SETUP | log CAP strengths (%.8e); `bk_coefficients.bin` restart parameters "…OptFactor (effective)" (double) |
| 11 initial packet | stdout WAVE FUNCTION (E15.8, I); `fort.21` record 2 kz, z0, α (double), record 3 E_trans | log ROTATIONAL BASIS AND INITIAL WAVE PACKET (%.16e); INPUT kz, α, E_trans (%.16e) |
| 12 analysis set-up | stdout ANALYSIS (E15.8, I); per-energy headers (D14.8) | log ASYMPTOTIC ANALYSIS; `scattering_summary.tsv` E_coll, E_total, E_z (%.16e) |
| 13 averages | `averages.out` (E15.8) | `averages.dat` (%.16e) |
| 13 analysis line | stdout Z (E15.8), step pinned via `fort.21` | log Z and t |
| 13 BK | `fort.21` (double) | `bk_coefficients.bin` (double) |
| 13 ψ final | `fort.9` (double; record 3 = surface state) | `checkpoint.bin` (double) |
| 14 probabilities | stdout per-energy blocks (E14.7, D14.7, F 7 dec., av. energies F 6 dec.) | `scattering_summary.tsv`, `probabilities_rovib.tsv` (%.16e) |
| 15 diffraction | stdout "probabilities as function of diffraction order" (E14.7) | `probabilities_diffraction.tsv` (%.16e) |
| 16 S-matrix | `fort.25` (Re/Im/\|S\|² D13.6, phase F10.6) | `s_matrix.tsv` (%.16e) |
| 18a fine analysis | (C++ `onlyAnalysis` on `fort.21`) | (C++ `onlyAnalysis` on `bk_coefficients.bin`) |

No counterpart: Fortran nmoder, nstatw, Δφ (C++ prints nφ only), per-v average energies,
"total momentum", the CPU TIMES breakdown; C++ θ weights, symmetry self-check, memory estimate,
the fixing time t of the analysis line.
