# Contributing

How to work on the code: the build variants, formatting and linting, the test suite in detail and
the repository layout. Installing, building and running the program: the [README](../README.md).

Bug reports and questions: [GitHub Issues](https://github.com/ewfsmeets/Cpp-SPO-DVR/issues).

## Build variants

The development builds live side by side in git-ignored directories.
`tools/build_all.sh [variant ...]` configures and builds them (all by default; the CUDA ones only
when nvcc is found). It does not use PGO; for that, see
[performance.md](performance.md#profile-guided-optimisation-1116--less-time-on-the-medium-grid-2026-10-26):

| Directory | Options | Purpose |
|---|---|---|
| `build` | (defaults) | the default (strict) build: oneMKL, no FMA (fused multiply-add), LTO |
| `build-fftw` | `SPODVR_FFT_LIBRARY=FFTW` | oneMKL LAPACK with the real FFTW |
| `build-nomkl` | `SPODVR_USE_MKL=OFF` | system LAPACK + FFTW: the fallback that must keep working |
| `build-poison` | `SPODVR_POISON_UNINITIALISED=ON` | NaN-fills new large arrays: finds reads before the first write |
| `build-fast` | `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fast CPU build (results differ at the 1e-12 to 1e-11 level; add PGO by hand for the fastest) |
| `build-cuda` | `SPODVR_USE_CUDA=ON` | GPU backend, no FMA (results comparable with the CPU build) |
| `build-cuda-fast` | `SPODVR_USE_CUDA=ON`, `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fastest GPU build (FMA in the kernels; no PGO) |

- Native restart files only work between builds with the same `SPODVR_FAST_MATH` setting
  ([README: restart](../README.md#continuing-a-run-restart)).
- Developer options besides those in the README's table: `-DSPODVR_BUILD_TOOLS=OFF` (default `ON`)
  skips the developer tools (operator timing, wave dump and compare); `-DSPODVR_WARNINGS_AS_ERRORS=ON`
  (default `OFF`) builds with `-Werror`; `-DSPODVR_CUDA_SYNC_DEBUG=ON` checks every CUDA kernel right
  after its launch ([cuda.md](cuda.md#building)); `SPODVR_PGO_DIR` sets the directory of the PGO
  profile ([performance.md](performance.md#profile-guided-optimisation-1116--less-time-on-the-medium-grid-2026-10-26)).
  The build is warning-free: keep it so.

## Formatting and linting

Git pre-commit hooks check every commit: clang-format (C++), ruff (Python: lint and format,
settings in `pyproject.toml`, `[tool.ruff]`), pymarkdown (Markdown lint, settings in `pyproject.toml`, `[tool.pymarkdown]`) and the
standard checks of pre-commit-hooks (trailing whitespace, end of file, YAML and JSON syntax,
merge-conflict markers); clang-tidy is run on demand. All versions are pinned in
`.pre-commit-config.yaml` and installed by pre-commit itself (no system packages needed):

```bash
pipx install pre-commit                                      # once per machine
pre-commit install                                           # once per clone
pre-commit run --all-files                                   # format and lint everything
pre-commit run clang-tidy --hook-stage manual --all-files
```

- The clang-tidy run needs `build/compile_commands.json`.
- clang-format, ruff and the whitespace hooks fix what they can and then stop the commit: stage the
  changed files and commit again. pymarkdown only reports.
- Markdown needs a blank line around headings, code blocks and lists (not around a code block inside
  a list item), a language on every code block (`text` for program output), and no inline HTML other
  than comments, `<a>`, `<br>`, `<details>`, `<summary>`, `<sub>` and `<sup>`. Outside code, write a
  literal `<` as `\<`.

## Testing

The commands for the whole suite, the quick tier and the GPU tests, with their run times and test
counts (measured on 2026-11-19, with the plot tool's tests), are in the [README](../README.md#testing). Run ctest serially (no
`-j`), with `OMP_NUM_THREADS` set to the number of physical cores.

- **Quick tier:** `ctest --test-dir build -E 'RegressionPropagation|PlotRuns'` runs everything but
  the propagation gates (gates 3 and 4 against the Fortran, 23–86 s each; the regular expression
  also leaves out their GPU twins `CudaRegressionPropagation.*`) and the plot tool's runs and the tests
  that read them (`PlotRuns.*`, about 160 s for the runs). Run it while working on a change;
  run the full suite before every commit that touches `src/`, `tests/`, the build files or reference
  data (the quick tier is enough for a commit of docs or Python tools only; for a change to
  `tools/plots` add `ctest -R '^(PlotRuns|Plots)\.'`), and always before merging into `main`.
- **Worktrees:** a fresh checkout or `git worktree` needs its own `uv sync` before the configure with
  tests (or `-DPython3_EXECUTABLE=` the main checkout's `.venv/bin/python`).
- **Partial builds:** with only some test targets built, a plain `ctest` counts every unbuilt target
  as a failure. Select with `-R`, or build everything.
- **Python tests:** every Python test of the repository runs under pytest in the repository's
  environment (`.venv/`, `uv sync`), one ctest per test module or test function, label `python`
  (`ctest --test-dir build -L python`). Each tool has its own `pytest.ini`; run one tool per pytest
  process (each tool has its own `pytest.ini` and import path):

  ```bash
  uv run python -I -B -m pytest -c tools/compare/pytest.ini --rootdir=tools/compare tools/compare \
      --compare-wave-tool=build/tools/spodvr_compare_wave
  uv run python -I -B -m pytest -c tools/fortran_reference/pytest.ini \
      --rootdir=tools/fortran_reference tools/fortran_reference
  uv run python -I -B -m pytest -c tools/wave_packet_designer/pytest.ini \
      --rootdir=tools/wave_packet_designer tools/wave_packet_designer
  uv run python -I -B -m pytest -c tools/plots/pytest.ini --rootdir=tools/plots tools/plots \
      --fixture-runs=build/tests/plot_fixtures
  ```

  Registration (`add_pytest_test` in `tests/CMakeLists.txt`): the designer is registered per module, so
  a new test function needs nothing and a new module one line; `tools/compare` and
  `tools/fortran_reference` per test function (their ctest names predate pytest), so a new function
  needs one line; the plot tool's modules are found by a glob (`tools/plots/tests/test_*.py`), so a new
  module needs nothing; a registration inside an `if()` also needs a `note_pytest_not_built` line in its
  `else()`. `Python.EveryTestRegistered` fails while a collected test has no ctest or two.

### What is tested

- unit tests, one program per module (`tests/test_<module>.cpp`);
- regression tests against the original Fortran code (`tests/test_regression_*.cpp`, test names
  `Regression…`/`Gate…`): potential, gas-phase states, propagation (norms and averages), on-the-fly
  analysis (probabilities, BK coefficients) and fine analysis. The Fortran reference output is
  stored in `tests/reference/` (the CRP and SCM driver outputs in `tests/crplib/reference/` and
  `tests/scm/reference/`). The gates, their tolerances and cases: [validation.md](validation.md);
- GPU tests (`tests/cuda/`, names `Cuda…`): only in builds with `-DSPODVR_USE_CUDA=ON`. They skip on
  a machine without a GPU;
- the wave-packet designer (`tools/wave_packet_designer/`, Python, names `WavePacketDesigner.…`): its
  copies of the program's rules against the C++ (input keys, checks, grid indices, packet values), its
  design rules and its figure. They run in the repository's Python environment (`.venv/`, made by
  `uv sync`; [README: Dependencies](../README.md#dependencies)).
  Two C++ tests (`MemoryEstimate.MatchesTheWavePacketDesignerCases` and its `CudaBackend`
  twin) check the program's grid indices and memory estimates against the designer's
  `expected_cases.json`: when they fail, the designer's copy needs updating
  ([how](../tools/wave_packet_designer/README.md#the-drift-guard-expected_casesjson)).
- the plot tool (`tools/plots/`, Python): `Plots.<module>` for its tests without run data, and
  `PlotRuns.<module>` for those that read real runs: the C++ test `PlotRuns.Fixtures`
  (`tests/test_plot_fixtures.cpp`) writes small runs with the current program into
  `build/tests/plot_fixtures` (ctest runs it first), so a changed output format fails the tool's
  readers ([tests](../tools/plots/README.md#tests)).

### Selecting tests

ctest `-R`/`-E` take a regular expression on the test name; `ctest --test-dir build -N` lists them
all:

```bash
ctest --test-dir build -E 'RegressionPropagation|PlotRuns'   # quick tier: no propagation gates, no plot runs
ctest --test-dir build -E Regression              # everything except the Fortran regression tests
ctest --test-dir build -R Gate3                   # gates 3 and 4 (propagation, on-the-fly analysis)
./build/tests/test_input --gtest_filter='Input.*'
```

- The last line runs one test program directly (GoogleTest's own filter).
- The wave-packet designer's tests take about 47 s (`ctest -R WavePacketDesigner`, 2026-11-19).
- **Suites run as one test:** ctest runs most GoogleTest cases as tests of their own, but five
  suites whose cases share an expensive set-up (a function-local static: program runs, potential
  grids) run as one test each, in one process, so the set-up is computed once:
  `Restart` and `FineAnalysisRuns` (`tests/test_restart.cpp`; its other suites stay one test per
  case), `RegressionPes` (every case of `tests/test_regression_pes.cpp` but the negative tests) and
  `RegressionPesScmNegative`, and `PlotRuns.Fixtures` (`tests/test_plot_fixtures.cpp`: the runs that
  the plot tool's tests read). Such a test fails when one of its cases fails or skips, or when its
  filter selects no case. To run one case, run the program with GoogleTest's filter, e.g.
  `./build/tests/test_regression_pes --gtest_filter='Scm/*'`; `-R Gate1` selects only the hard
  wall's gate-1 tests, `-R RegressionPes` the PES gates. A new suite with such a set-up gets an
  `add_one_process_test` line in `tests/CMakeLists.txt`.

### Long tests

The full 7000-step regression runs (hard wall and DZ, `cu111_dailight_small`; the latter alone takes
about 10 minutes on 8 cores), the negative tests of the low-Z and the SCM gates (1000- and 3000-step runs) and the CRP
cluster-pruning check on ten million points (`CrpLib.ClusterPruningLong.TenMillionPoints`) are
skipped unless `SPODVR_LONG_TESTS` is set:

```bash
SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build -R 'Long($|\.)' --output-on-failure
SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R 'Cuda.*Long($|\.)' \
  --output-on-failure
```

The regular expression selects the names that end in `Long` (or whose suite does); a plain
`-R Long` also runs the short `ScmFiles.ScmD6LongPathsAreNotTruncated`.

### Skipped regression tests and reference data

- A regression test skips (it does not fail) when its Fortran reference file is missing; the tests
  of the SCM cases `scm_cu111_srp48_small*` fail instead, and so do the PES gates, which ctest runs
  as one test that fails on a skip ([Suites run as one test](#selecting-tests)). A missing optional
  file (the BK file, the restart output) only drops that comparison.
- The scripts in `tools/fortran_reference/` regenerate the files; they need a Fortran reference
  install ([tools/fortran_reference/README.md](../tools/fortran_reference/README.md); how to get the
  Fortran code: [there](../tools/fortran_reference/README.md#getting-the-fortran-code)). Running the
  tests needs only the committed reference files.
- `ctest -R ReferenceManifest` fails when a file listed in a case's `MANIFEST.tsv` is missing or
  edited, or a case has no `MANIFEST.tsv` ([MANIFEST.tsv](../tools/fortran_reference/README.md#manifesttsv)).
- The test cases in `tests/reference/` include `hardwall_small` (also in `inputs/`), the calculation
  of `inputs/scm_example` (`scm_cu111_srp48_small`), the corrugated cases (`cu111_dailight_*`,
  `cu111_persson_small`) and the CRP cases (`crp_*`). Their `input.json` files run like any input
  ([README: Inputs](../README.md#inputs)).

### Other build variants

Run the same ctest command in their build directory ([Build variants](#build-variants)). Use
`build-poison` after touching the uninitialised-memory code. The fast-math builds (`build-fast`,
`build-cuda-fast`) are expected to fail a few tests that compare exact values; their gates pass. The
list, with the reason for each: [performance.md](performance.md#the-fast-math-build-spodvr_fast_math).

### Docs consistency

`python3 -I -B tools/check_docs.py [--spodvr build/spodvr]` (under 1 s) checks the docs against the
code: input keys, output files, links and anchors, decision IDs, local paths, command-line options
(with `--spodvr`), gates, environment variables, and the upstream commit of `docs/upstream_sync.md`
against the pin in `tools/fortran_reference/reference_config.sh` (and the install the runs use
against it); it also lists numbers repeated across docs (report-only). ctest runs it as
`Docs.Consistency` (`ctest --test-dir build -R Docs`). Accepted findings, each with a reason, are in
`tools/check_docs_allowlist.tsv`; a row that no longer suppresses anything fails the check.

### Numerically neutral changes and rounding

- Refactors must leave the results bitwise unchanged:
  `tools/check_bitwise_refactor.sh [--cuda] [--base-build DIR] [base-ref]` (default `HEAD`). It
  builds both trees in a scratch directory (or takes the base's builds from `DIR`), runs every run
  mode with both (potential only, a short propagation, a restart, an analysis-only run,
  `alsoAnalysis`, both checkpoint formats, a symmetry-reduced Dai–Light run, CRP and SCM set-ups)
  and compares every output file byte by byte, the printed output without paths, thread counts and
  timings, and `run_info.json` without the git commit and paths; `--cuda` repeats the propagating
  cases on the GPU. Its header lists the cases.
- For changes that may alter rounding: `tools/compare_propagation.sh [base-ref] [tolerance]` (its
  comparison, `tools/compare/compare_averages.py`, is tested by `ctest -R ComparePropagation`).
- Both build through `tools/scratch_build.sh`; to configure without network access, set
  `FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON` to the nlohmann_json sources.
- The final wave packet of two builds, or of the CPU and the GPU: `spodvr_dump_wave` and
  `spodvr_compare_wave` ([performance.md](performance.md)).

### Additions to the time loop

Everything added to the propagation loop (a reduction of ψ, a copy, a file write) is switched on by an
input key, and switched off it costs nothing measurable: no allocation, no backend call, only a test per
analysis step or averages row. Two checks show it:

- `TimeLoopCalls.*` (`tests/test_time_loop_calls.cpp`; on the GPU `CudaTimeLoopCalls.*`) runs the
  program through `run_program` with a backend that counts its calls (`tests/counting_backend.hpp`) and
  asserts each count, without the plot keys (the counts of the plain time loop) and with each of them. A
  hidden pass over ψ fails it, whatever the timing noise; a new key adds its case.
- `tools/benchmark/compare_speed.sh --base-build DIR --new-build DIR [--cuda] [--steps N]
  [--repetitions R]` runs two existing builds alternately without plot keys (5 runs each by default)
  on the benchmark's medium grid (8 threads, 50 steps) or, with `--cuda`, on the thesis DZ input on the
  GPU (40 steps). It passes when the new build's median "Propagation time" is at most max(1 %, the
  spread of the runs) slower and its host peak memory at most 16 MiB larger; GPU memory is not measured
  (`tools/compare/compare_speed.py`, tested by `ctest -R CompareSpeed`). It builds nothing: give it the
  build of the base commit and the new one.

### Production-grid validation

Beyond the gates, the full thesis calculation (`inputs/thesis_fig5_3_dz`, 7000 steps) was run with
the Fortran, the C++ CPU and the C++ GPU code and compared state by state:
[validation_thesis_fig5_3_dz](validation_thesis_fig5_3_dz/README.md). It is regenerated by its own
scripts, not by ctest.

## Repository layout

| Path | Content |
|---|---|
| `src/` | the program: `input/`, `grids/`, `potential/`, `symmetry/`, `gasphase/`, `propagation/`, `analysis/`, `io/`, `fft/`, `linalg/`, `core/`, `report/`, `app/`, `cuda/` (optional GPU backend), `crplib/` (the CRP library port), `fortran_text/` (the Fortran list-directed reader shared by `crplib/` and the SCM code) |
| `tests/` | unit tests (one file per module) and regression tests against Fortran reference output (`tests/reference/`, `tests/crplib/reference/`, `tests/scm/reference/`) |
| `tools/` | Fortran reference scripts, benchmark and comparison tools, the docs consistency checker, the wave-packet designer (`wave_packet_designer/`, Python), the plots of a run (`plots/`, Python) |
| `docs/` | the documentation (the [README's table](../README.md#documentation)), including the production-grid validation record (`validation_thesis_fig5_3_dz/`) |
| `inputs/` | example inputs |
| `potentials/` | data files of the Fortran potentials, one directory per system (verbatim copy of `potential.dir` of the upstream Fortran repository) |
