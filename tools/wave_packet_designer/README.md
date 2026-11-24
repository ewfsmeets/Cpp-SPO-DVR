# Wave-packet designer

A standalone Python tool for setting up the Z and r part of an input: the Z grids, the initial wave
packet, the absorbers and the r grid. It tells you, before any run, whether `spodvr` accepts the input
and what the packet, grids and absorbers look like. It needs no build of the program: it carries copies
of the program's rules (marked as such, see [Keeping up with the C++](#keeping-up-with-the-c)).

The physics behind the choices (packet width, momentum margins, absorber lengths, …) is explained in
[inputs/README.md](../../inputs/README.md); this tool applies it.

## Running

Two ways, both from inside the repository (uv finds its `pyproject.toml` upward from where you start; paths
are relative to where you start them):

- **The window** (matplotlib): the figure below with sliders for the recipe's knobs. It needs
  [uv](https://docs.astral.sh/uv/) (`pipx install uv`, like pre-commit in docs/CONTRIBUTING.md), which
  installs numpy, matplotlib and PyQt6 from the header of `gui_mpl.py` the first time. Under WSL the window
  opens through WSLg.
- **The command line**: text reports; runs in the repository's Python environment (`uv sync` once,
  [README: Dependencies](../../README.md#dependencies); numpy, and matplotlib for `check --png`).

```bash
uv run tools/wave_packet_designer/gui_mpl.py                                   # the DZ input
uv run tools/wave_packet_designer/gui_mpl.py inputs/thesis_fig5_3_er/input.json
uv run python -I -B tools/wave_packet_designer/wave_packet_designer.py check inputs/thesis_fig5_3_dz/input.json
```

![The designer window on the DZ input: six panels, the verdict, the sliders, the buttons and the help box](images/window.png)

The window on the DZ input (the first command above), with the mouse over the Z panel: the help box at
the bottom right explains that panel. The memory is the worst case, because the command gives no
`--representatives`; the RAM and the free GPU memory are those of the machine it ran on. To refresh the
picture after changing the window, open it the same way, keep the mouse over the Z panel and press `s`
(matplotlib's save key); the empty strip below the buttons was cropped.

The window shows six panels, a verdict and a help box:

- **Z**: the two Z grids as lanes (the scattering grid with its interaction region and absorber, the
  specular grid with its absorber and the packet's allowed region), and the initial packet |ψ(Z)|².
- **Z momentum**: |ψ(k)|² against the grid's limits ±π/dZ, with the outward fraction of the norm.
- **Absorbers**: reflection R and transmission T of each absorber over the energies it must absorb, with
  the levels of the README's examples.
- **r**: for each analysed vibrational level, the norm of the bound state beyond r, against the 1e-6
  limit and ROptStart.
- **Memory**: the peak memory on the CPU and the GPU against this machine's RAM and free GPU memory.
- **Time**: the classical paths Z(t) of the packet's fast, central and slow parts, in from Z0, reflected at
  a turning point Z_turn and, elastically, back out across Zana; and an outgoing wave of the lowest normal
  energy E_n,min that tEnd should converge (section 8). The slow part's elastic recrossing is a lower bound
  for tEnd; E_n,min (0.03 eV by default) sets what tEnd needs for the slow inelastic and diffractive
  channels.

The verdict (the panel at the bottom right) says what `spodvr` would do with the input, then each recommendation of
inputs/README.md it misses, worst first (the run fails, the results suffer, it costs more), with its
value, limit and what to change; then the grid sizes, the absorbers as placed and tEnd, and what nothing
here checks. Hovering over a panel, a slider or a button explains it in the help box; **Explain** lays the
explanations over the panels (the texts are in `explanations.py`).

The window opens on the input as read. The sliders are EMinWP, EMaxWP, dZ, Z0, ZSpOptStart, the lengths
of the Z absorbers, dR, nR, ROptStart, tEnd, the target energies of the absorbers tuned by energy, and
E_n,min; each keeps a black tick at the input's value. When dZ, ZSpOptStart or an absorber length moves, the Z grids follow as
the recipe places them (Zana on the first grid point at or beyond Z_asymptotic; nZ and nZSpecular the
smallest valid FFT sizes that hold the absorbers, which end on their last points); tEnd moves in whole
averaging intervals dtAve. **Suggest** applies the recipe, **Reset** returns to the input, and **Save**
writes the changed input to the file named in the box (never over an existing file). Options:
`--zasymp` and `--point-group` (CRP data sets), `--threads`, `--representatives` and `--operations`
(memory, as for `memory` below), `--ram` and `--gpu` (GiB, instead of this machine's), `--zturn` (bohr, the turning point of the time
budget; 2 by default, ZStart for the worst case).

| Command | What it does |
|---|---|
| `check <input.json>` | Runs the program's input checks and set-up stages on the input and prints the packet, the Z grids, the absorbers, the r grid, the gas-phase channels, the design checks of inputs/README.md, the time budget and the absorbers' R and T. |
| `check <input.json> --png figure.png` | The same, and the window's figure in a PNG file. |
| `check <input.json> --zasymp Z --point-group G` | The same for a CRP data set (SystemID `CRPlibrary`, `CRPDistrtSurface`), whose Z_asymptotic and point group come from the data set (the program's set-up output prints both). |
| `memory <input.json> [--threads N]` | The sizes that set the memory (θ and φ points, rotational basis, grid points) and the program's CPU and GPU memory estimates, array by array. `--threads` (default 8): the run's OpenMP threads, each with its own work arrays. |
| `memory <input.json> --representatives N [--operations N]` | The same for a run with `UseSymmetry`: the (x, y) points at which the program stores the potential, and the operations it uses (default: all of the point group). Its SURFACE SYMMETRY block prints both: "(x, y) points computed / total", and the names under "Operations used" (count them, identity included). Without them the estimate assumes no reduction, an upper bound. |
| `suggest <input.json> [--dz D] [--output new.json]` | The Z grids, Z0, Zana and the Z absorbers by the recipe of inputs/README.md (section 2): prints the suggested values next to the input's and checks the changed input; `--output` writes it to a new file (never over an existing one). `--dz` sets the step; by default the input's if it holds the momenta, else the largest that does (the local momenta and the time step of section 4 are not checked). |

Exit status: 0 `spodvr` accepts the input (for `suggest`: the changed input), 1 it would stop on it (or
`suggest` needs other keys fixed first), 2 the file or the options are unusable (also a CRP data set
without `--zasymp` and `--point-group`). The tool writes only the file that `suggest --output` names, and
never replaces an existing file, and the PNG that `check --png` names (replaced if it exists).

## What `check` reports

**Program checks.** What `spodvr` would do with the input, stage by stage as in its set-up: reading
(known keys, value types, the checks made while reading), the consistency checks (`validate_input`),
the grids, and the initial wave packet's TolSp check. Each finding names, in brackets, the C++ function
that raises it. Where `spodvr` stops at the first error, the tool lists every error of that stage.
Warnings the set-up prints (an absorber outside the strength table, a grid point the Fortran puts on the
other side of an absorber's end) are listed too.

**The input as the program sees it.** The packet (mean momentum kz, widths σ_k and σ_Z, the part
moving outward (the Gaussian tail of |ψ(k)|² at k > 0), the allowed region, the norm and the part of it
outside the region, |ψ(k)| at the grid edge), the Z grids (sizes, k_grid = π/dZ and the tail beyond it,
which wraps round to outward motion, the interaction region nzpot, the analysis line), each absorber
(grid indices, target energy, L/λ, strength), the r grid, and the gas phase: the open vibrational levels
and the highest open j, which set the analysed channels. An excerpt for the DZ input:

```text
INITIAL WAVE PACKET (Z)
-----------------------
  Energy range (EMinWP .. EMaxWP) .............. 0.4 .. 0.95 eV
  Mean momentum kz (towards the surface) ....... 13.2052 1/bohr
  Width in k, sigma_k = sqrt(alpha) ............ 1.4060 1/bohr
  Width in Z, sigma_Z = 1/(2 sigma_k) .......... 0.3556 bohr
  Centre Z0 .................................... 18 bohr
  Packet moving outward, k > 0 ................. 2.9e-21
  Allowed region (specular indices) ............ 125 .. 183
  Allowed region (Z) ........................... 15.0000 .. 21.9600 bohr
  Norm on the specular grid .................... 0.9999999999999998
  Norm outside the region, norm - norm inside .. 0.000e+00
  1 - norm inside (must be <= TolSp) ........... 2.220e-16 (TolSp 1e-07)
  |psi(k)| at the grid edge / maximum .......... 5.674e-10
```

"1 − norm inside" is what the program compares with TolSp; here it is the rounding of the whole norm,
not norm outside the region.

The gas phase is computed as the program does: the potential's V(r) at Z_asymptotic on the input's r
grid, solved for the input's masses and JMax. For the Dai-Light and hard-wall potentials V(r) is copied
from the C++ (`program/gas_phase_potentials.py`), so the states are the program's; a second solution on a
four times finer r grid gives the converged states, and where the two differ the r grid is too coarse. For
another potential (Persson, the CRP data sets) the DZ input's printed V(r) (`dz_gas_phase.py`) stands in
and the report marks the states as indicative.

**Design.** The recommendations of [inputs/README.md](../../inputs/README.md) applied to the input, each
with its value, the recommended limit and the README section: the momentum margins of the packet and of
dZ, Zana on the first grid point, Z0 and the packet's outer tail, the absorber lengths and ends, the r
absorber beyond the analysed bound states; for information, the lowest outgoing energy converged by tEnd
and the energy resolution. Then the time budget (section 8 of inputs/README.md, "TIME"): the time per
bohr of the slowest analysed part (the momentum of EMinWP, M / k_min), when that part, reflected
elastically, is back at Zana, and the tEnd that converges the lowest outgoing normal energy. Last R + T,
R and T of each absorber over the energies it must absorb (1D transfer matrix through the exact −iV_opt,
with the program's strengths; from 0.05 eV or EMinWP up). These are advice: they never change the exit
status.

## What `memory` reports

The program's MEMORY ESTIMATE and GPU MEMORY ESTIMATE tables, with the same rows and sizes: for the DZ
input with `--representatives 19` they equal the program's (7.9 GiB while setting up, 7.3 GiB during the
propagation, 8.8 GiB on the GPU). The CPU estimate covers run modes 0, 1 and 6 (a potential-only run
holds the potential only), the GPU estimate the propagating runs.

## Architecture

The code is layered; each layer imports only from the layers below it (checked by
`tests/test_structure.py`):

```mermaid
flowchart TB
    cli["wave_packet_designer.py<br/>command line"]
    window["gui_mpl.py<br/>window (matplotlib)"]
    plots["plots.py<br/>the figure, one function per panel"]
    report["report.py<br/>text tables"]
    evaluation["evaluation.py<br/>evaluate(parameters, settings) -> Evaluation"]
    program["program/<br/>copies of what spodvr computes or rejects"]
    design["design/<br/>advice from inputs/README.md"]
    data["settings.py, input_file.py, dz_gas_phase.py<br/>options; input files; the gas-phase curve"]
    cli --> report --> evaluation
    window --> plots --> report
    plots --> evaluation
    window --> evaluation
    window -- "input files, settings" --> data
    window -- "parse_input, key_of_field" --> program
    cli --> data
    cli -- "parse_input" --> program
    report -- "findings, units" --> program
    evaluation --> design --> program
    evaluation --> program
    evaluation --> data
```

- **One computation.** `evaluate()` computes everything once and returns an `Evaluation`; the front end
  and the report only print it. A new quantity belongs in `evaluate()` (or `program/`), never in the
  presentation. The one step before it: the front end reads the file (`input_file.py`) and its keys
  (`program/json_input.py`, `parse_input`), and calls `evaluate()` only when `spodvr` would get through
  the reading.
- **The program's rules live in `program/`.** Each module's docstring has a `Shadows:` line naming the
  C++ files and functions it copies. Most modules copy one C++ file and have its name
  (`validate.py` ↔ `src/input/validate.cpp`); to find the copy of any C++ file, search those lines:
  `grep -rn "spatial_grids.cpp" tools/wave_packet_designer`. `program/` imports only the standard
  library and numpy.
- **The parameters are the program's.** `program/config.py` has `InputParameters`, flat like the C++
  struct, one field per input key (the key and unit in the comment). A front end changes a value with
  `dataclasses.replace(parameters, dz=0.15)` and calls `evaluate()` again; it saves the change under the
  field's key (`key_of_field`, `input_file.save_changed_copy`).
- **Stages like the program.** `evaluate()` runs the set-up stages in the program's order; a stage runs
  only when the ones before it passed, and `Evaluation.stopped_at` names the stage where `spodvr` would
  stop. One stage is the designer's own: `settings`, when a CRP data set lacks its Z_asymptotic and point
  group, or the symmetry counts cannot hold for the grid (the command line then exits 2).

| Module | Responsibility | Shadows |
|---|---|---|
| `wave_packet_designer.py` | command line | — |
| `gui_mpl.py` | the window: sliders and buttons wired to `evaluate()` and `plots.py` | — |
| `plots.py` | the figure of an `Evaluation`, one function per panel | — |
| `explanations.py` | the texts of the help box and of Explain: per panel and per control | — |
| `report.py` | text report of an `Evaluation` | `src/report/format.cpp` (format_bytes), `sections.cpp` (memory table) |
| `evaluation.py` | `evaluate()`: the stages, the absorber and packet summaries, the memory, the gas phase, the design checks and R + T; `suggest()`; `fitted_parameters()` and `fitted_time()` for the window | `src/app/run.cpp` (stage order) |
| `settings.py` | values that are not input keys (Z_asymptotic of a CRP data set, …) | — |
| `input_file.py` | read input.json; save changed values into a copy of its text | — |
| `design/rules.py` | the recommendations of inputs/README.md, one function each | — (README sections) |
| `design/transmission.py` | R and T of an absorber (1D transfer matrix) | — |
| `design/timing.py` | the time budget: elastic recrossing, E_n,min, the paths of the time panel | — (README section 8) |
| `design/suggest.py` | the recipe of README section 2 for the Z keys | — |
| `program/config.py` | `InputParameters`, `DerivedParameters` | `src/input/config.hpp` |
| `program/json_input.py` | key table, value types, checks while reading | `src/input/json_input.cpp` |
| `program/validate.py` | the consistency checks | `src/input/validate.cpp` |
| `program/process.py` | derived parameters: masses, absorber strengths, kz, α; the absorber-table warning | `src/input/process.cpp`, `src/report/sections.cpp` |
| `program/wave_packet_params.py` | strength table, L/λ, α | `src/input/wave_packet_params.cpp` |
| `program/spatial_grids.py` | grids, nzpot, absorber indices, the grid-point trap | `src/grids/spatial_grids.cpp`, `uniform_grid.hpp`, `momentum_grid.cpp` |
| `program/initial_wave_packet.py` | ψ(Z), \|ψ(k)\|, the TolSp check | `src/propagation/initial_wave_packet.cpp` |
| `program/angular_grid.py` | θ and φ points, the J range, the rotational basis size | `src/grids/angular_grid.cpp`, `rotational_basis.cpp` |
| `program/memory_estimate.py` | the CPU and GPU memory, array by array | `src/propagation/memory_estimate.cpp`, `src/cuda/cuda_backend.cpp`, … |
| `program/gas_phase.py` | the FGH states of V(r), the open channels | `src/gasphase/fourier_grid_hamiltonian.cpp`, `gas_phase.cpp` |
| `program/gas_phase_potentials.py` | V(r) of the Dai-Light and hard-wall potentials at Z_asymptotic | `src/potential/dai_light_cu111_potential.cpp`, `hard_wall_potential.cpp`, `potential_grid.cpp` |
| `program/point_group.py` | the φ multiple and order of each point group | `src/symmetry/point_group.cpp` |
| `dz_gas_phase.py` | the DZ input's printed V(r): the stand-in for other potentials, and test data | — |
| `expected_cases.json` | the drift guard: thirteen cases with the designer's grid indices, angular sizes and memory estimates (written by `write_expected_cases.py`) | read by `tests/designer_cases.hpp` for the C++ tests |
| `program/potential_systems.py` | Z_asymptotic, point group and (for the analytic potentials) the lattice constants per SystemID | the analytic potentials, `src/potential/potential_factory.cpp` |
| `program/fft_sizes.py`, `real_compare.py`, `unit_conversion.py` | valid FFT sizes; the Fortran's tolerant comparisons and NINT; units and π | `src/fft/fft_sizes.cpp`, `src/core/real_compare.hpp`, `src/core/unit_conversion.hpp`, `src/core/constants.hpp` |
| `program/findings.py` | an error or warning of the program | `src/core/error.hpp` |

### Keeping up with the C++

The copies in `program/` must follow the C++. They reproduce its numbers bitwise where the C++ uses
plain double arithmetic (same operations in the same order: e.g. `x / 27.2114`, not `x * (1/27.2114)`;
NINT rounds half away from zero, unlike Python's `round`). To change …

- **a C++ file named in a `Shadows:` line:** edit the Python modules whose `Shadows:` lines name it
  (grep, see above). The tests in `tests/` compare with the program:
  - `test_validate.py`: the key and message of every check, in order, with `validate.cpp`; the analytic
    potentials' lattice constants with the literals of their C++ sources;
  - `test_input_reading.py`: the key table with `json_input.cpp`;
  - `test_packet.py`, `test_grids.py`, `test_absorbers.py`: the set-up output of the DZ run (masses, α and kz bitwise, the grid
    indices, the packet's norms and |ψ(k)|, the absorber strengths) and the L/λ the Fortran printed;
  - `test_memory.py`: every row of the DZ run's two memory tables, and the hand counts of
    `tests/test_memory_estimate.cpp` (the repository's); `test_expected_cases.py`: the drift guard (below);
  - `test_gas_phase.py`: the DZ run's state energies (to 1e-13 eV) and open channels;
  - `test_design.py` checks the rules against the numbers inputs/README.md prints (sections 3, 4, 7, 8).

  The command of the DZ run is in `tests/helpers.py`; the Fortran values come with the C++ test that
  the test's docstring names.
- **an input key:** a field in `program/config.py` and a row in `program/json_input.py`
  (`test_input_reading.py` fails until both exist and match `json_input.cpp`).
- **a consistency check:** a line in `program/validate.py` at the same place as in `validate.cpp`, and a
  violating row in `tests/test_validate.py` (which fails until both exist).
- **what the report shows:** `report.py` (and `evaluation.py` if a new quantity is needed).
- **a recommendation of inputs/README.md:** a function in `design/rules.py`, added to `RULES`, and a failing
  row in `tests/test_design.py` (which fails for a rule without one). Its message cites the README
  section; the physics stays in the README.
- **a renamed C++ function:** the `Shadows:` line (`tests/test_structure.py` fails otherwise).

#### The drift guard: `expected_cases.json`

Most comparisons above use numbers recorded from the program once, so they do not notice a later change
of the C++. For the grid indices, the angular sizes and the memory estimates, which change whenever the
program gains or resizes an array, `expected_cases.json` closes that gap. Each of its thirteen cases holds a
definition written by hand (a base input, key changes, the threads, and the symmetry counts the program
computes) and the designer's results, and both sides are tested against it:

- `tests/test_expected_cases.py` checks the designer, and that a planted slip (plain `>=` for the tolerant
  comparison that places nzpot) fails exactly the one case built for it;
- `MemoryEstimate.MatchesTheWavePacketDesignerCases` (`tests/test_memory_estimate.cpp`) checks the program's
  grid indices, angular sizes, symmetry counts and CPU estimate, and
  `CudaBackend.MemoryEstimateMatchesTheWavePacketDesignerCases` (`tests/cuda/`, CUDA builds only) its GPU
  estimate.

When one of the C++ tests fails, the program has changed something the designer copies; the message names
the case and the field. Update the copy in `program/` (its `Shadows:` lines say where), then regenerate the
file and run both sides (`-R WavePacketDesigner` matches the C++ tests too):

```bash
uv run python -I -B tools/wave_packet_designer/write_expected_cases.py
OMP_NUM_THREADS=8 ctest --test-dir build -R WavePacketDesigner --output-on-failure
```

`write_expected_cases.py` keeps each case's definition and recomputes the rest. To add a case, write its
definition (for `UseSymmetry` the counts from the run's SURFACE SYMMETRY block, or from the C++ test's
failure), then run it.

### Conventions

- Write for a reader who knows Python and [inputs/README.md](../../inputs/README.md): short functions
  with one job, a docstring saying what and why (and the units), descriptive names (the input key in a
  comment where it differs), plain loops where they read better than numpy.
- Atomic units inside; eV, amu and degrees only where the input keys have them.
- Type hints on every function (optional local check: `uvx mypy tools/wave_packet_designer`).
- Run with `uv run python -I -B` (the repository's environment); entry files set
  `sys.dont_write_bytecode` and put their own directory on `sys.path`. Formatting and lint: ruff (pre-commit, settings in the repository's `pyproject.toml`), lines up
  to 120 characters.

## Tests

```bash
OMP_NUM_THREADS=8 ctest --test-dir build -R WavePacketDesigner --output-on-failure
uv run python -I -B -m pytest -c tools/wave_packet_designer/pytest.ini \
    --rootdir=tools/wave_packet_designer tools/wave_packet_designer                 # all, directly
uv run python -I -B -m pytest -c tools/wave_packet_designer/pytest.ini \
    --rootdir=tools/wave_packet_designer tools/wave_packet_designer/tests/test_plots.py   # one module
```

The tests run under pytest (`pytest.ini`, `tests/conftest.py`). Each module in `tests/` holds small
`test_…()` functions, each with a docstring saying what it pins down (a test that needs a scratch directory
takes `work`, pytest's `tmp_path` under the name the tests always used). ctest runs one module per test
(`WavePacketDesigner.<Module>`, registered in the repository's
[tests/CMakeLists.txt](../../tests/CMakeLists.txt)). To add a test, write the function: its module's ctest
runs it. A new module also gets a line there (the ctest `Python.EveryTestRegistered` fails until it has
one). Run the designer's tests in a pytest process of their own (its `pytest.ini` puts the tool's
directory on the import path, and its `tests` package is not unique). `tests/conftest.py` selects matplotlib's Agg backend (ctest
also keeps its cache in the build tree), so the module `test_plots.py` draws every panel and drives the
window's sliders and buttons without opening a window. Its window tests skip the idle redraws (on Agg
each one redraws the whole figure at once, and every slider change asks for one) and draw the states they
check (`draw_window` in `tests/test_plots.py`).

## Limits

- No surface potential is evaluated: Z_asymptotic and the point group of a CRP data set must be given, the gas
  phase of the Persson potential and the CRP data sets is the DZ input's stand-in, and the surface-symmetry
  reduction of the memory must be given
  (`--representatives`).
- No fine-analysis memory estimate (an analysis-only run).
- No plot data (`writeDensities`, `writePotential`, `PlotSites`, `dtFrames`, `FrameExtras`, `NormBudget`,
  `EnergyResolvedSites`, `EnergyResolvedEnergies`): their memory and disk use are not in the estimates; reading warns
  when a plot-data key is set.
- `write2DPES`: the bounds x2DPES ≤ alat and y2DPES ≤ blat are checked only for the three analytic potentials
  (the designer does not know the lattice of a CRP data set, or of an SCM slab); the other six checks always.
  The program writes the cut after the grids and stops: the memory estimates and the gas phase shown are
  those of the same input without write2DPES.
- Not checked, as the program does after the grids: that some operation of the point group is usable on
  the grid with `UseSymmetry` (`src/potential/potential_symmetry.cpp`), that the GPU has the memory of the
  GPU estimate (`--device cuda`, `src/app/backend_factory.cpp`), and the physical memory against the CPU
  estimate.
- Not copied:
  - the SCM key checks (`check_scm_keys`);
  - the flux-analysis input checks (a propagating run with `doFluxAnalysis` stops anyway; the tool
    reports that);
  - the restart-directory check;
  - the checks of the plot-data keys (`src/input/validate_outputs.cpp`, e.g. `writeDensities` only in a
    propagating run, the `PlotSites` values); the designer reads these keys (their types and the fields of a site)
    but does not model the plot data;
  - the checks of `make_angular_grid` (`src/grids/angular_grid.cpp`): those on JIn, JMax, MJIn and MJMax
    repeat checks of `validate_input`, and the two on the number of φ points hold by construction (it is
    `nearest_fft_size(2·MJMax+1)` for the φ multiple of the point group).
- |ψ(k)| comes from numpy's FFT, the program's from FFTW or oneMKL: equal to rounding.
