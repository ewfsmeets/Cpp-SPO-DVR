# Plots of a spodvr run

Figures and a report of one `spodvr` run, read from the files the run wrote: the run verdict (is the run
complete and converged?), the results against energy, the averages, the gas phase and the Balint-Kurti
coefficients. A run split into restart segments is one run: give every segment's directory and the plots
join them. The wave-packet designer ([tools/wave_packet_designer](../wave_packet_designer/README.md)) checks
an input before a run; this tool looks at what the run produced.

## Running

Once per checkout, `uv sync` (it builds the repository's Python environment, `.venv/`;
[README: Dependencies](../../README.md#dependencies)). Then, from the repository root:

```sh
uv run python -I -B tools/plots/spodvr_plot.py list                       # the figures and the files each needs
uv run python -I -B tools/plots/spodvr_plot.py verdict RUN...             # the verdict as text
uv run python -I -B tools/plots/spodvr_plot.py plot results RUN... --out DIR
uv run python -I -B tools/plots/spodvr_plot.py report RUN... --out DIR    # report.md + every figure
```

- `RUN...` is the run's output directory, or every segment's directory of a restarted run, in any order.
  Segments that do not join (one does not start where another ended) are refused; a segment whose
  RestartDir names another directory gives a `warning:`. The last segment alone also works: its
  bk_coefficients.bin holds the whole run, its averages.dat only that segment.
- Each figure is written as `DIR/<name>.png` (no metadata: reproducible bytes) and `DIR/<name>.pdf`. The report
  opens with the verdict, then shows every figure whose files the run has, each with its title and what good
  looks like (`explanations.py`), and lists the figures it could not draw and why.
- `plot` and `report` write only into `--out`, and there replace only files the tool wrote before (listed in
  `DIR/written_by_spodvr_plot.json`): any other file of the same name stops them before they write anything.
- Runs made before spodvr wrote `run_info.json` are refused.

Exit status: 0 done (the verdict is advice and never changes it); 1 the report could not draw a figure or make
the verdict (report.md says which); 2 bad arguments, an unreadable run, a figure whose files the run does not
have (`plot`), or a file in `--out` that the tool did not write (`error: <file>: <what>; <fix>`). A traceback
is a bug of the tool.

| Figure | Shows | Needs |
|---|---|---|
| `verdict` | the run verdict as a text panel (the report prints it as text) | whatever the run has |
| `results` | P_scattered and P_reaction, P(v') (when two or more levels are scattered into) and the mean final energies against the collision energy | scattering_summary.tsv |
| `averages` | norm, 1 − norm, rest norm, ⟨Z⟩ and ⟨r⟩ against time | averages.dat |
| `gas_phase` | V_asy(r) with the levels of the initial j and the r absorber; \|ψ(r)\| of the lowest levels on a log axis | gas_phase_states.dat, asymptotic_potential.dat |
| `bk` | \|c(t)\|² of the strongest channels; P_scattered(E; t) / P_scattered(E; tEnd) | bk_coefficients.bin (one record at least) |

Every figure also reads `run_info.json` (the run's metadata). The file formats:
[docs/output_files.md](../../docs/output_files.md).

### The verdict

Each finding has a severity, its value against a limit, what it means and the section of
[inputs/README.md](../../inputs/README.md) that explains it, worst first. Thresholds are named constants in
`compute/diagnostics.py`, each marked *derived* (from the physics or inputs/README.md) or *heuristic* (a
judgement); `verdict.py` reads the run's files for them.

| Finding | Severity | When |
|---|---|---|
| run incomplete | error | averages.dat, scattering_summary.tsv or s_matrix.tsv missing: the program writes them after the propagation, so the run stopped early or is still running (bk_coefficients.bin, appended every analysis step, shows how far it got). Then no other check runs (but "propagation diverged") |
| propagation diverged | error | a NaN or inf in averages.dat or bk_coefficients.bin |
| unanalysed flux, upper bound | warning | more than 1e-3 of the norm still inside the analysis lines at tEnd (`checkpoint.bin`), summed over all energies |
| P_scattered still moving | warning | P_scattered changed by more than max(1e-3 P_s, 1e-4) over the last 10 % of the run at some energy: the check per energy |
| unitarity error | warning | P_reaction below −1e-4 |
| slowest converged energy | warning above 1e-3, else information | E_n,min of inputs/README.md section 8 (Z_turn 2 bohr), and the probability already in the channels that leave with a normal energy below it |
| norm left | information, always | the norm at tEnd: it includes flux already analysed but not yet absorbed, and bounds only the error averaged over the packet's spectrum |
| not checked: … | information | a check the run's files do not allow: no checkpoint.bin (`CheckpointFormat` `fortran`), fewer than two BK records, or a run that does not propagate (potential only, analysis only) |

E_n,min is a lower bound on the time a run needs: "P_scattered still moving" is the empirical check.

## Architecture

The tool is plain numpy and matplotlib in layers; a module imports only from lower layers (or its own
package), which `tests/test_structure.py` checks. An arrow means "is used by":

```mermaid
flowchart LR
    copies["copies/ (0)"] --> readers["readers/ (1)"] --> compute["compute/ (2)"]
    compute --> style["style.py, explanations.py, verdict.py (3)"] --> panels["panels/ (4)"]
    panels --> figures["figures.py (5)"] --> report["report.py (6)"] --> cli["spodvr_plot.py (7)"]
```

| Module | Responsibility | Shadows (the C++ it follows) |
|---|---|---|
| `copies/s_matrix.py` | the S matrix from the BK coefficients: time sums, prefactor, phases, the packet's weight of k0 | `src/analysis/s_matrix.cpp` |
| `copies/real_compare.py` | the tolerant comparison `rne` | `src/core/real_compare.hpp` |
| `readers/errors.py` | `ReadError`, the one error the tool raises for a run it cannot read | — |
| `readers/tables.py` | averages.dat and the tab-separated tables (read by column name), the gas-phase files | `src/report/propagation_report.cpp`, `analysis_files.cpp`, `sections.cpp` |
| `readers/binary.py` | bk_coefficients.bin, checkpoint.bin (memory-mapped); the byte layouts as named struct formats | `src/io/checkpoint.cpp` |
| `readers/run_info.py` | run_info.json as `RunInfo` | `src/report/run_info_file.cpp` |
| `readers/segments.py` | restart segments in time order; records joined (a repeated record must be bitwise equal) | — |
| `compute/analysis.py` | a run's analysis energies, analysis time step and S-matrix factors, from run_info.json | — |
| `compute/results.py` | probabilities, mean energies and the unitarity error against energy | — |
| `compute/time_series.py` | BK coefficients per channel, P_scattered when the time sums stop at t | — |
| `compute/gas_phase.py` | the levels and which are bound | — |
| `compute/diagnostics.py` | the verdict's findings and their thresholds | — |
| `verdict.py` | the verdict of a run: reads what its files allow, a missing file is a finding | — |
| `style.py` | colours, rcParams, the scale rule, the info box, saving | the designer's colours (`tools/wave_packet_designer/plots.py`) |
| `explanations.py` | each figure's title and what it shows (printed under it in the report) | — |
| `panels/*.py` | one `draw_<panel>(ax, …)` per panel; draws (and masks below a floor), computes nothing | — |
| `figures.py` | one `figure_<name>(runs, info, options)` per entry, and `PLOTS`, the only list of entries | — |
| `report.py` | report.md with the verdict first; the overwrite guard of `plot` and `report` | — |
| `spodvr_plot.py` | the command line | — |

To change… edit:

- a threshold of the verdict: its constant in `compute/diagnostics.py` (the texts are formatted from it; the
  tests in `tests/test_diagnostics.py` fire just above each limit);
- which files the verdict reads: `verdict.py`;
- a colour, a font or the scale rule: `style.py` (the only place with colours);
- a figure's title or the text under it: `explanations.py`;
- a file format the program writes: the reader whose `Shadows:` line names that writer
  (`grep -rn 'Shadows:.*<file>' tools/plots`); the fixture runs (written by the current program) fail the
  readers' tests until it follows.

### Conventions

- Write for a reader who knows Python and the physics: short functions, type hints and a docstring on every
  function (ruff's D1 and ANN rules for this tool, set in the repository's `pyproject.toml`), descriptive names
  with a unit suffix on every dimensioned value (`_ev`, `_hartree`, `_bohr`, `_au`, `_per_bohr`; none =
  dimensionless); unit factors from `RunInfo`, never module constants.
- Readers, `copies/` and `compute/` never print or draw; panels never compute; figures are
  `matplotlib.figure.Figure` objects (no pyplot). These are kept by review.
- `tests/test_structure.py` checks: the layers; numpy only in layers 0–2 and matplotlib from layer 3 up (no
  other third-party package); functions of at most 40 code lines (except in `copies/`); a docstring on nested
  functions too; decorators only `@dataclass`, `@pytest.fixture` and `@pytest.mark.parametrize`; no classes
  but dataclasses and `ReadError`; no lambda or `partial` in `figures.py`; no `matplotlib.use`; no test
  module importing a name that starts with `Test`; every `Shadows:` line names an existing C++ file and a
  function it defines or calls.
- `copies/` reproduce the C++'s formulas and, where they are plain double arithmetic, its operation order; the
  time sums are a matrix product (the C++ adds term by term), equal to rounding.
- Values below a panel's floor are masked (not drawn); axes are logarithmic when a series spans more than
  1e3 (`style.scale_and_limits`); every axis is labelled with its unit; the info box under each figure names
  the run, its grids, the files read and the program version (a figure with notes asks `new_figure` for
  room for them).
- Errors: `ReadError("<file>: <what>; <fix>")`.

### How to add a plot

1. A reader and its test, only for a new file.
2. A compute function and a test against a hand-computed value or the program's files.
3. The panel: `draw_<panel>(ax, data…)` in `panels/<topic>.py`.
4. A `figure_<name>` in `figures.py` (the first line of its docstring is what `list` prints), its line in
   `PLOTS` (from there `list`, `plot` and `report` find it), and its title and text in `explanations.py`.
5. A content test (the texts, marks and curves the panel must show) in `tests/test_runs_figures.py`; if the
   figure needs a file that a stopped or potential-only run lacks, nothing else; otherwise add its name to
   those kinds in that test's `EXPECTED` table.

The averages figure, as `figures.py` has it (its info box also counts the segments):

```python
def figure_averages(runs: list[Path], info: RunInfo, options: argparse.Namespace) -> Figures:
    """The norms and mean positions against time, the segments joined."""
    averages, starts_au = read_joined_averages(runs)                  # read (readers/segments.py)
    fig = new_figure()                                                # style.py: the look, the info strip
    left, right = fig.subplots(1, 2)
    draw_norm(left, averages.times_au, averages.norm, averages.rest_norm, starts_au)   # one line per panel
    draw_mean_positions(right, averages.times_au, averages.z_mean_bohr, averages.r_mean_bohr, starts_au)
    add_info_box(fig, info, ["averages.dat", "run_info.json"])        # the files it read
    return [("averages", fig)]                                        # (file stem, figure)
```

and in `PLOTS`: `("averages", figure_averages, ["averages.dat"])`, the files it needs in the run directory.

## Tests

pytest, one ctest per test module, registered by a glob in the repository's `tests/CMakeLists.txt`: the file
name decides how. `tests/test_runs_<name>.py` becomes `PlotRuns.<name>`: it reads the fixture runs, which the
C++ test `PlotRuns.Fixtures` ([tests/test_plot_fixtures.cpp](../../tests/test_plot_fixtures.cpp)) writes with
the current program into `build/tests/plot_fixtures` (about 160 s; ctest runs it first), and gets
`--fixture-runs`. Any other `tests/test_<name>.py` becomes `Plots.<name>` and runs without them.
`tests/conftest.py` also makes runs that stopped early from one fixture. The `PlotRuns.*` tests are not in the
quick tier: after changing the tool, run them too ([docs/CONTRIBUTING.md](../../docs/CONTRIBUTING.md#testing)):

```sh
OMP_NUM_THREADS=8 ctest --test-dir build -R '^(PlotRuns|Plots)\.'
uv run python -I -B -m pytest -c tools/plots/pytest.ini --rootdir=tools/plots tools/plots \
    --fixture-runs=build/tests/plot_fixtures
```

The fixtures are small and not fully converged; `tests/test_runs_diagnostics.py` asserts the findings each
one raises. The S matrix rebuilt from `bk_coefficients.bin` with `copies/s_matrix.py` must equal
`s_matrix.tsv` to 1e-12 of the largest |S| at each energy (`tests/test_runs_s_matrix.py`): the check that the
BK reader, the time sums, the prefactor and the phases are right.

## Limits

- One run at a time (its restart segments joined); no comparison of runs.
- Only what the run's files hold: the wave packet itself is not plotted yet.
- The results are the on-the-fly analysis (`scattering_summary.tsv`); the fine analysis (`fine_analysis/`,
  also of an `onlyAnalysis` run) is not plotted yet.
