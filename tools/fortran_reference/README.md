# Fortran reference tools

The scripts here run the Fortran SPO-DVR code to produce the reference data of the regression
tests (`tests/reference/`, `tests/crplib/reference/`; the gates are described in
[docs/validation.md](../../docs/validation.md)). They write straight into those directories: run
them in a scratch clone or worktree when you only want to experiment.

| Script | What it does |
|---|---|
| [reference_config.sh](reference_config.sh) | The one home of the pins, paths, thread count and environment, and the functions the scripts share (sourced, not run) |
| [build_fortran_reference.sh](build_fortran_reference.sh) | Builds the pinned upstream commit into a read-only reference install |
| `run_setup_reference.sh [case]` | Potential-only run: `fortran_setup_stdout.txt`, `debug_asyp.dat`, `pes_reference.txt` (gates 1, 2) |
| `run_propagation_reference.sh [case]` | Full run: `fortran_stdout.txt`, `fortran_averages.out`, `fortran_fort21.bin`, `fortran_fort25.txt` (gates 3, 4) |
| `run_restart_reference.sh` | hardwall_small as two runs, the second a Fortran restart: `fortran_restart_stdout.txt` |
| `run_fine_analysis_reference.sh [case] [variant ...]` | Analysis-only runs on the committed `fortran_fort21.bin`: `fine_<variant>/` (gate 5) |
| `run_crp_reference.sh [--generate] kernels\|<case>` | The CRP library driver ([crp_driver/](crp_driver/)) on the committed point lists and data sets of `tests/crplib/reference/`; `--generate` first rewrites them (needs the C++ build) |
| [regenerate_all.sh](regenerate_all.sh) | Regenerates every reference file, in dependency order ([below](#regenerating-all-references)) |
| [refdata.py](refdata.py) `manifest` | Writes (`record`) and deletes (`delete`) the lines of `MANIFEST.tsv` ([below](#manifesttsv)) |
| [ReferenceManifest.cmake](ReferenceManifest.cmake) | Checks every file a `MANIFEST.tsv` lists against its sha256 (`ctest -R ReferenceManifest`) |
| `test_reference_tools.sh` | The negative tests of these tools (`ctest -R FortranReference`; no Fortran run) |

Each script's header gives its arguments, outputs, environment and run time. Every run script
records the files it writes in the case's `MANIFEST.tsv`.

## Environment variables

| Variable | Meaning |
|---|---|
| `SPODVR_FORTRAN_REF_DIR` | The reference install that the run scripts use and that `build_fortran_reference.sh` installs into. Default: `~/work/SPO-DVR_ref/<id>`, with the id pinned in `reference_config.sh` (`FORTRAN_REF_ID` for the runs, `FORTRAN_BUILD_ID` for the build) |
| `SPODVR_FORTRAN_MIRROR` | The upstream git mirror the build archives the sources from. Default: `~/work/SPO-DVR` (read only: the build never writes there) |
| `SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT` | `1`: run although an installed run-time package differs from the install's `BUILDINFO.tsv` (a warning instead of a refusal) |

The variables that named the executable or the object directory directly are retired: the scripts
refuse to run while one is set. Point `SPODVR_FORTRAN_REF_DIR` at an install instead.

## The reference install

An install is a directory `~/work/SPO-DVR_ref/<id>/` with a fixed layout:

| Path | Content |
|---|---|
| `bin/spodvr.x` | The executable |
| `OBJS/` | Its objects and modules (`*.o`, `*.mod`), linked into the CRP driver |
| `BUILDINFO.tsv` | Provenance (below) |
| `SHA256SUMS` | sha256 of the files above (`sha256sum -c` in the install) |
| `build.log`, `sources.tar` | The DoBuild log and the exact sources (installs made by the build script) |

`v5.19` is the 5.19 archive (`~/work/SPO-DVR_v5.19`, symlinked; it must never be rebuilt or
deleted). Installs are read-only; to rebuild one: `chmod -R u+w <install> && rm -r <install>`.

`BUILDINFO.tsv` is a TSV file: `#` comment lines, the header `key<TAB>value`, rows sorted by key,
`-` for a value that was not recorded. It is parsed, never sourced. A copy of each install's file
is committed in [buildinfo/](buildinfo/).

| Key | Value |
|---|---|
| `id`, `version`, `commit`, `svn`, `tree` | The install id, program version, upstream commit, svn revision and git tree hash |
| `source`, `tarball.sha256` | Where the sources came from, and the sha256 of `sources.tar` |
| `defines`, `flags` | The preprocessor defines and compiler flags of the build |
| `exe.sha256` | sha256 of `bin/spodvr.x` |
| `package.build.<name>` | dpkg versions of `cpp-15`, `gfortran-15`, `libc6-dev`, `libfftw3-dev` |
| `package.runtime.<name>` | dpkg versions of `libc6`, `libfftw3-double3`, `libgfortran5`, `libgomp1` |
| `zgvb` | The SSE libmvec functions the executable imports (`_ZGVb*`, sorted) |
| `date` | The build date (UTC) |

## What every run checks

- **The install** (`check_install`): every file against `SHA256SUMS`, the executable against
  `exe.sha256`, and the installed run-time packages against `package.runtime.*`: a libm,
  libgfortran, libgomp or FFTW upgrade can move results by an ulp.
- **The environment** (`run_clean`): `env -i` with a fixed `PATH`, `LC_ALL=C` and
  `OMP_NUM_THREADS=4` (`FORTRAN_REF_THREADS`), whatever the caller's environment holds.
- **The scratch directory** (`new_scratch_dir`): a new `mktemp -d` directory, or a named one that
  must not exist yet (a stale `.pes` or `pesmorphparameters.dat` would be read).
- **The run** (`run_fortran`): exit status 0 and every required output file present (a Fortran
  `STOP` on bad input exits with status 0).
- **The CRP driver build** (`build_fortran_driver`): the install as above, and the installed
  gfortran equal to `package.build.gfortran-15` (its `.mod` files are compiler-version specific).

Every refusal prints `error: <what>; <how to fix>` and exits with status 1.

## Building a reference install

```bash
tools/fortran_reference/build_fortran_reference.sh               # about 40 s (measured 2026-11-12)
tools/fortran_reference/build_fortran_reference.sh --check-log LOG
tools/fortran_reference/build_fortran_reference.sh --check-binary FILE...
```

The build checks the pinned commit and tree hash in the mirror, builds from `git archive` of
`sources.dir` with `GNU.mk`'s `-Ofast -mtune=native -march=native` replaced by `-O3` (the upstream
default enables fast-math and FMA), and refuses:
- a log with an error line, or with any line outside the allow-lists: PreProcess, compile and
  link lines with exactly the pinned defines and flags, and DoBuild's progress messages;
- an executable or object with VEX/EVEX instructions, AVX libmvec imports (`_ZGVc`, `_ZGVd`,
  `_ZGVe`) or crtfastmath (`set_fast_math`).

It then installs with one `mv`, makes the install read-only and writes `buildinfo/<id>.tsv`, which
is committed with the reference data made from it. `--check-log` and `--check-binary` run only
those checks, e.g. on a build made elsewhere.

## MANIFEST.tsv

Each reference case directory has a `MANIFEST.tsv`: one line per file a Fortran run generated
there (not the inputs: the case's `DVR.inp` and `input.json`, the CRP `dataset/` and `*.in`).
The run scripts write it through `refdata.py manifest record`, which replaces the lines of the
files it records and refuses the inputs (and `fortran_unformatted/`, hand-made test records), so
that a regeneration can never delete one. It is a TSV file: `#` writer lines, a header row, `-` for a value that was not
recorded, rows sorted by path. Hashes other than the file's own are shortened to 12 hex digits.

| Column | Value |
|---|---|
| `path` | The file, relative to the case directory |
| `sha256` | Its sha256 |
| `script` | The script that generated it |
| `commit` | `git describe --always` of this repository, with `-dirty` when a tracked file outside the reference trees (`tests/**/reference/`) differs from it: those change in every regeneration |
| `fortran_ref` | The install id (`BUILDINFO.tsv`) and `<program>:<sha256>` of the program that ran: `spodvr.x`, or the CRP driver linked against the install's objects |
| `input_sha` | `<input>:<sha256>` of each input as run, comma-separated: the `DVR.inp` with the script's changes, a staged CRP data set (`dataset`, hashed as a tree: the sha256 of its sorted `sha256sum` lines), the fine analysis's `fort.22` (equal to the case's `fortran_fort21.bin`), the driver's `*.in` |
| `runtime_drift` | The run-time packages whose installed version differs from `BUILDINFO.tsv` (`package=version`), when `SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1` let the run go ahead |
| `cpu` | The CPU model and `flags:<sha256>` of its feature flags (libm picks code by CPU feature) |
| `date` | The UTC date of the run |

The lines of the files generated before the MANIFESTs existed carry `fortran_ref`
`5.19 (pre-MANIFEST)`, `-` in the columns not recorded then, and the script that made them
(`refdata.py manifest record --pre-manifest`).

`ctest -R ReferenceManifest` checks that every case directory has a `MANIFEST.tsv` and that every
listed file exists with its recorded sha256: a missing or edited reference, or a case left without
MANIFEST by an interrupted regeneration, fails there, where the regression tests would skip
(`FilesMatch`; `ChangedOrMissingFileFails` shows that it catches one changed byte, one deleted file
and a deleted MANIFEST).

## Regenerating all references

```bash
tools/fortran_reference/regenerate_all.sh      # about 45 min on 4 threads (estimate, below)
git status --short tests/reference tests/crplib/reference
```

The script refuses unless the reference trees are clean in git, checks the install, then deletes
every file the MANIFESTs list and the MANIFESTs themselves, so that nothing stale survives. It
runs the steps in dependency order: setup, propagation, the restart, fine analysis (on the
regenerated `fortran_fort21.bin`) and the CRP driver on the committed point lists. Each step
checks its required outputs; the end runs the `ReferenceManifest` check. Its list of steps is the
list of reference cases: a new case adds its step in the commit that adds the case. The scratch
directories of all runs stay in one new `mktemp -d` directory, printed at the start. The run time
is estimated from single-case Fortran runs measured on 2026-11-11 (cu111_dailight_small, 7000
steps: 1212 s with 5.19, 1570 s with 5.24 `-O3`, 4 threads); the whole run is not yet timed.
