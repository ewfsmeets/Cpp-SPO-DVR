# Fortran reference tools

The scripts here run the Fortran SPO-DVR code to produce the reference data of the regression
tests (`tests/reference/`, `tests/crplib/reference/`, `tests/scm/reference/`; the gates are described in
[docs/validation.md](../../docs/validation.md)). They write straight into those directories: run
them in a scratch clone or worktree when you only want to experiment. Which Fortran the reference
is, and why it is built with `-O3`: [docs/validation.md](../../docs/validation.md#reference-data).
The recipe, in order: [build](#building-a-reference-install) an install,
[regenerate](#regenerating-all-references) the references, [compare](#comparing-two-reference-trees)
them with the committed ones; [moving to a newer upstream commit](#moving-the-reference-to-a-newer-upstream-commit)
strings the three together.

| Script | What it does |
|---|---|
| [reference_config.sh](reference_config.sh) | The one home of the pins, paths, thread count and environment, and the functions the scripts share (sourced, not run) |
| [build_fortran_reference.sh](build_fortran_reference.sh) | Builds the pinned upstream commit into a read-only reference install |
| `run_setup_reference.sh [case]` | Potential-only run: `fortran_setup_stdout.txt`, `debug_asyp.dat`, `pes_reference.txt` (gates 1, 2) |
| `run_propagation_reference.sh [case]` | Full run: `fortran_stdout.txt`, `fortran_averages.out`, `fortran_fort21.bin`, `fortran_fort25.txt` (gates 3, 4) |
| `run_restart_reference.sh` | hardwall_small as two runs, the second a Fortran restart: `fortran_restart_stdout.txt` |
| `run_fine_analysis_reference.sh [case] [variant ...]` | Analysis-only runs on the committed `fortran_fort21.bin`: `fine_<variant>/` (gate 5) |
| `run_crp_reference.sh [--generate] kernels\|<case>` | The CRP library driver ([crp_driver/](crp_driver/)) on the committed point lists and data sets of `tests/crplib/reference/`; `--generate` first rewrites them (needs the C++ build) |
| `run_scm_reference.sh [--generate] example\|synthetic` | The SCM driver ([scm_driver/](scm_driver/)) on the committed inputs of `tests/scm/reference/` ([below](#the-scm-driver)); `--generate` first rewrites them |
| [regenerate_all.sh](regenerate_all.sh) | Regenerates every reference file, in dependency order ([below](#regenerating-all-references)) |
| [refdata.py](refdata.py) `manifest` | Writes (`record`) and deletes (`delete`) the lines of `MANIFEST.tsv` ([below](#manifesttsv)) |
| [refdata.py](refdata.py) `compare` | Compares two reference trees against an expected-differences file such as [expected/r1713.tsv](expected/r1713.tsv) ([below](#comparing-two-reference-trees)) |
| [ReferenceManifest.cmake](ReferenceManifest.cmake) | Checks every file a `MANIFEST.tsv` lists against its sha256 (`ctest -R ReferenceManifest`) |
| `test_reference_tools.sh`, `test_refdata.py` | The negative tests of these tools (`ctest -R FortranReference`; no Fortran run) |

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
- **The CRP and SCM driver builds** (`build_fortran_driver`): the install as above, and the
  installed gfortran equal to `package.build.gfortran-15` (its `.mod` files are compiler-version
  specific).

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
there (not the inputs: the case's `DVR.inp` and `input.json`, the CRP `dataset/`, the driver
inputs `*.in`).
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

`refdata.py manifest record --pre-manifest` records files generated before the MANIFESTs existed
(it bootstrapped the MANIFESTs of the 5.19 data): `fortran_ref` `5.19 (pre-MANIFEST)`, `-` in the
columns not recorded then, and the script that made them. Since the data were regenerated with
5.24 r1713, every line names that install.

When a run-time package of the install is upgraded and the rebuilt install is byte-identical
(as for libc6 2.43-2ubuntu2.7 on 2026-11-14), its `BUILDINFO.tsv` and the tracked copy in
[buildinfo/](buildinfo/) are re-recorded with the new version, so new reference data need no
drift flag. The run-time packages of an older MANIFEST line are then those of the copy at the
commit its `commit` column names (git history).

`ctest -R ReferenceManifest` checks that every case directory has a `MANIFEST.tsv` and that every
listed file exists with its recorded sha256: a missing or edited reference, or a case left without
MANIFEST by an interrupted regeneration, fails there, where the regression tests would skip
(`FilesMatch`; `ChangedOrMissingFileFails` shows that it catches one changed byte, one deleted file
and a deleted MANIFEST).

## The SCM driver

[scm_driver/scm_reference_driver.f90](scm_driver/scm_reference_driver.f90) calls the SCM
distorted-surface correction of the reference build in SPO-DVR's order for SystemID
`CRPDistrtSurface` (`Set_PES_DataSet`, `InitialiseDistPESModule`, `GetCRPDistSurfParams`,
`InitializeCRPDistSurfModule`, `SetupDistPESModule`), one case per process (the Fortran sets its
module state up only once), and
prints per point the Cartesian centre of mass, the correction ΔV alone (Hartree), the CRP at the
ideal-lattice position and SPO-DVR's total (eV), plus the parameters and slab atoms the Fortran
read; every double as its 16 hex digits and in decimal. It stops with status 1 unless the total
equals CRP + ΔV / `Constants_eV2Ha` bitwise at every point (then its x, y and u/s, v/s are the
module's). [scm_driver/build.sh](scm_driver/build.sh) links it against the install's SCM and CRP
objects; [scm_driver/make_scm_inputs.py](scm_driver/make_scm_inputs.py) writes the inputs
(`--generate`). All runs use the CRP data set `potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v`.

| Case | Inputs | Runs (`<run>.out`) |
|---|---|---|
| `example` | `points.in` (384 points); `params.dat` and `examplePosfile.dat` of `potentials/H2Cu111SRP_C3v.dir/` | `rescaling_on`, `rescaling_off`: configuration 1, cut-off 16 bohr |
| `synthetic` | `points.in` (284 points); `params.in` (layer 2 unlike layer 1); `slab_v3.in`, `slab_v2.in` (three configurations each, random displacements; three version-3 atoms at the cut-off: 3 ulp of c² inside, on it, 1 ulp outside) | `v3_slab3_rescaling_on`: configuration 3, cut-off 9.5 bohr; `v2_slab2_rescaling_off`: configuration 2, cut-off 4 bohr (the Fortran applies none in version 2; a reader that does would keep 1 atom of layers 1–2 instead of 10) |

The points cover the four sites from Z = −1 bohr (below the CRP's Zmin, where only the correction
sees the unclamped Z) to 6 bohr in four orientations, an r sweep through P17 and P18, H atoms
1e-3 and 1e-9 bohr either side of p7 and p8 (the switch), the V_asy points
(0, 0, zasymp, r, π/2, 0), random points, and 20 random points from Z = 8 bohr up to the plane
Z = zasymp (6 on it). `ctest -R Scm.ReferenceFiles` checks that every run
is complete and that the C++ readers read the parameter and slab files bit for bit as the Fortran
did. `ReferenceManifest.FilesMatch` hashes the files the SCM MANIFESTs list; its check for a case
without MANIFEST, and `refdata.py compare`, cover only `tests/reference/` and
`tests/crplib/reference/` so far: after a regeneration, compare the SCM files with `git diff` (a
rerun is byte-identical).

## Regenerating all references

```bash
tools/fortran_reference/regenerate_all.sh      # about 45 min on 4 threads (see below)
git status --short tests/reference tests/crplib/reference tests/scm/reference
```

The script refuses unless the reference trees are clean in git, checks the install, then deletes
every file the MANIFESTs list and the MANIFESTs themselves, so that nothing stale survives. It
runs the steps in dependency order: setup, propagation, the restart, fine analysis (on the
regenerated `fortran_fort21.bin`), the CRP driver and the SCM driver on the committed inputs. Each step
checks its required outputs; the end runs the `ReferenceManifest` check. Its list of steps is the
list of reference cases: a new case adds its step in the commit that adds the case. The scratch
directories of all runs stay in one new `mktemp -d` directory, printed at the start. The whole
run takes about 45 min with the r1713 install on 4 threads: 2381–2538 s (2026-11-12) before the
case cu111_dailight_short_lowz was added, whose run takes 158 s (2026-11-12); the
longest case, cu111_dailight_small with 7000 steps, took 1570 s when run alone (measured
2026-11-11).
Then compare the regenerated tree with the committed one ([below](#comparing-two-reference-trees)).

## Comparing two reference trees

```bash
python3 -I -B tools/fortran_reference/refdata.py compare OLD NEW --expected tools/fortran_reference/expected/r1713.tsv
python3 -I -B tools/fortran_reference/refdata.py compare OLD NEW --expected FILE --no-provenance
```

`compare` reads the reference files (`tests/reference/`, `tests/crplib/reference/`) of two
repository checkouts: `OLD` with the committed references (e.g. the main checkout), `NEW` with
regenerated ones (e.g. a worktree after `regenerate_all.sh`). It classifies each file by the
row of the expected-differences file whose glob matches its path (a file that two rows match is
refused); a file that no row matches must be byte-identical, and a file added in `NEW` always
fails. It prints one line per file that is not identical (`ok` or `FAIL`, the class, what
differs: lines, numbers, the largest absolute difference and printed units, a changed line count,
or for BK files the fraction of the gate-4 scale), a note for each row that matches no file, and
a summary. Exit status 0 only if every difference is allowed; it takes 4–7 s on the full trees
(measured 2026-11-12).

It also checks the provenance of `NEW` (`--no-provenance` skips only this, e.g. for files copied
in by hand): every case has a `MANIFEST.tsv` that lists exactly its generated files with their
sha256; all lines name the install of the expected file's `# install: <id>` line, whose
`BUILDINFO.tsv` copy in [buildinfo/](buildinfo/) has the executable the lines name
(`exe.sha256`); no line records run-time drift; and every fine-analysis line names its input
`fort.22`, which is the case's `fortran_fort21.bin`.

The expected-differences file is a TSV file: `#` writer lines, among them exactly one
`# install: <id>` (the install whose output it describes), the header
`glob<TAB>allowed<TAB>limit<TAB>reason`, `-` for an empty field, rows sorted by glob. The globs
must be disjoint: in a glob, `*` and `?` match within one path component. The classes, a closed
list:

| Class | Allowed differences | `limit` |
|---|---|---|
| `identical` | None: the same as no row, written to document a measurement | `-` |
| `volatile` | Text lines that differ between two runs of one build: the version line (build and run dates), the thread count, the set-up timings, the `iz = N` lines (printed from a parallel loop), the CPU timing table at the end | `-` |
| `labels` | As `volatile`, plus the label whitelist `LABELS_519_TO_524` in refdata.py: the stdout lines whose text changed from 5.19 to 5.24 while their numbers did not. Every `labels` row uses it, so a later bump of the reference must review it | `-` |
| `abs` | As `labels`, plus every number within `limit` of its old value | absolute limit |
| `printed_units` | As `labels`, plus every number within `limit` units of its last printed digit | units |
| `signs_only` | As `labels`, plus a sign change of the numbers at the listed positions (1-based) of a table row: a line of numbers only, with at least as many numbers as the largest position (an S-matrix row of `fort.25`, not a block header). With the sign, their printed digits may change by one unit (the negated value is computed another way) | positions, e.g. `4,6` |
| `bk_scaled` | BK coefficients (`fortran_fort21.bin`): identical header records, every coefficient within `limit` times the gate-4 scale of `OLD` (its largest coefficient modulus, at least 1e-4) | fraction of the scale |
| `unread` | Any change of content (no test reads the file; the line still reports what changed, e.g. the line count of an emptied file); it must exist in both trees | `-` |
| `removed` | The file may be missing from `NEW` (if present, it must be identical) | `-` |

Text files are compared line by line, after the volatile lines are dropped and (for every class
but `volatile`) the whitelisted label lines are put in one form. A number may move inside its
right-aligned field (a minus sign takes a blank), but the field width, the number of numbers and
every other character must stay. A value that is not finite (NaN or out of the double range, a
printed number or a BK coefficient) never passes unless it is bitwise unchanged, and a BK file of
`OLD` with one has no gate-4 scale. In the CRP driver outputs the hex word of each value is skipped: the 17-digit decimal
on the same line determines the double.

[expected/r1713.tsv](expected/r1713.tsv) lists the differences of the switch from the 5.19 data
to the 5.24 r1713 `-O3` build, measured on 2026-11-11, each with its reason; the regenerated tree
of the switch passed it. `ctest -R FortranReference.Compare` checks that it matches the committed
files, that every class it allows
passes, and that each fault fails: an edited number, a deleted file, an unknown label, a change
of |S| or of a block header's sign in a fine `fort.25`, a BK change beyond the limit, NaN BK
coefficients, a printed number out of range, a changed number count, three MANIFEST faults, an
install other than the expected file's, a fine-analysis line without `fort.22`, overlapping rows
and a malformed expected file ([test_refdata.py](test_refdata.py)).

## Moving the reference to a newer upstream commit

The pin is in [reference_config.sh](reference_config.sh): `FORTRAN_BUILD_*` names the upstream
commit the build script builds, `FORTRAN_REF_ID` the install every run uses. The commit the port is
up to speed with, and the port status of each upstream change, are in
[docs/upstream_sync.md](../../docs/upstream_sync.md). Steps:

1. **Review the upstream changes.** After `git -C ~/work/SPO-DVR pull`, list the commits after the
   pin with `git -C ~/work/SPO-DVR log --oneline <pin>..master` and read
   `git -C ~/work/SPO-DVR diff -w <pin> <new> -- sources.dir/SOURCES.dir`. Give each change a row
   in the table of docs/upstream_sync.md (ported, a deviation, or not ported).
2. **Pin and build.** Set `FORTRAN_BUILD_COMMIT`, `FORTRAN_BUILD_TREE`
   (`git -C ~/work/SPO-DVR rev-parse <new>^{tree}`), `FORTRAN_BUILD_SVN` and
   `FORTRAN_BUILD_VERSION`, and check that `FORTRAN_BUILD_GNU_MK_LINE` is still line 19 of
   `sources.dir/GNU.mk`. Run `build_fortran_reference.sh` ([above](#building-a-reference-install));
   a log line outside its allow-lists stops the build: understand it before widening a list.
3. **Regenerate** in a scratch worktree with `SPODVR_FORTRAN_REF_DIR` set to the new install:
   `regenerate_all.sh` ([above](#regenerating-all-references)).
4. **Compare** with the committed references. Write `expected/<id>.tsv` (with
   `# install: <id>`): a first run with only the header row lists every difference as a `FAIL`
   line. Give each difference a row with its reason, or stop when one is not understood; review
   the label whitelist of refdata.py when stdout labels changed. `refdata.py compare` must then
   exit with status 0 ([above](#comparing-two-reference-trees)).
5. **Switch and test.** Set `FORTRAN_REF_ID` to the new id; run the full ctest, the long tests
   (`SPODVR_LONG_TESTS=1`) and the CUDA tests: every gate passes, with the same skipped tests as
   before. A failing gate is a finding, not a tolerance to raise (docs/validation.md "Tolerances").
6. **Commit together:** the regenerated data with their MANIFESTs, `buildinfo/<id>.tsv`,
   `expected/<id>.tsv`, the pin, docs/upstream_sync.md (the commit, its svn revision and date, and
   the new rows; `tools/check_docs.py` fails while its commit differs from `FORTRAN_BUILD_COMMIT`
   or `FORTRAN_REF_ID` from `FORTRAN_BUILD_ID`) and the "Reference history" of [docs/validation.md](../../docs/validation.md#reference-history).
