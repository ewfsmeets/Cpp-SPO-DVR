# Upstream sync

The Fortran SPO-DVR code continues upstream. This page names the upstream commit the C++ port is
up to speed with: every upstream change up to that commit has been reviewed, its port status is
in the table below, and the Fortran reference data of the regression tests come from it
([validation.md](validation.md#reference-data)). The port started from version 5.19.

| Item | Value |
|---|---|
| Upstream repository | the git-svn mirror of the SPO-DVR SVN trunk (Leiden), branch `master`; locally `~/work/SPO-DVR`, updated with `git pull` and never built in |
| Upstream commit | `dce0ecc44a401098528e9ef94c54e1086d393ebb` |
| svn revision | r1713 |
| Date | 2026-11-09 |
| Version | 5.24 plus the fine-analysis sign fix (`s1_real` → `s1_imag` in `Analyzor.f90 smatrix_calc`) |
| Reference install | `r1713-dce0ecc`, built with `-O3` ([tools/fortran_reference/README.md](../tools/fortran_reference/README.md)) |

The pin itself lives in
[tools/fortran_reference/reference_config.sh](../tools/fortran_reference/reference_config.sh)
(`FORTRAN_BUILD_COMMIT`, `FORTRAN_BUILD_SVN`); `tools/check_docs.py` (check `upstream-pin`) fails
when this page names another commit or svn revision, or when the runs use another install than
the pinned build (`FORTRAN_REF_ID` ≠ `FORTRAN_BUILD_ID`).

## Newer upstream commits

```bash
git -C ~/work/SPO-DVR pull                                               # update the mirror
git -C ~/work/SPO-DVR log --oneline dce0ecc..master                      # commits after the pin
git -C ~/work/SPO-DVR diff -w dce0ecc master -- sources.dir/SOURCES.dir  # their source changes
```

On 2026-11-12 the mirror's `master` was `dce0ecc` itself: no newer commit. Each newer change gets a
row below; moving the reference to a newer commit follows the recipe in
[tools/fortran_reference/README.md](../tools/fortran_reference/README.md#moving-the-reference-to-a-newer-upstream-commit).

## Port status of the upstream changes since 5.19

5.19 is git `601d833` (svn r1572, 2023-07-13). Up to `dce0ecc`, 43 commits changed
`sources.dir/SOURCES.dir` (`git -C ~/work/SPO-DVR log --oneline 601d833..dce0ecc -- sources.dir/SOURCES.dir`);
many source files changed only by the GPL 3 license header (f52dbf3, 2025-09-19). For inputs that
use only 5.19 features, the operators, absorber formulas, input defaults and file formats are
unchanged: built with `-O3`, 5.24 reproduces the 5.19 results bitwise on every gate case (checked
2026-11-11), except the signs that the `s1_imag` fix changes and the stdout labels. The status is
one of:
- **ported**: the C++ does what the reference does;
- **deviation**: the C++ does it differently on purpose; the ID is defined in
  [fortran_deviations.md](fortran_deviations.md);
- **not ported**: the last column says what such an input does in the C++ today, or why the
  change does not concern the C++ ("not relevant").

The last two rows are 5.19 features the port lacks; they are listed because later work ports them.

| Upstream change | Commit, date | Status | In the C++ today |
|---|---|---|---|
| `RStart = 0` rejected (`InputNew.f90 ValidateInput`: `.RLE.` instead of `<`) | 9761502, 2025-11-25 | **ported** | `RStart ≤ 0` is an input error (exit code 2), with the 1e-12 tolerance of `.RLE.` (`src/input/validate.cpp`) |
| Fine-analysis interpolation of Im S: `s1_real` → `s1_imag` (`Analyzor.f90 smatrix_calc`) | dce0ecc, 2026-11-09 | **ported** | fixed in the C++ before upstream, together with the two other interpolation bugs, which the reference keeps ([fortran_fine_analysis_bugs.md](fortran_fine_analysis_bugs.md)) |
| Optical potentials at the low end of the Z grid, scattering and specular (keys `ZOptLow…`, `ZSpOptLow…`) | c441b95 and fixes, 2025-09-12 to 09-14 | **ported** (its input checks: **deviations** low-Z D1–D4) | on the CPU and the GPU; gates 3 and 4 compare it with the reference in the case `cu111_dailight_short_lowz` ([validation.md](validation.md)); the Fortran's silent cases are input errors (low-Z D1–D3) |
| PES morphing (`Potential.f90 MorphPotential`): no input key, runs when `pesmorphparameters.dat` is in the run directory | 27d57d6, 2024-02-01 (tool acc06c9) | **not ported** | ignored without a warning: the C++ never opens `pesmorphparameters.dat`, so V is the unmorphed potential |
| Distorted-surface (SCM) CRP potential, SystemID `CRPDistrtSurface` (keys `DistortedPESMode`, `SCMParFilepath`, `PosFilepath`, `PosFileSlabNum`, `EnableSurfExpan`, `SCMCutoff`) | 41e39e5, 2023-09-05 (example files c472cfe) | **not ported** | input errors (exit code 2): the SCM keys are unknown keys, the SystemID an unknown system ID. Open physics questions: [scm_open_questions.md](scm_open_questions.md) |
| PESes `H2W110` (H₂/W(110), vdW CRP) and `HW110` (H/W(110)) | 947d613, 2024-06-24; 8ed9de3, 2025-09-17; fixes to 2025-11-28 | **not ported** | unknown system ID (input error, exit code 2); the data files are in `potentials/H2W110.dir` |
| Build: `GNU.mk` `O3FLAGS = -Ofast -mtune=native -march=native` instead of `-O3` (fast-math and FMA) | 3b3a8c0, 2026-02-05 | **not ported** | not relevant: the C++ uses strict floating point without FMA; the reference install is built with `-O3` ([validation.md](validation.md#reference-data)) |
| Averages: the norm, ⟨Z⟩ and ⟨r⟩ reductions with `SCHEDULE(STATIC,1)` (`Averages.f90`) | 8c6d5b6, 8abe04a, 2024-02 | **not ported** | not relevant: the C++ sums in its own fixed order; the 5.24 averages of the five propagation cases equal the 5.19 ones ([validation.md](validation.md#reference-data)) |
| OpenMP `DYNAMIC` → `STATIC,1` in `SplitOperator`, `KineticOp`, `WaveFunction`, `Potential`, `FluxAnalyzor` | 7ec514c, 2024-02-07 | **not ported** | not relevant: independent loops, results unchanged |
| Specular Z absorber loop over all Z points (`Absorption.f90`) | be75c78, 2025-09-13 | **ported** | the specular Z factors are applied at every Z point on the CPU and the GPU; for inputs without a low-Z absorber the results are unchanged (the factor is exactly 1 outside the absorbers) |
| FFTs through oneMKL (Intel compiler builds) | a614a65, 2023-09-06 | **not ported** | not relevant: FFTW builds are unchanged, and the C++ has its own FFT layer (oneMKL's FFTW3 interface or FFTW) |
| FFT wrappers: typed dummies, the in-place r2c/c2r and real-to-real paths of `fftw_2_scs.f90` now stop with an error, an allocation fix of the strided real FFTs in `ffts.f90` | 4834d19, 2025-09-21; 3862347, acd72f0, 2023-12 | **not ported** | not relevant: the C++ has its own FFT layer, and 5.24 built with `-O3` reproduces the 5.19 results bitwise |
| PES4Cu100 (`H2Cu100PES4`) and the NAG library removed; fixed-format compiles dropped; make files renamed (`GFORTRAN.mk` → `GNU.mk`) | 6b5c2cf, 2023-09-06; 11ed220, 5fe2c9d | **not ported** | not relevant: never ported (`H2Cu100PES4` is an unknown system ID) |
| Stdout labels and columns (absorber lines, lattice "U/V" labels, grid start and end indices, the NAG and MKL lines) | 293f410 and others, 2025-09 | **not ported** | not relevant: the gates read none of the changed lines; the reference comparison maps them (`LABELS_519_TO_524` in `tools/fortran_reference/refdata.py`) |
| Dead code removed (`MyFunctions.f90 LengthOfVectorSquared`, `TransToCoordSystem`) | ab47c48, 2026-10-25 | **not ported** | not relevant |
| Flux analysis (`FluxAnalyzor.f90`: reaction probability at r_c, its (x, y) density, flux restart files); already in 5.19 | — | **not ported** (its input checks: **deviations** flux Q1, flux Q2) | a propagating run (new or restart) with `doFluxAnalysis: true` stops with exit code 3 once the input checks have passed; an analysis-only run prints a note (flux Q6), a potential-only run ignores the flux keys |
| Other 5.19 features: the `.pes` potential cache, write2DPES, ProjOp 1, MPI, the data-file PESes | — | **not ported** | the potential is always computed (a `.pes` file is ignored); write2DPES and ProjOp 1 are rejected ([fortran_deviations.md](fortran_deviations.md)); OpenMP only; the other PESes are unknown system IDs (their data files are in `potentials/`) |

<details>
<summary>Details of the new upstream features (keys, defaults, quirks)</summary>

- **Low-Z optical potentials** (`InputNew.f90 ReadInput`, `Setup.f90`, `Absorption.f90`).
  - Keys: `ZOptLowStart`, `ZOptLowEnd`, and `ZOptLowEnergy` or `ZOptLowFactor`; the same four as
    `ZSpOptLow…` for the specular grid. They are not pre-registered, so a `DVR.inp` needs
    `Name : REAL : value;`.
  - Switched on only when Start and End are both given with Start ≤ End; then Start ≥ ZStart,
    End ≤ ZAna, Energy or Factor, and Energy ≥ 0 are required.
  - At the same Z the low-end value overwrites a high-end absorber value instead of adding to it.
- **PES morphing** (`Potential.f90 MorphPotential`, examples in upstream `tools.dir/morphing`).
  - It scales or shifts FBR coefficients (J, mJ, n, m) with Gaussians in Z and r, after V is
    computed, and needs MjMax = JMax.
  - From reading the code: it runs before the `.pes` cache is written, the cache name does not
    record it, and its loops assume the full (x, y) grid (`UseSymmetry`).
- **SCM** (`DistortedSurfacePES.dir`, SystemID `CRPDistrtSurface`): a CRP data set plus a
  parameter file and a file of surface atom positions. The keys are read only for that SystemID.
  In the CRP library it adds only an accessor for the interlayer spacings (`CRP_QD_Wrapper.f90`).
- **H2W110, HW110** (`H2W110PES-vdW.f90`, `HonW110PES-vdW.f90`): C2v; they read their data files
  from `potentials/H2W110.dir`. HW110 builds on the H₂ PES.

</details>
