# Deviations from the Fortran code

Every known difference between the C++ program and the Fortran reference (upstream 5.24 at svn
r1713, built with `-O3`: [validation.md](validation.md#reference-data)): Fortran bugs that the C++
fixes, Fortran quirks that it reproduces on purpose, features that were not ported, and places
where the Fortran computes in quad precision while the C++ uses double. The port started from 5.19,
which gives the same results as the reference for its features; where the two versions differ, the
entry says so. Add every new deviation here. The fine-analysis bugs are described in detail in
[fortran_fine_analysis_bugs.md](fortran_fine_analysis_bugs.md); the gates that bound the numerical
differences are in [validation.md](validation.md).

- `nZSpecular` read depends on `nproj` before `ProjOp` is read (always uses default 2).
- flux Q1: the Fortran stops when `EMaxFluxAna .RLT. EMaxWP` (message `'EMaxFluxAna < EMinWP'`):
  it rejects every flux energy range that ends inside the wave packet and accepts ranges beyond
  it. C++: error if EMaxFluxAna > EMaxWP (tolerant, like the EMinFluxAna check). Details and a
  demonstration: [fortran_input_bugs.md](fortran_input_bugs.md#flux-q1-emaxfluxana-is-checked-the-wrong-way-round).
- flux Q2: the Fortran rejects only `deltaEFluxAna < 0` (plain `<`), so ΔE = 0 passes and
  `FluxAnalyzor.f90 setupFluxAnalyzor` divides by it (a garbage energy count, no flux results).
  C++: error if DeltaEFluxAna ≤ 0 (tolerant). Details: [fortran_input_bugs.md](fortran_input_bugs.md#flux-q2-a-zero-deltaefluxana-passes-the-check).
- Not ported (yet): the flux analysis (`FluxAnalyzor.f90`). `doFluxAnalysis = true` stops a
  propagating run (new or restart) with exit code 3 once the input checks have passed;
  analysis-only (`onlyAnalysis`) and potential-only (`Potential`) runs accept the keys and ignore them.
  A port was validated against the Fortran and reverted (2026-11-14): with ProjOp 2 the Fortran's
  flux leaves out the initial channel. Findings, tests and what was tried:
  [flux_analysis_findings.md](flux_analysis_findings.md).
- `InitGenComplexTransform` silently returns on invalid nφ.
- `tolsp` (`TolSp`) has no default in `globaldata.f90` (uninitialised if absent). C++: required key.
- `nZSpecular < nZ` is not checked by Fortran (out-of-bounds zgrid access). C++: hard error.
- Kept on purpose (Fortran semantics reproduced): `MJMax <= 0` ⇒ the basis/grid uses MJMax = JMax
  (`setupdvr.f90`/`transform.f90` "MjMax < 1"), while ValidateInput checks |MJIn| ≤ |MJMax| with the
  raw value. `VAna < 1` and absent `JAna` are replaced by the program-determined open-channel limits
  (`GasPhase.f90`).
- Not ported: `.pes` file caching (`ReadPesFromFile`/`WritePesToFile`); the potential is always
  computed. `UseRCutOff` is silently forced off by Fortran `ProcessInput`; C++ rejects
  `UseRCutOff: true` / `RCutOff`. Also rejected (NotSupportedError): `SystemMix` (+`SystemID2`,
  `SystemMixCoeff`), `write2DPES` (+ coordinates), `PESMixingAlfa`.
- Surface symmetry, deviations: the Fortran switches `UseSymmetry` off with a warning for C2v and
  C4v (complete but never validated code) and has no code for C3v/C3; C++ supports all five
  (gates for those: C++ with vs without symmetry). nX or nY = 1 with UseSymmetry and a trivial
  usable group are errors (Fortran: warn and continue). Representatives differ from the Fortran
  regions (triangle A / 1, rectangle A): equal values up to PES rounding. The PES self-check (θ, Z and r at grid points, x, y, φ anywhere; until 2026-11-09 θ was sampled anywhere too, which could reach geometries the grid never uses)
  (`SymmetryTolerance`, default 1e-5 eV) is not in the Fortran. Native checkpoint format 2 adds
  the number of φ points to the restart parameters.
- Diffraction order: derived from the reciprocal lattice (nearest-neighbour step distance,
  `SurfaceLattice::diffraction_order`) instead of the Fortran per-point-group table
  (`getDiffractionOrder`, which stops for other groups); identical integers for all Fortran
  groups (exhaustive test).
- Persson PES quirk, kept: coordinates are converted with 0.529167 Å/bohr, the lattice constant
  with cboan 0.529177, so the PES has period cc/0.529167 ≠ alat (1.9e-5 relative). Symmetric about
  the origin, but not periodic on the grid: with UseSymmetry the self-check fails (hundreds of eV
  deep in the wall, ~4 meV with VMax 10 eV); the Fortran would use the map silently.
- Dai-Light PES property, kept: V jumps to the gas-phase value at Z ≥ zasymp = 15 (the surface
  terms are still ~1e-4..1e-3 eV at 14.99).
- Giving both `ZOptEnergy` and `ZOptFactor` (same for R/ZSp) is an error (Fortran silently ignores
  the factor).
- absorber D1: an absorber energy of 0 (`ZOptEnergy`, `ROptEnergy`, `ZSpOptEnergy`) means no
  absorption (strength 0). The Fortran computes the strength only for an energy > 0 and otherwise
  keeps the factor's default (`globaldata.f90 SetupGlobalDataModule`, e.g. 4.38e-3 hartree for
  ZOpt), so it absorbs without saying so. The low-Z absorbers' default strength is 0, so there both
  codes agree. Details: [fortran_input_bugs.md](fortran_input_bugs.md#absorber-d1-an-absorber-energy-of-0-keeps-the-default-strength).
- The low-Z absorbers of upstream 5.24 (keys `ZOptLow…`, `ZSpOptLow…`; `InputNew.f90 ReadInput`,
  `Absorption.f90 SetupAbsorpPot`). The Fortran's silent cases are errors here; details and
  demonstrations: [fortran_input_bugs.md](fortran_input_bugs.md#the-low-z-absorbers-of-upstream-524).
  - low-Z D1: Start > End is an error (the Fortran leaves the absorber off without a message).
  - low-Z D2: Start = End is an error (the Fortran switches on a zero-length absorber: a grid point
    equal to Start gets 0·(0/0) = NaN, and the wave function becomes NaN).
  - low-Z D3: an incomplete set of keys (only one of Start and End, or the energy or factor
    without them) is an error (the Fortran leaves the absorber off without a message).
  - low-Z D4: the set-up output shows the absorbers' own values. The Fortran's printout
    (`InputNew.f90 PrintInput`) chooses between energy and coefficient by the high-end
    `ZOptEnergy`, and its specular block prints the scattering grid's values. The index ranges
    printed here are the points inside (the Fortran's `Setup.f90` prints the last point below
    Start, and its index loop never tests the first grid point).
  - Kept (not deviations): Start ≥ ZStart and End ≤ Zana are compared exactly, with the Fortran's
    plain `<` and `>`; the energy check is tolerant (`.RLT.`), as in the Fortran. Where the
    specular low-end and high-end absorbers overlap, the low-end value replaces the high-end one
    (the Fortran's assignment order). The specular Z factors are applied at every Z point, as 5.24
    does (`Absorption.f90 PerformEXPOpticalOpSpec`); outside the absorbers the factor is exactly 1.
  - C++ only: the set-up warns when a grid point lies on different sides of an absorber's
    strongest end in the two codes (low end: Start; high end: End). The grid points differ by an
    ulp (the quad-precision entry below), so there one code multiplies by exp(−1.5·A·dt) and the
    other by 1. The set-up recomputes the Fortran's point in quad precision (`fortran_grid_point`,
    `src/grids/spatial_grids.cpp`, with `__float128`).
- The check `ZSpOptStart < ZAsymptotic` is unreachable (implied by the Z0 checks); kept for parity.
- `-fdefault-real-8` also promotes `DOUBLE PRECISION`, `DFLOAT`, `DBLE` and `…d0` literals to
  REAL(16) (quad precision). Expressions such as `rmin + dfloat(i-1)*delr` (Setup.f90 grids) are
  therefore evaluated in quad and rounded once to double, i.e. correctly rounded (verified: the
  Fortran r grid differs from plain-double `rmin + i*dr` by 1 ulp at 13 of 56 points). Accepted by
  the project as within numerical limits: C++ keeps plain double grids. Effect on gate 1: vasy table
  agrees to ≤1e-14, the full 6D potential to ~5e-14 relative (gate tolerance 1e-13).
  The same applies to every momentum quantity in Setup.f90 (lines ~129–131, 454–525): delkz,
  delksp, delkr, the kinetic energies ekinz/ekinr/ensp used in exp(−iT·dt), and zk/rkg/zksp. These
  are ~1 ulp off in C++ too (≈1e-12 phase after 7000 steps): first suspect if a propagation gate
  drifts. When porting any routine, look for `DFLOAT`, `DBLE`, `…d0` literals and
  `DOUBLE PRECISION`; all are quad precision in this build.
- gfortran `-O3` vectorises some loops with glibc libmvec (`_ZGVbN2v_exp/sin/cos/acos`, up to a
  few ulp from the scalar functions): `angulargrid.f90` acos (θ grid), `legendre.f90`
  AssociatedLegendre exp, and the cos/sin of exp(−iV·dt) in `Potential.f90`
  PerformEXPInteractionOp, and the cos/exp of `PESForCu111Persson.f90` (Dai-Light and the hard
  wall are scalar). C++ uses the scalar functions, so these agree only to a few ulp.
  Real`**`integer is libgcc `__powidf2` (repeated squaring): use `fortran_integer_power`, not
  `std::pow`.
- More quad-precision sites (same policy, plain double in C++): `erot = dfloat(j(j+1))/(...)`
  (GridEnergies.f90; in practice bitwise equal), the CAP profiles `vimag/vimar/vimas` and their
  low-end counterparts (complex × `dfloat(nabpow+1)` → complex(16), Absorption.f90), `scalFacX/Y = sqrt(delx/dfloat(nx))`
  (KineticOp.f90), `kinx(nx/2+1)` (GridEnergies.f90) and `conx/conxy` in TwoGrids.f90
  (MomentumShifting.f90 computes the same phase with integer `(ix-1)`, i.e. in double).
- Analysis quad-precision sites (plain double in C++): `the_time = tottym + dfloat(propStep -
  nstep)*dtlan` (Main.f90) and `ensctr = enrgi + dfloat(ien-1)*denrgy` (Analyzor.f90 SetupAnalysis).
- Kept on purpose (`Analyzor.f90 BalintKurtiCoeffs`): with a fixed specular analysis surface the BK
  coefficient of the initial channel is *replaced* by the specular-grid projection (the
  scattering-grid part of that channel is discarded). The phase of S is set to 0 when Re S = 0
  (`smatrix_combine`), even if Im S ≠ 0.
- Fortran restart bug, fixed: `anlyTol` and `oldValZanly` are neither saved in fort.9 nor
  recomputed on a restart (uninitialised while the dynamic specular surface is not yet fixed).
  C++ native checkpoints store the full surface state; reading Fortran files recomputes anlyTol
  from the input's initial wave packet and sets oldValZanly = izanlySpecChan (what a continuous
  run has). A Fortran restart reconstructs the BK times as istep·deltat.
- Fine analysis, fixed (`docs/fortran_fine_analysis_bugs.md`): `smatrix_calc` (nread 2, non-NAG)
  uses ATAN instead of ATAN2 and takes the signs of Re/Im S separately from the linear
  interpolation. Only S and phases at the analysis energies are affected (probabilities: P
  interpolated linearly, unaffected). C++: P linear, phase = arg of the linearly interpolated S.
  5.19 also interpolated Im S with `s2_imag - s1_real` (should be `s1_imag`); the reference has
  that typo fixed (svn r1713), so gate 5 emulates only the ATAN and sign handling. Fortran
  nen95 = 1 reads out of bounds (C++: no interpolation for nen95 < 2); ndatol < 2 gives an FFT
  length 1 (C++: error).
- Fine analysis, kept: truncation to ndused = largest FFT length <= records (trapezoid, end
  weight 1/2, so not equal to the on-the-fly rectangle sum); open channel if kz^2 >= 0 (on the
  fly: > 0); interpolation index by plain `>=` while the range check is tolerant; a channel
  closed at E_z but open at a neighbouring grid point keeps its small interpolated P (in the
  sums); nen95 = 0 gives a warning and no results. On a grid point the Fortran ensctr (quad
  precision) may choose the neighbouring interval. Deviations: Fortran ignores the fort.22
  header (C++ checks it, and that records lie every dtAna from t = 0); `alsoAnalysis` in
  Fortran replaces the on-the-fly results, C++ writes both (`fine_analysis/`). The transform is
  one zero-padded FFT with trapezoid weights instead of hlftrf's two (anti-)Hermitian FFTs
  (equal to rounding, tested).
- `NrOfEs = 1` divides by zero in `denrgy = (EMaxWP-EMinWP)/(nenrgy-1)` (Fortran and C++ alike).
- ProjOp = 1 (Fortran `CheckOneGridContinuation`/`ProceedWithOneGrid`) is not ported (legacy
  option); the input reader rejects it.
- Kept on purpose (`Analyzor.f90 FindNewAnalysisPoint`): the clamp of `izanlySpecChan` to `nschmi`
  is overwritten on the next line, and reaching `izanly` fixes the surface without leaving the
  loop. The rest norm of the averages depends on this dynamic analysis surface.
- Reference-case property (not a bug, same in Fortran): the FGH grid is periodic, so the decay of
  ψ_v below r_min wraps around to r_max (ψ_00(r_max) ≈ 7e-4, a ~1e-7 plateau in between). The r
  absorber (ROptStart = 4.55) removes this tail, ≈ 5e-8 of norm per 10 a.u. of propagation.

## CRP library (`src/crplib/`, SystemID `CRPlibrary`)

The port of the Fortran CRP library (`CRPlibrary.dir/`, docs/crp_library.md). V and the analytic
gradients are bitwise equal to the Fortran library on the synthetic test data sets
(`tests/crplib/`), except where listed.

- Kept on purpose (Fortran semantics reproduced):
  - The C6v atom-surface Fourier table (`CRP_C6v_AtomSurface.f90 FOUR3D`) has the truncated DATA
    literals `0.666666666666667` and `1.333333333333`: it is the inverse of its site basis only to
    8.4e-14, so the atom PES misses its own data at the sites by ~1e-13 relative.
  - `ORIGEN`/`ORIGEN2` fold with `INT(ori)` and `INT(ori - 1)` for negative values: a point at an
    exact negative lattice coordinate lands on the far edge of the cell (u = 1, not 0). Equal
    values by periodicity, different rounding.
  - The library's own unit factors (`UnitConversion.f90`: 27.2113961 eV/Hartree, 0.529177249 Å/bohr)
    and their composition `x * ((1/c_in) * c_out)`; SPO-DVR converts the PES (eV) back with
    chaev = 27.2114, so data given in Hartree are scaled by 27.2113961/27.2114.
  - C4v `MURO`: a negative first interlayer spacing is a relative value (z = −IS(1)·a/√2), but it is
    unit-converted as a length like every spacing.
  - The cos/sin(π/3) and cos(π/2) that gfortran folded at compile time are hexadecimal literals
    (`crplib/constants.hpp`); they equal the run-time values on this machine (tested).
- Errors (DataError) where the Fortran misbehaves silently or crashes: a missing file of a used
  molecular cut (the Fortran skips it and later evaluates an unallocated spline), wrong numbers of cuts
  (atom: exactly 7/10 for C6v/C3v, else uninitialised or out-of-bounds sums), spline tables with
  fewer than 4 points (`DSPLIN` writes `cm(0)`) or more than 999 (`SPLIN2D` stops), grids that are
  not strictly increasing, an unknown symmetry string (the atom PES continues with garbage),
  ZsoMin > ZsoMax, cut file names longer than 50 characters (truncated by `CHARACTER(50)`), a
  `<data set>.dat` file (multi-PES branch, empty in the Fortran), C3v with fewer than 4 layers or
  Disp1 ≠ 0 and C4v with fewer than 2 layers. A folded coordinate beyond ±2·10⁹ lattice constants
  (`INT` overflow) is an EvaluationError. Equal ZsoMin and ZsoMax are accepted (as in the Fortran).
- Kept on purpose: a listed extra cut (beyond the ones the interpolation uses) may lack its file.
  The Fortran skips it (`INQUIRE EXIST`): no data, no limit on Rmin, Rmax, Zmin; so does the C++,
  with a warning. The real SRP48 data set (`potentials/H2Cu111SRP_C3v.dir`) lists two such cuts.
- Few-ulp difference (decision D4): with CutsInterpolation = F the Fortran subtracts the atom
  potentials from the cuts using cos θ, cos φ computed by gfortran's vectorised libmvec
  (`_ZGVbN2v_cos`); the C++ uses `std::cos`. That is a few ulp of the data values, carried into
  every result. On the synthetic reference sets: |ΔV| ≤ 4.4e-16 Hartree, |ΔF| ≤ 2.7e-15
  Hartree/bohr at under 12 % of the points; on the real SRP48 data set (data up to 72 Hartree)
  up to 1.9e-14 at a few points. The tests compare to 4e-15·max(unit·S, |value|), S the largest
  |data| in Hartree (at least 1). A scratch build calling `_ZGVbN2v_cos` was bitwise equal
  everywhere, also on the real data set (2026-11-11). CutsInterpolation = T data sets are bitwise
  equal.
- In the program (SystemID `CRPlibrary`, `src/potential/crp_potential.*`): `CRPlibDataSet` is
  required and is an error with another SystemID (the Fortran reads and ignores it). Data-set
  paths longer than 50 characters work: the Fortran opens `CRPData` with the full path but every
  other file through the `CHARACTER(50)` `PES%DataSet`, so it fails there. The library's STOPs
  (exit status 0 in the Fortran) and data-set errors are input errors naming `CRPlibDataSet`
  (exit status 2); an evaluation error names the grid point. The setup output reports the clamped
  r and Z grid points and the (Z, r) points where the library silently returns V = 0 (shifted to
  −VH2(r_last) by the wrapper; reproduced, decision D2). With `UseSymmetry` the potential of C3v and C4v data sets
  is computed on the symmetry-reduced grid (the Fortran switches the symmetry off for them); it
  differs from the full potential by rounding (docs/crp_library.md "Speed").
- Kept on purpose: for a vertical molecule (θ = 0 exactly) RRINTII/RRINT set dθ/dX = dθ/dY = 0,
  but the sin 2θ terms of the hollow and t2 site functions give V a lateral slope there; the x and
  y gradient components miss it (z is right). Irrelevant for SPO-DVR (V only; no θ grid point is 0).
- Molecular cut files whose rows repeat r or Z inconsistently: the Fortran keeps the last row's
  values for the grid but subtracts the atom potentials (CutsInterpolation = F) at each row's own
  (r, Z); the C++ uses the final grid values for both. Identical for consistent files.
- Not ported: the C4v 1D splines of the atom cuts (computed but never evaluated by the Fortran).

## SCM distorted-surface PES (`src/potential/scm/`, SystemID `CRPDistrtSurface`)

The port of the Fortran's distorted-surface correction (`DistortedSurfacePES.dir/`). Ported so far:
the readers of its two data files, the parameter file and the slab position file
(`scm_files.*`, tests `tests/test_scm_io.cpp`); the potential itself cannot be selected yet
([upstream_sync.md](upstream_sync.md)). Values and atoms read from the example files and from
synthetic multi-configuration files are bitwise equal to what the reference's `SCM_IO.f90` reads
(checked 2026-11-14 with a one-off gfortran driver against the r1713-dce0ecc objects).
Open physics questions: [scm_open_questions.md](scm_open_questions.md).

- SCM D5: position-file versions. The Fortran reads versions 1, 2 and 3 (`SCM_IO.f90 ReadSCMSlab`);
  version 1 generates the ideal lattice from the order of the atoms, with bugs (an INTEGER bulk
  layer spacing, `uCelly` used for `uCellx`, two layers hard-coded). The C++ reads version 3 (the
  example file) and version 2; version 1 is an error naming `PosFilepath` that advises version 3.
  Test: `ScmPositionFileErrors.ScmD5OnlyVersions3And2`.
- SCM D6: where the Fortran continues silently, the C++ stops with an error or does not cut:
  - a configuration number (`PosFileSlabNum`) ≤ 0: the Fortran reads configuration 1 (its skip
    loop in `ReadSCMSlab` runs zero times). Test: `ScmPositionFileErrors.ScmD6ConfigurationAtLeast1`;
  - `MaxNumSurfAtoms` ≤ 0 in a version 2 file: the Fortran uses no atom, ΔV = 0 everywhere. Test:
    `ScmPositionFileErrors.ScmD6MaxNumSurfAtomsAtLeast1`;
  - a version 2 configuration without an atom of layers 1–2: the Fortran's `ReadSCMSlab` (CASE(2))
    has no zero-atom check, so ΔV = 0 (version 3 stops in both, "No atoms in SCM range"). Test:
    `ScmPositionFileErrors.ScmD6Version2ConfigurationWithoutLayer1Or2Atom`;
  - a layer index ≤ 0: the Fortran keeps the atom and indexes `Pars(0)`, out of bounds. The C++
    checks every atom of the configuration read, also one outside the cut-off. Test:
    `ScmPositionFileErrors.ScmD6LayerAtLeast1`;
  - P17 ≥ P18 in the parameter file: `SCM_Types.f90 Setup_3Body_Rydberg_Par` clamps r with
    IF/ELSE IF and no check, so every r ≥ P17 becomes P18. Test:
    `ScmParameterFileErrors.ScmD6ClampRangeMustBeIncreasing`;
  - the file paths: the Fortran cuts them at 100 characters (the `CHARACTER(100)` arguments of
    `InitialiseSCMModule`); the C++ uses them in full. Test: `ScmFiles.ScmD6LongPathsAreNotTruncated`.
- SCM D8: the Fortran's checks of the parameter file's first line, kept. It must start with
  `VERSION` in columns 1–7 (case-sensitive); the version is read list-directed from column 8.
  Without the tag the Fortran assumes version 1 and rejects it ("not (yet) supported"), as every
  version other than 2. The C++ rejects the same files, naming `SCMParFilepath`. Tests:
  `ScmParameterFileErrors.ScmD8VersionTag`, `.ScmD8OnlyVersion2`.
- SCM D14: slab size. The Fortran requires nz ≥ 2 (`UsedLayers .GT. SlabAtomSize(3)` in
  `ReadSCMSlab`); so does the C++, which also requires nx, ny ≥ 1 in version 3 (the Fortran
  multiplies the three sizes: two negative ones give a positive atom count). Not checked by the
  Fortran either, and needing the CRP data set: an FCC111 slab with a CRP set that is not C3v or
  C6v; that check belongs to the set-up of the potential, not ported yet. Test:
  `ScmPositionFileErrors.ScmD14NzAtLeast2NxNyAtLeast1`.
- Reader differences without an ID (tests in `ScmPositionFile.*` and `ScmParameterFileErrors.*`):
  - Not reproduced: version 3's temporary arrays in the Fortran hold nx·ny·2 atoms, so a file with
    more atoms of layers 1–2 inside the cut-off overflows them; the C++ has no such limit.
  - Version 2: the Fortran keeps atoms of layers below the second and skips them in the potential
    loop (`SCM_Potential.f90 SCM_SurfInteract_3Body`, before any arithmetic); the C++ drops them
    while reading. The same ΔV; only the atom count differs.
  - Kept: the atoms start on line 4 even when a header statement spans lines (the Fortran reads
    the header, re-opens the file and skips three records); the cut-off test is
    ((x² + y²) + z²) < c² on the distorted position, as in the Fortran (other forms, among them
    the FMA contractions of a fast-math build, decide boundary atoms differently).
  - The shared list-directed reader (`src/fortran_text/`) rejects input that gfortran accepts:
    repeat counts, null values, `/`, NaN and Infinity, and an exponent without its letter
    (`1.0+3`, which gfortran reads as 1000). It cuts neither line 1 of the parameter file at 100
    characters nor `SlabFileTag` and `SlabSurfStruct` at 80 (`CHARACTER(80)`), and compares the
    structure without blank padding: a quoted `'FCC111 '` passes the Fortran's check
    (`GetSurfaceParameters`) but not the C++ one.

## Upstream changes after 5.19

The port started from 5.19; the reference is now 5.24 at svn r1713. The port status of every
upstream change since 5.19 (ported, a deviation listed here, or not ported, and what such an input
does in the C++ today) is in [upstream_sync.md](upstream_sync.md).
