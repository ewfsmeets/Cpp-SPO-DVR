# Deviations from the Fortran code

Every known difference between the C++ program and the Fortran reference (version 5.19, upstream
svn r1572, in `~/work/SPO-DVR_v5.19`): Fortran bugs that the C++ fixes, Fortran quirks that it reproduces on purpose,
features that were not ported, and places where the Fortran computes in quad precision while the
C++ uses double. Add every new deviation here. The fine-analysis bugs are described in detail in
[fortran_fine_analysis_bugs.md](fortran_fine_analysis_bugs.md); the gates that bound the numerical
differences are in [validation.md](validation.md).

- `nZSpecular` read depends on `nproj` before `ProjOp` is read (always uses default 2).
- Error message `'EMaxFluxAna < EMinWP'` for the check `EMaxFluxAna < EMaxWP`.
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
  (GridEnergies.f90; in practice bitwise equal), the CAP profiles `vimag/vimar/vimas`
  (complex × `dfloat(nabpow+1)` → complex(16), Absorption.f90), `scalFacX/Y = sqrt(delx/dfloat(nx))`
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
  interpolates Im S with `s2_imag - s1_real` (should be `s1_imag`), uses ATAN instead of ATAN2,
  and takes the signs of Re/Im S separately from the linear interpolation. Only S and phases at
  the analysis energies are affected (probabilities: P interpolated linearly, unaffected). C++:
  P linear, phase = arg of the linearly interpolated S. Upstream fixed the `s1_real` typo on
  2026-11-09 (svn r1713, after 5.24); the ATAN and sign handling are unchanged there. Fortran nen95 = 1 reads out of bounds
  (C++: no interpolation for nen95 < 2); ndatol < 2 gives an FFT length 1 (C++: error).
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

## Upstream changes after 5.19 (not ported)

The Fortran code continued upstream, to version 5.24 (2025) and a fix on 2026-11-09 (svn r1713).
`git diff -w 601d833 master -- sources.dir/SOURCES.dir` in the upstream repository lists them
(601d833 = svn r1572 = 5.19, identical to the reference sources). The C++ follows 5.19 throughout.
The changes below come from reading that diff; the 5.24 code has not been built or run here.

**Changes for inputs that use only 5.19 features.** The operators, absorber formulas, FFT paths
of the build, input defaults and the formats of the averages, probability and BK files are
unchanged. What differs:
- Build: the GNU make file (`GNU.mk`, formerly `GFORTRAN.mk`) uses
  `O3FLAGS = -Ofast -mtune=native -march=native` instead of `-O3`, i.e. fast-math and FMA. Results
  of a default 5.24 build are therefore not comparable at gate level and depend on the CPU. Build
  with `O3FLAGS = -O3` to compare.
- `InputNew.f90`: `RStart <= 0` is rejected (5.19: `< 0`). The C++ keeps the 5.19 check.
- `Averages.f90`: the norm, ⟨Z⟩ and ⟨r⟩ reductions got `SCHEDULE(STATIC,1)`, so their partial
  sums are split differently (last-digit changes in `averages.out`). The C++ uses a fixed order.
- `Analyzor.f90`: the `s1_real` typo of the fine analysis is fixed (r1713). The C++ already fixed
  it, see above.
- OpenMP `DYNAMIC` scheduling replaced by `STATIC,1` everywhere: no effect on independent loops.

**New features.**
- Optical potentials at the low end of the Z grid, scattering and specular:
  - Keys: `ZOptLowStart`, `ZOptLowEnd`, `ZOptLowEnergy` or `ZOptLowFactor`, and `ZSpOptLow…` for
    the specular grid. They are off by default.
  - The value overwrites any high-end absorber at the same Z.
- PES morphing (`Potential.f90 MorphPotential`). It has no input key and runs when
  `pesmorphparameters.dat` is in the run directory. It scales or shifts FBR coefficients with
  Gaussians in Z and r and needs MjMax = JMax. From reading the code: it runs before the `.pes`
  cache is written, the cache name does not record it, and its loops assume the full (x, y) grid
  (`UseSymmetry`).
- Distorted-surface (SCM) CRP potential, SystemID `CRPDistrtSurface`. Keys:
  `DistortedPESMode`, `SCMParFilepath`, `PosFilepath`, `PosFileSlabNum`, `EnableSurfExpan`,
  `SCMCutoff`. In the CRP library this adds only an accessor for the interlayer spacings
  (`CRP_QD_Wrapper.f90`).
- PES `H2W110` and `HW110` (vdW CRP, W(110)). PES4Cu100 and every NAG dependency were removed.
- FFTs through oneMKL (Intel compiler only). In-place r2c/c2r and the real-to-real FFTs of
  `fftw_2_scs.f90` now stop with an error (unused in 5.19 runs).
