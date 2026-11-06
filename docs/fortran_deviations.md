# Deviations from the Fortran code

Every known difference between the C++ program and the Fortran reference (`~/work/SPO-DVR`,
version 5.19): Fortran bugs that the C++ fixes, Fortran quirks that it reproduces on purpose,
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
  `SystemMixCoeff`), `write2DPES` (+ coordinates), `PESMixingAlfa`, `CRPlibDataSet`.
- Surface symmetry, deviations: the Fortran switches `UseSymmetry` off with a warning for C2v and
  C4v (complete but never validated code) and has no code for C3v/C3; C++ supports all five
  (gates for those: C++ with vs without symmetry). nX or nY = 1 with UseSymmetry and a trivial
  usable group are errors (Fortran: warn and continue). Representatives differ from the Fortran
  regions (triangle A / 1, rectangle A): equal values up to PES rounding. The PES self-check
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
  P linear, phase = arg of the linearly interpolated S. Fortran nen95 = 1 reads out of bounds
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
