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
- flux Q3, kept: the number of flux energies is `INT((EMaxFluxAna − EMinFluxAna)/DeltaEFluxAna) + 1`
  with the energies in hartree (`FluxAnalyzor.f90 setupFluxAnalyzor`), so EMaxFluxAna is dropped for
  most inputs, by rounding even when DeltaEFluxAna divides the range: 0.2–0.6 eV in steps of 0.1 gives
  4 energies (0.2–0.5 eV), 0.40–0.95 eV in steps of 0.025 gives 22 (the last 0.925 eV). C++: the same
  count (`flux_energy_count`); `make_flux_energies` adds a warning when the last energy lies below
  EMaxFluxAna (`src/analysis/flux_analysis.cpp`; test `FluxEnergies.CountDropsTheEndPointAsTheFortran`).
- flux Q4, fixed: the Fortran weighs every analysis step of the flux sums with the input dtAna
  (`FluxAnalyzor.f90 storeTimeIndependentWF`, coeff = … · dtAnalysis), but a step lasts nanspo·dtStep
  with nanspo = INT(dtAna/dtStep) (`InputNew.f90 ProcessInput`), and nothing checks that dtAna is a
  multiple of dtStep. C++: the step's own length nanspo·dtStep (`derived.analysis_time_step`), as the
  S-matrix sums use (`Analyzor.f90 BalintKurtiCoeffs`, dtmul). The two are equal to the bit when
  nanspo·dtStep rounds to dtAna, as in every gate case (dtStep 2, dtAna 10). They differ otherwise
  (test `FluxPropagation.FluxQ4WeightIsTheEffectiveAnalysisStep`, the first case):
  - by a fraction of dtStep when dtAna is not a multiple of it: dtAna = 10.5 at dtStep 2 gives steps
    of 10 a.u., whose Fortran weights are 1.05 times too large;
  - by a whole dtStep when the quotient truncates: 0.3/0.1 = 2.9999999999999996, so dtAna 0.3 at
    dtStep 0.1 gives nanspo = 2, steps of 0.2;
  - by rounding only: dtAna 0.9 at dtStep 0.3 gives nanspo = 3 and 3 · 0.3 = 0.8999999999999999.
- flux Q5, kept: `FluxAnalyzor.f90` converts its energies, and the wave packet's EMinWP and EMaxWP
  for b(k), with `Constants_eV2Ha` = 0.036749308286 (`units::precise`); the initial wave packet and
  every other energy use chaev = 27.2114 (`units::legacy`). The two differ by 1.27e-7 relative, and
  |b(k)|² differs from the momentum distribution of the C++ initial packet by 2.9e-7 of its peak
  (hardwall_small). C++: the same factors (tests `FluxEnergies.FollowTheFortranExpressions` and
  `FluxMomentumAmplitude.IsTheMomentumDistributionOfTheInitialWavePacket`).
- flux Q6: with `onlyAnalysis` (nread = 2) the Fortran sets the flux analysis up
  (`FluxAnalyzor.f90 setupFluxAnalyzor` allocates φ and φ′ and, when NewCalc is false, reads the
  flux restart files) but computes and writes nothing (`Main.f90` calls `computeFluxProbability` only
  without onlyAnalysis). C++: an analysis-only run with `doFluxAnalysis` prints one note and has no
  flux arrays and no flux files (`src/app/run.cpp`; test
  `FineAnalysisRuns.FluxQ6AnalysisOnlyRunPrintsANoteAndNoFlux`).
- flux Q7, kept: the flux analysis leaves out the initial channel at r_c. At every point Z < zasymp,
  every r included, each potential step moves the wave packet's projection on the initial channel
  (rotational state J_in, mJ_in times the initial (x, y) plane wave, the specular diffraction order)
  to the specular grid (`Potential.f90 PerformEXPPotentialOp` calls `TwoGrids.f90 SubSpecFromScat`
  for every izr ≤ nzpot·nr), so the scattering wave holds none of it there; and
  `FluxAnalyzor.f90 storeTimeIndependentWF` takes ψ and ∂ψ/∂r at r_c from the scattering wave only
  (`getWaveFunctionAtZr`, `getWaveFuncDerivAtZr`). After the collision that channel also carries
  scattered and reacting parts. P(E) misses that channel's own flux (its cross terms with the rest
  vanish in the (x, y) sum, the plane waves being orthogonal); P(x, y; E) misses its own flux and the
  cross terms. The incoming packet χ_vin(r)·G(Z) carries no flux in r of its own (χ_vin is real), so
  of it only cross terms are missed, in P(x, y; E). C++: the same (`PotentialOperator`, its
  SubSpecFromScat block; gate 6 compares with the Fortran). The FLUX ANALYSIS set-up block prints
  |χ_vin(r_c)| / max |χ_vin| and warns above 1e-6; this ratio bounds only the incoming packet's part
  of what is left out. The size of the rest is not measured. Tests
  `FluxInitialStateRatio.IsTheSizeAtTheFluxSurfaceOverTheLargest`,
  `FluxReport.SetupBlockShowsTheMemoryRuleAndWarnsAboveTheFluxQ7Limit`.
- flux Q8, checked: the Fortran's flux restart files `flux_coeffs.dat` (φ) and `flux_coeffs_deriv.dat`
  (φ′) have no header: one gfortran record φ(:, iz, n) per plane Z < zasymp and flux energy, energy
  outer (`FluxAnalyzor.f90 writeFluxCoeffs`). `readFluxCoeffs` reads them unchecked: longer records
  and extra records are read silently, so a file written for another energy grid or rotational
  basis is read shifted; only too few or too short records stop it (a gfortran run-time error).
  C++: `read_fortran_flux_coefficients` (`src/io/fortran_files.cpp`) requires exactly
  nE × nz_potential records of (rotational functions × nX × nY) complex values each, and otherwise
  stops with an error naming the expected sizes (test
  `FortranFluxFiles.RecordCountAndLengthMustFitTheInput`).
- flux Q9, kept: b(k) (`MyFunctions.f90 getBkCoeffWavePacket`) has the phase e^(+i·δk·z0), δk = k0 − k,
  the initial wave packet's momentum amplitude e^(−i·δk·z0): b is its complex conjugate. P(E) does not
  depend on it (a phase per energy cancels in conj(φ)·φ′), but the energy-resolved φ and the Fortran's
  flux restart files do. C++: the same phase (`wave_packet_momentum_amplitude`; test
  `FluxMomentumAmplitude.IsTheMomentumDistributionOfTheInitialWavePacket` checks b = conj of the
  packet's Fourier transform).
- Not ported (yet): the restart of the flux analysis (`FluxAnalyzor.f90 readWriteFluxCoeffs`). A
  restart with `doFluxAnalysis = true` stops with exit code 3 once the input checks have passed; a
  new run computes the flux analysis (gate 6, [validation.md](validation.md));
  analysis-only (`onlyAnalysis`) runs accept the keys and print a note (flux Q6), potential-only
  (`Potential`) runs accept them and ignore them.
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
  wall are scalar; so are b(k) and the r-derivative table of the flux analysis: `pow`, `cexp` and
  `tan` in the 5.24 objects, where `(-1.0)**distance` is a parity select, ±1 as with `__powidf2`). C++ uses the scalar functions, so these agree only to a few ulp.
  Real`**`integer is libgcc `__powidf2` (repeated squaring): use `fortran_integer_power`, not
  `std::pow`.
- Signed zeros, accepted (not reproduced; values equal): gfortran multiplies a complex number by a
  real x as by the complex number (x, 0); the C++ multiplies its two parts by x. The values are the
  same; only the sign of an exactly zero part can differ: (−0, −1)·1 has the real part +0 with
  gfortran 15.2 (bits 0000000000000000) and −0 in the C++ (8000000000000000), 2026-11-12. This holds
  for every real × complex product of the port, e.g. the FFT scaling (`v = scale_ * v` in
  `src/fft/fftw_fft.cpp`, `Y = Scale * Y` in `fftw_2_scs.f90`) and b(k), the r-derivative and c_n of
  the flux analysis. Compare such data by value (==), not by bytes.
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

## Upstream changes after 5.19

The port started from 5.19; the reference is now 5.24 at svn r1713. The port status of every
upstream change since 5.19 (ported, a deviation listed here, or not ported, and what such an input
does in the C++ today) is in [upstream_sync.md](upstream_sync.md).
