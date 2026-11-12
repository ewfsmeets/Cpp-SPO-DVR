# CRP library (`src/crplib/`)

A C++ port of the Fortran CRP library of M. Wijzenbroek (`CRPlibrary.dir/` of the Fortran code):
potential energy surfaces of a diatomic molecule at a surface from ab-initio data, interpolated
with the Corrugation Reducing Procedure (CRP). In SPO-DVR it is the PES of `SystemID` `CRPlibrary`;
the input key `CRPlibDataSet` names the data set ([input_schema.md](input_schema.md)).

- **Scope:** the atom-surface PES and the molecule-surface PES of C3v, C4v and C6v surfaces, with
  V and analytic gradients; the general entry points with unit conversion (`CRP_Wrapper.f90`); the
  quantum-dynamics wrapper (`CRP_QD_Wrapper.f90`, V only) that SPO-DVR uses.
- **Not ported** (unused by SPO-DVR): the 7D wrapper (surface coordinate Q), multi-PES data sets
  (a `<data set>.dat` file; that branch is empty in the Fortran), MPI synchronisation, PES mixing,
  and the routines no active code path calls (`CHECKR`, `DPL2D`, `DPLCOFT`, the C4v 1D splines).
- **Self-contained:** CMake target `spodvr_crplib`, namespace `crplib`, standard library only.
  It keeps the Fortran library's own constants and unit factors (below). A CTest
  (`CrpLib.IncludeBoundary`) checks that the library includes only its own headers and that only
  the adapter `src/potential/crp_potential.cpp` includes the library.
- **Exact:** the library is compiled with strict floating-point flags in every build (also the
  fast-math one). V and gradients are bitwise equal to the Fortran library on the synthetic test
  data sets; the one exception (a few ulp, data sets with CutsInterpolation = F) and the Fortran
  quirks that are kept are listed in [fortran_deviations.md](fortran_deviations.md).

## Public API

All evaluation is `const`, allocation-free and thread-safe; tables are built by the constructors.
Internal units are bohr and Hartree.

| Header | Content | Fortran |
|---|---|---|
| `crplib/dataset.hpp` | `DataSet` (an exact image of the files, in their own units), `read_data_set`, `write_data_set`, `validate_data_set`, `data_set_warnings` | the reading parts of `SetUp_CRP_*_PES`, `Read_*_Data` |
| `crplib/atom_surface.hpp` | `AtomSurface::evaluate(position)` → V and gradient of one atom; `make_atom_surface(data)` | `CRP_C{3,4,6}v_AtomSurface.f90` |
| `crplib/molecule_surface.hpp` | `MoleculeSurface(data)`, `evaluate(a, b)` → V and the gradients for both atoms, `nullopt` where the Fortran gives up silently; the data ranges, the gas-phase potential, the switching function | `CRP_C{3,4,6}v_MoleculeSurface.f90` (`POT6D`), `CRP_GasPhase.f90` |
| `crplib/wrapper.hpp` | `calculate_molecule_surface`, `calculate_atom_surface`: positions in, V and gradient out, in any unit set | `CRP_Wrapper.f90` `Calculate_CRP_*_PES` |
| `crplib/qd_wrapper.hpp` | `QdWrapper::potential_ev(u, v, Z, r, θ, φ)` (skewed lattice coordinates, bohr; eV), `parameters()`, `range_report(z, r)` | `CRP_QD_Wrapper.f90` |
| `crplib/units.hpp` | unit codes and conversions | `UnitConversion.f90` |
| `crplib/errors.hpp` | `DataError` (data set unusable; names the file and line), `EvaluationError` (a geometry where the Fortran STOPs) | the STOP sites |

The Fortran passes positions and gradients as arrays whose component orders differ:
R = (Za, Xa, Ya, Zb, Xb, Yb) but F = (dXa, dYa, dZa, dXb, dYb, dZb). The C++ uses named
`Vec3` fields instead. Cartesian x runs along the first lattice vector, z along the surface normal,
and the top layer is at z = 0.

## Data-set format

A data set is a directory. Every file is read with Fortran list-directed input (`READ(u, *)`): one
READ takes the next items from one or more lines and skips the rest of its last line. Separators
are blanks, commas or tabs; `D` exponents are accepted. Repeat counts (`3*1.0`), null values
(`,,`) and `/` are rejected with an error (the Fortran would accept them).

| File | Content (one READ per item) |
|---|---|
| `CRPData` | symmetry `C3v`, `C4v` or `C6v` (exact spelling); length and energy unit codes of all files; lattice constant (nearest-neighbour distance); number of layers N; N − 1 interlayer spacings (one READ, also when empty); optionally `Disp1 Disp2` (atom PES only, else 0) |
| `CRPData.info` (optional) | three lines: short description, long description, reference |
| `CRPData.3D` | number of atom-surface cuts; per cut `x y file` (lattice units, Cartesian) |
| `CRPData.6D` | number of molecule-surface cuts and CutsInterpolation (`T`/`F`); per cut `x y θ φ file` (lattice units: skewed (u, v) for C3v and C6v, Cartesian for C4v; degrees) |
| `Atom/intrep.dat` | dV/dρ at the first and at the last point; number of points; per point `ρ V` (the corrugation function of the atom-atom distance ρ) |
| `Atom/<cut>` | number of points; dV/dz at the first and at the last point; per point `z V` |
| `Molecule/extrapol.dat` | ZsoMin; ZsoMax; number of points; per point `Z V` (the Z dependence of the gas-phase potential; its last Z is Zasy) |
| `Molecule/VH2.dat` | number of points; per point `r V` (the free molecule) |
| `Molecule/<cut>` | number of r points; number of Z points (≤ 999); per point `Z r V`, r running fastest |

Unit codes: length 0 = bohr, 1 = Å; energy 2 = Hartree, 3 = eV, 4 = kcal/mol, 5 = kJ/mol,
6 = J. The library's factors are 0.529177249 Å/bohr and 27.2113961 eV/Hartree; kcal/mol 627.50956,
kJ/mol 627.50956·4.184 and J 4.35974434e-18 per Hartree. The program converts the PES (eV) to
Hartree with its own factor 27.2114, as the Fortran does.

The data must satisfy these conditions; `read_data_set` checks them and names the file:
- Atom-surface cuts: exactly 7 (C6v) or 10 (C3v), at least 6 (C4v; all with the same z grid).
- Molecule-surface cuts: at least 18 (C6v), 29 (C3v) or 16 (C4v). Further cuts are not used, but
  their grids still limit the covered range (a warning). A further cut may lack its file: it is
  then skipped and does not limit the range, as in the Fortran (a warning).
- Every table has at least 4 strictly increasing points; cut file names have at most 50
  characters; ZsoMin ≤ ZsoMax; C3v needs N ≥ 4 layers and Disp1 = 0, C4v N ≥ 2.
- A molecule-surface cut file repeats r and Z in every line; the last line with a given index wins.

The cut order is fixed: the interpolation expects cut k at the site and orientation of the
appendix. With CutsInterpolation = F the cuts hold V and the library subtracts the atom-surface
PES of both atoms at every point; with T they hold the interpolation function I = V − Va − Vb.

## Evaluation

V = Va + Vb + I below ZsoMin. Between ZsoMin and ZsoMax it is blended with the gas-phase potential
VH2(r) + Zdep(Z) by F = ½ + ½ cos(π (Z − ZsoMin)/(ZsoMax − ZsoMin)); above ZsoMax it is the gas-phase
potential, and Zdep is constant beyond Zasy. The atom-surface PES is the sum of the corrugation
function over a fixed cluster of surface atoms plus a smooth interpolated remainder:

| Symmetry | Cluster | Coverage radius | Remainder |
|---|---|---|---|
| C6v | 30 atoms of the top layer | 2.60 a | Fourier, 7 sites |
| C3v | 8 × 8 atoms in 4 ABC-stacked layers | 2.89 a (in-plane) | Fourier, 10 sites |
| C4v | 32 atoms in 2 layers | 2.06 a | periodic 3D spline, 6 sites |

The coverage radius is the largest atom-atom distance for which every surface atom near the cell
is in the cluster. If the corrugation function reaches further (its last ρ), V is not periodic
and jumps where the point is folded into the cell: the library warns.

### Outside the data

| Where | Fortran | C++ |
|---|---|---|
| r outside [Rmin, Rmax], Z < Zmin (QD wrapper; the ranges covered by all molecule-surface cuts) | evaluated at Rmin + 10⁻¹², Rmax − 10⁻¹², Zmin + 10⁻¹² bohr | the same; the setup output counts the clamped grid points |
| (r, Z) of the molecule outside a cut grid, Z < ZsoMax; r outside VH2.dat above it | V = 0 and F = 0, silently (the QD wrapper then returns −VH2(r_last)) | the same (`evaluate` returns `nullopt`); the setup output counts these (Z, r) points |
| an atom outside the z range of the atom-surface cuts (C3v, C6v; below it for C4v); an atom closer to a surface atom than the first ρ of the corrugation function (C3v, C6v; C4v extrapolates) | STOP | `EvaluationError`; in SPO-DVR an input error naming `CRPlibDataSet` and the grid point |
| C4v atom above the last z of its cuts | V = 0 | the same |

So the atom-surface cuts must cover z from Zmin − Rmax/2 to ZsoMax + Rmax/2 (for a molecule at
any orientation), and the molecular cuts should cover the r and Z grid of the dynamics.

## Speed

On the dynamics grid the PES is evaluated tens of millions of times (microseconds each; timings in
[performance.md](performance.md) "CRP library").
- Use `UseSymmetry` (with `nPhisOptimal`): the PES is then evaluated only at the symmetry-unique
  (x, y) points, 4 to 8 times fewer on a 12 × 12 grid. The Fortran switches the symmetry off for
  C3v and C4v. The interpolation is symmetric to rounding (the symmetric images run through
  different fold and Fourier arithmetic): on the synthetic data sets the gathered potential equals
  the full one bitwise at the representatives and to ≤ 1.6e-13 of max(|V|, 1 eV) elsewhere, and
  the self-check finds ≤ 6e-13 eV. Leave it off only for bitwise comparisons with the Fortran.
- Gate 1 also runs with `UseSymmetry` against the full Fortran potential
  (`CrpSymmetric/RegressionPes`, [validation.md](validation.md)).

Inside the library, the speed-ups leave every result bitwise unchanged (`crplib::EvaluationOptions`;
the tests compare with them switched off, `tools/benchmark/crp_timing.cpp` times both):
- **Cluster pruning** (S2). The repulsive sum (MURII, MURO) runs over a fixed cluster (C6v 30
  atoms, C3v 256, C4v 32), but only atoms within the last point of the corrugation function
  contribute. The library sums over the atoms that can be in reach of a point of the fold cell
  (chosen once) and skips, per call, those whose height differs by more than the reach. This is
  exact: an atom out of reach adds +0 to V and ±0 to the gradient sums, which start at +0 and never
  become −0 (in round-to-nearest x + y = −0 only if both are −0), so leaving it out changes no bit;
  the other atoms are summed in the cluster order. The computed distance is never below the
  computed height difference (monotone rounding), so the z test is exact; the in-plane selection
  keeps a relative margin of 10⁻¹⁰ and 2·10⁻⁸ lattice constants around the cell, and a point that
  the fold left outside the cell (rounding far from the origin) uses the full cluster. The reduced
  atom cuts are built with the full cluster. On the synthetic data sets (cut-off 6 bohr) C3v keeps
  54 of 256 atoms. Verified bitwise (`tests/crplib/test_cluster_pruning.cpp`): dense points in and
  on the fold cell for cut-offs from 0.35 a to 10 a (where nothing is pruned) and every cluster
  variant; points at the cut-off distance of an atom, to the last ulp, with a corrugation function
  that does not vanish there; the PES up to 10⁸ lattice constants from the origin; the QD wrapper
  on a dynamics grid; 4·10⁷ random points (`SPODVR_LONG_TESTS`); and the mutation check below.
- **Shared grid lookup** (S4). Every molecular cut is a bicubic spline on its own (r, Z) grid, and
  the Fortran locates (r, Z) in each of them (interval search and cell weights). Cuts with bitwise
  identical grids (usually all of them) share one lookup per evaluation; only the cell values
  differ. A point outside any grid still gives the silent V = 0, whichever cut it is. Grids are
  compared bit by bit (not with ==, which would merge −0 and +0 knots). Verified bitwise against
  the per-cut lookup on the synthetic data sets (one grid), with a cut with its own grid (and
  points outside only that one), and with grids that differ only in the sign of a zero knot.
- **V without the derivatives** (S3). The Fortran computes V and the gradient together and the QD
  wrapper discards the gradient. `MoleculeSurface::value` and `AtomSurface::value` run the same
  evaluation with `Want::value`: every statement of V is executed as written, the derivative
  statements are skipped (the d/dZ, d/dr and angular parts of the site functions and their sines,
  Vinter 2–7, the Fourier basis derivatives, DPLP, the gradient of the repulsive sum, the chain rule,
  TRANSF). Which sines a value needs was checked per routine (the hollow site's sin 3φ and the
  C3v sine stars of the atom PES enter V). Verified: V-only equals the full V bitwise at every
  Fortran reference point, per site function, and for every combination of S2, S3 and S4 in the QD
  wrapper on a grid that includes the silent branch.
- **θ trig once** (S5, built in). cos/sin 2θ and 4θ are computed once per evaluation and passed
  to the site functions (the Fortran recomputes them in each); the same expressions, so the same
  values (tested).

One thread, per QD-wrapper evaluation (`spodvr_crp_timing`, synthetic data sets; every V bitwise
equal):

| Data set | plain | S2 | S3 | S4 | all |
|---|---|---|---|---|---|
| C6v | 2.41 µs | 2.33 | 1.97 | 1.94 | 1.42 (1.7×) |
| C3v | 4.46 µs | 3.05 | 3.52 | 3.65 | 1.77 (2.5×) |
| C4v | 1.97 µs | 1.79 | 1.47 | 1.54 | 1.03 (1.9×) |

With 8 threads on a 12 × 12 grid the potential of the C3v set took 38.6 s in the Fortran, 11.6 s
here, 2.5 s with `UseSymmetry` ([performance.md](performance.md) "CRP library").

## Testing without a real data set

The library tests need no real data set (the one in the repository is checked separately, see
"The real data set" below). The tests (`tests/crplib/`, see its README) use:
- Properties that need no data: splines reproduce cubics, the Fourier tables invert their site
  bases, every site function returns its cut values at the cut orientations, homonuclear and
  surface symmetry, periodicity, gradients against finite differences.
- Synthetic data sets in the file format (`tests/crplib/crp_synthetic.hpp`): an exactly symmetric
  pairwise model for each symmetry, and "impulse" sets (one molecular cut equal to 1, everything
  else 0), whose response is known analytically.
- Bitwise references from the Fortran library: a driver
  (`tools/fortran_reference/crp_driver/`) linked against the objects of the Fortran reference
  install (`OBJS/`, [tools/fortran_reference/README.md](../tools/fortran_reference/README.md))
  evaluates point lists on the synthetic data sets
  (`tools/fortran_reference/run_crp_reference.sh`).
- Gates 1 and 2 of the program on the same data sets (`tests/reference/crp_*_synthetic_small`,
  [validation.md](validation.md)).

### Mutation check

`tools/crp_mutation_check/run.sh` applies deliberate bugs to a git worktree and runs the library
tests; each must make a test fail, except the removal of a redundant safety margin (marked in
`mutations.py`). Result (2026-11-11, 245 library tests; the numbers are the failing tests):

| Mutation | Caught by (number of failing tests) |
|---|---|
| bridge images swapped | impulse responses, Fortran comparisons (18) |
| 0.5 instead of the folded cos(π/3) in the QD wrapper's skew transform | QD wrapper on the impulse sets, bitwise (10) |
| 0.5 instead of the folded cos(π/3) in the C3v MURII | C3v atom-surface PES vs the Fortran (3) |
| SPLINT slope without its curvature term | spline tests, Fortran comparisons (47) |
| SPLINT instead of DPL for the atom-surface cuts | Fortran comparisons (6) |
| x / c instead of x · (1/c) (Å → bohr) | Å/eV cases vs the Fortran (10) |
| (x · (1/a)) · b instead of x · ((1/a) · b) (energy units) | unit tests (1) |
| `floor` instead of the ORIGEN truncation | fold tests, Fortran comparisons (9) |
| C3v t2h and t2f cuts swapped | impulse responses, real data set vs the Fortran (15) |
| sign flip in a chain-rule term | Fortran comparisons, finite differences (26) |
| TRANSF regions 4 and 5 swapped | C4v atom-surface PES vs the Fortran (4) |
| dFso term dropped from the gradient | Fortran comparisons, finite differences (13) |
| S2: the reachable atoms in reverse order | Fortran comparisons, pruning tests (17) |
| S2: in-plane selection with half the reach | pruning tests, representable model (40) |
| S2: per-call z test with half the reach | pruning tests, representable model (28) |
| S2: no fallback to the full cluster outside the fold cell | pruning test outside the cell (1) |
| S4: every cut evaluated in the cell of the first group | cut with a moved knot (1) |
| S4: grids compared with == (merges −0 and +0 knots) | grids differing in a zero knot (1) |
| S4: only the first group located | cuts with their own grid (2) |
| S3: sin 3φ of the hollow site skipped for V | Fortran comparisons (15) |
| S3: the C3v sine stars of the atom PES skipped for V | C3v atom-surface PES vs the Fortran (11) |
| S3: V-only blend without the gas-phase term | Fortran comparisons (16) |
| a missing file of an unused extra cut rejected | the real data set, the missing-cut tests (16) |
| validation accepts a used cut without a file | data-set test (1) |
| S2 safety margins removed: both; the relative one; the cell enlargement; `<` for `≤` | survive, as expected |

Three mutations first survived and led to new tests. The skew constant: the pairwise data sets
have CutsInterpolation = F, whose Fortran comparison allows a few ulp (decision D4), so the impulse
sets (compared bitwise) got QD-wrapper points too, including centres with u = 0, where
x = cos(π/3)·v alone. The S2 fallback: no fold output left the cell by more than the margin, so a
test now evaluates points outside the cell directly. The S4 cell of the first group: the test cut
with its own grid differed only at the end of its grid, so a cut with a moved interior knot was
added. The S2 margins are redundant by construction: the z test needs none (the computed distance
is never below the computed height difference), and in the plane either margin alone covers the
rounding; no test case can tell them apart.

## The real data set (H₂/Cu(111) SRP48)

`potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v`, a copy of upstream `potential.dir` (A. Mondal et
al., J. Phys. Chem. A 117, 8770 (2013); upstream input `input.dir/v0j0_WPC_H2Cu111SRP_C3v_01.inp`):
C3v, Å and eV, a = 2.6015 Å, 4 layers, 10 atom cuts, CutsInterpolation = F (the cuts hold V, up to
1950 eV at the wall). `CRPData.6D` lists 31 molecular cuts, 29 used; the two extra ones have no file
and are skipped, as in the Fortran (two warnings). Checked 2026-11-11 against Fortran 5.19:
- **Reading** (`tests/crplib/test_real_data_set.cpp`): every file as gfortran reads it, including the
  trailing comments after the items (`<-- …`, `# …`, `3.4d0  ZCM …`).
- **The library on it**: V equals the data at 3132 cut knots (cut site and orientation, Z below
  ZsoMin) to 6.2e-15 relative; C3v symmetry and atom exchange to 2.3e-14; the gradient matches
  finite differences; S2, S3, S4 bitwise; 118 of the 256 cluster atoms in reach.
- **Against the Fortran library** (reference case `cu111_srp48`, driver on the 5.19 objects; the
  point lists follow the data set's ranges): atom PES bitwise in bohr/Hartree and Å/eV; molecule
  PES, QD wrapper, gas phase and CAMBCORD within D4 (≤ 1.9e-14 Hartree); with libmvec's cos the
  molecule PES is bitwise equal too.
- **Against the Fortran program** (gates 1/2, `tests/reference/crp_cu111_srp48_small`, the upstream
  input on coarse grids with VMax imposed): V_asy 7.7e-16, 6D 1.5e-13 of max(|V|, 1 eV) with and
  without `UseSymmetry` (D4; 1.5e-14 with libmvec's cos), rovibrational states pass.
- **On the upstream production grid** (C++, potential only, `UseSymmetry` + `nPhisOptimal`;
  nX = nY = 20, JMax 18, MJMax 12): 77 of 400 (x, y) points computed, self-check 9.2e-14 eV,
  40 of 64 r points beyond Rmax = 2.3 Å (clamped), no (Z, r) point without data; 82 s on 4
  threads and 1 GB (indicative: measured while another job ran).

## Checking a new data set

1. Run the program in potential-only mode (`"Potential": true`) and read the `CRP DATA SET`
   block: the data ranges, the clamped grid points, the (Z, r) points without data (should be 0)
   and the warnings.
2. With `UseSymmetry` the PES is compared with its symmetry images at 1000 points before use
   (`SymmetryTolerance`); a data set that is not symmetric fails there.
3. For a comparison with the Fortran, run `tools/fortran_reference/run_setup_reference.sh` on a
   case whose `input.json` names the data set (it is copied as `dataset`: the Fortran fails for
   data-set paths longer than 50 characters).

## Appendix: cuts, sites and orientations

Sites in lattice units, skewed (u, v) for C3v and C6v, Cartesian for C4v; θ/φ in degrees. The
"images" are the symmetry-equivalent sites the interpolation also uses, with the rotation of φ.

**C6v molecule (18 cuts)**

| Cuts | Site | Orientations θ/φ, in cut order |
|---|---|---|
| 1–3 | top (0, 0) | 0; 90/0; 90/30 |
| 4–7 | bridge (½, 0) | 0; 90/0; 90/60; 90/90 |
| 8–12 | hollow (⅓, ⅓) | 0; 90/0; 90/30; 45/30; 45/210 |
| 13–18 | t2h (⅙, ⅙) | 0; 90/30; 90/120; 45/30; 45/210; 45/120 |

Images: bridge (½, ½) φ+60 and (0, ½) φ−60; hollow (⅔, ⅔) φ+180; t2h (⅔, ⅙) φ−120,
(⅙, ⅔) φ+120, (⅚, ⅚) φ+180, (⅓, ⅚) φ+60, (⅚, ⅓) φ−60.

**C3v molecule (29 cuts)**

| Cuts | Site | Orientations θ/φ, in cut order | Images |
|---|---|---|---|
| 1–7 | top, bridge | as C6v | as C6v |
| 8–12 | hcp (⅓, ⅓) | as the C6v hollow | — |
| 13–18 | t2h (⅙, ⅙) | 0; 90/30; 90/120; 45/30; 45/120; 45/210 | (⅔, ⅙) φ−120, (⅙, ⅔) φ+120 |
| 19–24 | t2f (⅓, ⅚) | 0; 90/240; 90/330; 45/150; 45/240; 45/330 | (⅚, ⅚) φ+120, (⅚, ⅓) φ−120 |
| 25–29 | fcc (⅔, ⅔) | as hcp | — |

**C4v molecule (16 cuts)**

| Cuts | Site | Orientations θ/φ, in cut order | Images |
|---|---|---|---|
| 1–3 | top (0, 0) | 0; 90/0; 90/45 | — |
| 4–7 | bridge (½, 0) | 0; 90/0; 90/45; 90/90 | (0, ½) φ+90 |
| 8–10 | hollow (½, ½) | 0; 90/0; 90/45 | — |
| 11–16 | t2h (¼, ¼) | 0; 45/225; 45/45; 90/45; 90/135; 45/135 | (¾, ¼) φ−90, (¾, ¾) φ−180, (¼, ¾) φ+90 |

**Atom-surface sites** (`CRPData.3D` takes Cartesian lattice units)
- C6v (7, skewed): (0, 0), (0, ½), (⅓, ⅓), (⅙, ⅙), (⅙, ⅓), (0, ⅙), (0, ⅔).
- C3v (10, skewed): (0, 0), (½, ½), (⅓, ⅓), (⅙, ⅔), (⅙, ½), (0, ⅙), (0, ⅓), (⅓, ⅚), (½, ⅔), (⅔, ⅔).
- C4v (6, Cartesian): top (0, 0), hollow (½, ½), bridge (½, 0), t2h (¼, ¼), b2h (½, ¼), t2b (¼, 0).
