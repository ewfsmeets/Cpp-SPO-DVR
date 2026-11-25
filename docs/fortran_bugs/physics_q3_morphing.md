# physics Q3: PES morphing changes V where nothing is morphed

A bug of the Fortran SPO-DVR code (upstream r1727, git `3d00331`, and earlier versions as stated below), found in an
audit of the active build for errors that make a computed quantity wrong. The candidate was read by an auditor, checked
by an independent reviewer who saw the code but not the auditor's conclusion, and run where possible. Line numbers refer to
`sources.dir/SOURCES.dir/` at `3d00331`. The short entry is in [fortran_deviations.md](../fortran_deviations.md) under
the same ID; the other bugs of the audit are listed in [README.md](README.md).

**Status:** open in the Fortran (5.20 to r1727); not yet reported upstream (2026-11-25). Not ported to the C++.

## The demonstration runs

- **Executable:** the reference install `r1727-3d00331` (upstream git `3d00331` built with `O3FLAGS = -O3`;
  [tools/fortran_reference/README.md](../../tools/fortran_reference/README.md#the-reference-install)).
- **Inputs:** the `DVR.inp` of a reference case in `tests/reference/` with the keys of each case changed or appended
  (and its data files, as `tools/fortran_reference/reference_config.sh` `stage_inputs` stages them).
- **How:** each run in a new scratch directory through `reference_config.sh` (`run_clean`: a clean environment, 4
  OpenMP threads).

## The Fortran code

When the run directory holds `pesmorphparameters.dat`, `Potential.f90 MorphPotential` (lines 1553–1876, called at
line 720 after the PES is computed and before the switch and the `.pes` cache) changes chosen expansion coefficients of
V. Each line of the file names a coefficient (n, m: the diffraction indices in x and y; J, mJ: the spherical harmonic)
and Gaussian weights in Z and r; per (Z, r) point the routine transforms V(θ, φ, x, y) to (J, mJ, n, m) with the
program's own angular transform (`DVRToFBRGenComplexTransform`, basis J ≤ JMax, |mJ| ≤ MJMax = JMax) and FFTs in x
and y, applies

```fortran
                  irot = reverseRottable(J(iCount),mJ(icount))
                  ixy = xy_table(n(icount),m(icount))
                  vfbr_red(irot,ixy,threadNum) = vfbr_red(irot,ixy,threadNum)*(1+cons(icount)*RWeight*ZWeight)+eps(icount)*RWeight*ZWeight
```

to the eight sign combinations (±n, ±m, ±mJ; lines 1724–1758), and transforms back (lines 1764–1785).

## What goes wrong

- **(a) The round trip is a projection.** The angular grid has (JMax + 1)·nφ points and the basis (JMax + 1)²
  functions, so DVR → FBR → DVR removes every part of V outside J ≤ JMax, |mJ| ≤ JMax, at every (Z, r), whatever the
  file asks. The propagation uses V on the full grid, so this changes the dynamics. Near the repulsive wall the
  projection spreads its large values over all orientations.
- **(b) Zero indices compound `cons`.** For n = 0, m = 0 or mJ = 0 several of the eight updates hit the same
  coefficient (`xy_table(0,m) = xy_table(-0,m)`): 2, 4 or 8 times. The additive `eps` part is then counted k times,
  which the authors' own model expects (upstream `tools.dir/morphing`, a factor 2 per zero index), but the factor
  becomes (1 + cons·W)^k instead of 1 + cons·W.
- **(c) `UseSymmetry`:** the loops run over nx·ny (x, y) points (lines 1670–1700, 1770–1785), but with symmetry the
  potential holds only the representatives, so they read and write beyond each (Z, r) block and beyond the array.
  Already noted from the code in [upstream_sync.md](../upstream_sync.md); not run.

## Demonstration

`cu111_dailight_small` (JMax = MJMax = 4) as a potential-only run, without a morph file (the base), and with each of
these files (format: count, then `n m J mJ maxZ wdthZ maxR wdthR cons eps`):

```text
1                                         1
0 0 0 0 4.0 0.5 1.4 0.2 0.0 0.0           0 0 0 0 4.0 0.5 1.4 0.2 0.05 0.0
(zero: changes nothing)                   (iso: the isotropic coefficient, cons 0.05)
```

Each run prints `* Morphing the potential...` and exits with status 0. The stored V (the `.pes` file of each run):

- **(a)** the zero file changes V at 13 439 996 of the 13 440 000 grid values. At r = 1.1–1.7 bohr: max |ΔV| 0.08 meV
  for Z ≥ 6, 0.53 meV for Z ≥ 5, 4.9 meV for Z ≥ 4, 41 meV for Z ≥ 3 (barrier region) and 292 meV for Z ≥ 2 bohr.
- **(b)** at the centre of the Gaussians (Z = 4.0, r = 1.4 bohr, weight W = 1.5915) the iso file shifts V by the same
  amount at all 1920 (x, y, θ, φ) points (spread 3e-15 of the shift), relative to the zero file: −3.775473 eV on a mean
  of −4.467249 eV, a factor 1.845145. That is (1 + 0.05·W)^8 = 1.845145; one application would give 1.0796.

## The C++ behaviour

Not ported: the C++ never opens `pesmorphparameters.dat` and uses the unmorphed potential
([upstream_sync.md](../upstream_sync.md)).

## Proposed fix

- (a) Apply the change as a difference: transform V, compute the new coefficients, and add the back-transformed
  *change* to the original grid values, so that coefficients left alone keep V unchanged.
- (b) Apply each update once per distinct (irot, ixy) pair, e.g. by skipping the sign combinations that repeat an
  index; keep the k-fold `eps` if that is the intended convention, and document it.
- (c) Refuse morphing with `UseSymmetry`, or morph before the potential is reduced to the representatives.
