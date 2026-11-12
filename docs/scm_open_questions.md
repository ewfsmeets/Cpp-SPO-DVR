# SCM distorted-surface potential: open questions

Upstream SPO-DVR 5.24 adds a potential for a thermally distorted surface: the **static corrugation
model** (SCM), SystemID `CRPDistrtSurface`. It is the CRP 6D potential of an ideal surface plus a
pairwise correction for the displaced metal atoms of one snapshot of a hot surface. Reading the
Fortran closely and measuring it raised questions that only its authors can settle. This document
explains each one for a reader who knows SPO-DVR and the CRP potentials but not SCM: what the code
does, why it is questionable, how large the effect is, what to do in your own runs meanwhile, and
what the C++ code does. It is also the basis of our message to the authors.

**Status.** The C++ code does not implement SCM yet: it rejects `SystemID: CRPDistrtSurface` as an
unknown system. The planned port reproduces the Fortran exactly (it is the reference code), with
three exceptions that reject an input instead of guessing: HD (Q12), `UseSymmetry` (Q10) and
position-file version 1 (Q7). It adds reports and warnings where noted below. If the authors answer
a question with a change, the change becomes an explicit, documented option with a test; the
default stays the Fortran's behaviour.

**Where the numbers come from.** Every number describes **one example snapshot**,
`potentials/H2Cu111SRP_C3v.dir/examplePosfile.dat` with `params.dat` from the same directory and the
SRP48 CRP data set `H2Cu111SRP_C3v`; production snapshots are not available (Q6). Each number is
marked with its source:
- [P] measured 2026-11-11 with a scratch C++ prototype of the 5.24 SCM correction whose ΔV is
  bitwise equal to the Fortran's (4048 of 4048 test points, against the 5.24 objects built with
  `-O3`), combined with the C++ CRP library;
- [T] computed 2026-11-12 with a Python transcription of the same Fortran formulas (ΔV only,
  accurate to well below 1 meV; it reproduces the [P] values of ΔV), as in
  [make_figures.py](scm_open_questions/make_figures.py);
- [lit] from the literature or general knowledge, not checked here.

The numbers will be re-measured with the C++ port once it exists. Fortran references name the file
and routine in upstream `sources.dir/SOURCES.dir/DistortedSurfacePES.dir/`.

## Contents

- [SCM in one page](#scm-in-one-page) and [notation](#notation)
- [The questions by priority](#the-questions-by-priority)
- [Q2. The jump at the cell edges, and what a 1×1 cell can describe](#q2-the-jump-at-the-cell-edges-and-what-a-11-cell-can-describe)
- [Q4. The lattice rescaling](#q4-the-lattice-rescaling)
- [Q1. The units of the bond-length clamp](#q1-the-units-of-the-bond-length-clamp)
- [Q11. Re-centring a snapshot, and its drift](#q11-re-centring-a-snapshot-and-its-drift)
- [Q3. Deep holes where the CRP clamps Z or r](#q3-deep-holes-where-the-crp-clamps-z-or-r)
- [Q6, Q5, Q9. The data: production snapshots, parameters, provenance](#q6-q5-q9-the-data-production-snapshots-parameters-provenance)
- [Q12. HD](#q12-hd)
- [Technical questions](#technical-questions) (Q7, Q8, Q10, Q13)
- [Literature](#literature)

## SCM in one page

**The problem.** A normal SPO-DVR run puts H₂ on an ideal, rigid surface: every Cu atom sits at
its lattice site. The experiments SCM is compared with are at a surface temperature T_s; the
state-resolved H₂/Cu(111) data at 925 K were obtained from associative desorption and converted to
reaction probabilities by detailed balance [lit]. At T_s every Cu atom is displaced from its site
and the lattice is expanded. Measured reaction-probability curves broaden with T_s; a rigid, ideal
potential cannot broaden them.

**The idea: freeze the surface (the sudden approximation).** H₂ crosses the barrier region in
roughly 10–20 fs, short compared with the period of the surface vibrations (about 140 fs for the Cu
Debye frequency, θ_D ≈ 343 K [lit]). So: take one snapshot of the hot surface from a classical
molecular-dynamics (MD) run, freeze it, run an ordinary 6D wave-packet calculation on the potential
of that frozen surface, and repeat for many snapshots. The observable at T_s is the average over
snapshots, P(E; T_s) = (1/N) Σᵢ Pᵢ(E). SPO-DVR does one snapshot per run; the averaging happens
outside it. What this misses: energy exchange with the lattice (recoil, phonons) and surface
motion during the collision.

**The potential** (`DistortedSurfaceWrapper.f90 GetCRPDistSurfPotential`):

```
V_SCM(X, Y, Z, r, θ, φ) = V_CRP(X/s, Y/s, Z, r, θ, φ)
                        + Σ_{H atoms a = 1, 2} Σ_{Cu atoms j} [ f_L(|h_a − R_j^displaced|) − f_L(|h′_a − R_j^ideal|) ]
```

- The second term, **ΔV**, moves each Cu atom from its ideal site to its snapshot position: it
  removes that atom's H–Cu interaction at the ideal site and adds it at the displaced one.
- **The ideal positions** R_j^ideal (from the position file) lie exactly on the CRP set's lattice,
  nearest-neighbour distance a₀ = 4.916154 bohr (2.6015 Å) [T]. **The displaced positions** come
  from MD on a lattice with the snapshot's own constant a_T (here 4.88967 bohr, per the file
  header; a fit to the displaced top layers cannot tell a_T from a₀ at this noise level).
- **Lattice rescaling** (Q4): s = a_T/a₀ (here 0.99461). The CRP term is evaluated at X/s, Y/s, so
  its ideal cell is mapped onto the snapshot's cell; in the subtracted term, the molecule's centre
  is scaled laterally by 1/s (bond vector unchanged), giving h′_a. So V_CRP and the subtracted
  term describe the same ideal surface, and ΔV swaps in the real atoms.
- For an undistorted snapshot on the a₀ lattice without rescaling, ΔV is exactly zero.
- **Which atoms count** (`SCM_IO.f90 ReadSCMSlab`, file version 3): the Cu atoms of layers 1 and 2
  whose displaced position lies within `SCMCutoff` (default 16 bohr) of the origin: 73 atoms here
  (37 + 36).

**The pair function** f_L(R) (`SCM_Potential.f90 Switch_Rydberg_Potential`) is a modified Rydberg
form, −(p3 + p4·q + p5·q² + p6·q³)·exp(−p2·q) with q = R − p1, switched smoothly (cosine) to a
constant C between p7 ≈ 7.44 and p8 ≈ 10.46 bohr. C cancels in ΔV, so an atom farther than p8 from
both H atoms (in both of its positions) contributes exactly 0. **Each parameter depends on the bond
length** (`SCM_Types.f90 Setup_3Body_Rydberg_Par`): p_k = P_k + P_(8+k)·t, with t = r clamped to
[P17, P18] = [0.725, 1.847] (Q1). This makes f an effective three-body (H–H–Cu) term. Values at
r = 1.4 bohr [T]: f(1) = 45.0 eV, f(2) = 3.13 eV, f(2.5) = 0.82 eV, zero at 5.5 bohr, a well of
−34.5 meV at 6.75 bohr, C = −4.2 meV.

![The SCM pair function f(R) of the example params.dat for three values of the clamped bond length t. As coded, every bond longer than 1.847 bohr uses the orange curve, which stays repulsive at the 3.3-bohr H–Cu distance of H in the fcc hollow; with the clamp read in Å, the green curve (r = 3.49 bohr) is attractive there (−0.26 eV).](scm_open_questions/pair_function.png)

**What it captures:** a raised atom adds local repulsion, so the barrier above it rises; a lowered
or shifted atom softens or moves the local corrugation. **What it misses:** every non-pairwise
electronic effect of the distortion, atoms below layer 2, and any change of the H–H interaction.

**The example snapshot** [T]: a 7 × 7 × 7 slab, one configuration, file version 3. Top-layer
displacements about the layer mean: rms 0.35, 0.25 and 0.45 bohr in x, y and z (0.62 bohr in 3D);
the deepest atom is 1.16 bohr below its site, the highest 1.24 bohr above (outside the cut-off
sphere; the highest kept atom is 0.94 bohr above the surface plane). The deeper layers have
0.24–0.35 bohr per component, close to the classical Debye estimate for Cu at 925 K, about
0.25 bohr [T, with θ_D from lit]: consistent with a 925 K snapshot.

**Size of ΔV in the origin cell** (12 × 12 points of the cell, r = 1.4 bohr, θ = 90°, φ = 0) [T]:

| Z (bohr) | ΔV min … max (eV) | cell mean (eV) |
|---|---|---|
| 2.5 | −0.70 … −0.05 | −0.24 |
| 4 | −0.17 … −0.01 | −0.12 |
| 6 | −0.008 … +0.029 | |
| 10 | ≈ 10⁻⁴ | |

At Z ≤ 4 bohr ΔV is negative everywhere in this cell (at Z = 6 bohr, 117 of the 144 points are
positive) [T]. A plausible reason: the cell lies in a local depression of the snapshot, with its
four corner atoms 0.01–0.57 bohr below their sites and a neighbour 1.16 bohr below, while the layer
means are only −0.026 (layer 1) and −0.004 bohr (layer 2) [T]. Another origin gives other signs. ΔV is exactly 0 once both H atoms are more than p8 above every kept atom (the highest is at
z = 0.94 bohr, so from z = 11.40 bohr up).

**One SCM run** uses one snapshot (`PosFileSlabNum`); an ensemble is N runs that differ only in
the snapshot number. The C++ input will use the same key names (documented with the port).

### Notation

| Symbol | Meaning |
|---|---|
| a₀ | nearest-neighbour distance of the CRP set's ideal lattice: 4.916154 bohr (bulk lattice constant √2·a₀ = 3.679 Å) |
| a_T | nearest-neighbour distance of the snapshot's MD lattice, from the position file: 4.88967 bohr (bulk 3.659 Å) |
| s | a_T/a₀, the lateral rescaling factor (here 0.99461) |
| u, v | lattice coordinates along the two 60° cell axes; X = u + v/2, Y = (√3/2)·v |
| L | the layer of a Cu atom (1 = top); f_L uses parameter line L of `params.dat` |
| P_1 … P_18 | the 18 numbers of a line of `params.dat` (lengths in bohr, energies in hartree; for P17/P18 see Q1) |
| p_1 … p_8 | the pair parameters at a given bond length: p_k = P_k + P_(8+k)·t |
| t | the bond length r clamped to [P17, P18] |
| "lowest V on the slice" | the minimum of V over 12 × 12 sites, θ = 0, 22.5, …, 90° and 12 φ at fixed (Z, r) = (2.25, 1.95) bohr, near the SRP48 minimum barrier; a coarse sample, not a saddle-point search |

## The questions by priority

The labels Q1–Q13 are stable names (used in our message to the authors), not ranks.

| Priority | Question | Size (example snapshot) | Effect on ensemble-averaged observables | In your own runs until answered | C++ port |
|---|---|---|---|---|---|
| 1 | Q2: the jump at the cell edges; what a 1×1 cell can describe | steps of 0.3–0.9 eV near the barrier | **systematic** in diffraction and specular scattering; reaction: unknown | treat per-snapshot diffraction as unconverged | reproduces it; reports the jump |
| 2 | Q4: the lattice rescaling | +14 meV mean ΔV near the barrier; CRP values unchanged (mapped); open: is the CRP set the 0 K or the 925 K lattice? | **systematic** | run as coded | as coded |
| 3 | Q1: is the bond-length clamp in bohr or Å? | lowest V on the slice +0.034 eV (limits in Å) or −0.29 eV (r dependence in Å); larger beyond the barrier | **systematic** | allow ≈ 0.03 eV (limits in Å) to 0.3 eV (r dependence in Å) of systematic uncertainty | as coded (bohr) |
| 4 | Q11: may a snapshot be re-centred; remove its drift? | lowest V on the slice: sd 0.061 eV over the 49 cells of one snapshot; drift 12–52 meV | statistics, and a small systematic shift | do not count re-centred cells as independent | the example only |
| 5 | Q3: holes where the CRP clamps Z or r | V down to 255 eV below the gas phase below Z = 0.47 bohr, 357 eV for r > 4.35 bohr | none unless amplitude reaches them | keep the Z grid above 0.47 bohr; set `ImposeVMin` | reproduces them; warns |
| 6 | Q12: HD (all these PESes, not only SCM) | V off by up to 0.4–0.5 eV at accessible tilted geometries near the barrier height | systematic (HD only) | avoid HD | rejects HD with SCM |
| 7 | Q6, Q5, Q9: the production data and provenance | — | not needed: the port is validated against the Fortran on the example; the data would only allow comparison with published results | use the example for method tests only | the example only |

## Q2. The jump at the cell edges, and what a 1×1 cell can describe

**What the code does.** The x, y grid covers one surface unit cell, [0, a_T)² in (u, v) with a
top-layer atom at the origin, and the FFTs make it periodic. ΔV is computed from the actual atoms
around the origin (`SCM_IO.f90 ReadSCMSlab` keeps those within the cut-off sphere around the
origin), and the periodic grid then repeats this one distorted cell everywhere. A 1×1 cell holds
one atom per layer, so the only distortion that is periodic in it is a rigid shift of a whole
layer: every real snapshot breaks the periodicity.

**Why that makes a jump.** Near one edge, an H atom feels the actual displaced atoms beyond that
edge. Its periodic copy just inside the opposite edge feels the atoms on the other side of the
origin cell, which are displaced differently. So the potential steps where the two copies meet.

![ΔV along a line across the cell edge v = 0 at u = (2/3)·a, Z = 2.25, r = 1.95 bohr, θ = 90°, φ = 0: within 0.01 eV of the largest jump along this edge at θ = 90° (24 points, 12 φ). Dashed: the snapshot itself, smooth across the edge. Solid: what the run uses, the origin cell repeated periodically, with a 0.68 eV step at the edge.](scm_open_questions/edge_jump.png)

**How large.** The jump is |V(u, 0) − V(u, a_T)| along the edge v = 0 (and the same along u = 0):
the two values that the periodic grid joins. The CRP term is periodic, so the jump is that of ΔV.
- Near the barrier, (Z, r) = (2.25, 1.95) bohr, θ = 90°, 12 φ, 12 points per edge, counting only
  points where both sides lie below the gas phase + 1.5 eV (the region a molecule with up to about
  1 eV collision energy can reach): 0.64 eV on average over the 49 possible cell origins of the
  snapshot (the maximum per origin; range 0.32–0.92 eV); 0.68 eV for the example's own origin [P].
- For comparison, at Z = 2.5 bohr, r = 1.4 bohr: the jump is 0.54 eV at most and 0.20 eV on
  average along both edges (16 points per edge, 6 φ), while ΔV changes by 0.011 eV (median) and 0.094 eV (max) between
  neighbouring points a/16 apart inside the cell [T]. The step is a real discontinuity, not a steep
  slope.
- It fades with height (r = 1.4 bohr, maximum along the edges): 0.54 eV at Z = 2.5, 0.13 eV at
  Z = 4, 0.035 eV at Z = 6 and 2 × 10⁻⁴ eV at Z = 10 bohr [P].
- Two of the three bridge sites of the cell lie on its edges, and the bridge is where the
  H₂/Cu(111) barrier is lowest. This snapshot's lowest V on the slice lies on an edge, at
  (0.42, 0)·a [P].

**What it changes.**
1. The barrier at an edge site has two values, one on each side.
2. A step in a periodic potential has Fourier components that fall off only as 1/k, and in SPO-DVR
   those components couple to the diffraction channels. The step therefore feeds high diffraction
   orders in every snapshot. This adds |V_G|² ≥ 0, so, to first order and unless the step is
   correlated with the snapshot's own Fourier components, it does **not** average out over an
   ensemble. We do not know how much of the computed diffraction or specular attenuation comes from
   it; this is a hypothesis to check, not a finding.
3. The x, y grid converges slowly (the 1/k tail), so per-snapshot diffraction probabilities may not
   converge with nX, nY. Not tested.

**The underlying question: what a 1×1 cell can describe.** With one distorted cell repeated, the
surface is perfectly periodic, so it scatters only into Bragg (diffraction) channels. Thermal
attenuation of the specular peak in experiment goes mostly into diffuse, non-Bragg scattering
(the Debye–Waller picture [lit]), which a 1×1 cell cannot produce. Which observables is an SCM run
meant to give: reaction probabilities, or also diffraction and specular attenuation?

**Possible answers.** (a) It is intended: the ensemble average is taken as the observable (then:
was the size of the artefact checked?); (b) choose the cell origin so that no bridge site lies on
an edge; (c) make ΔV periodic over the cell (this changes the model: only a rigid layer shift is
periodic); (d) use a larger cell (a supercell: more grid points, more diffraction channels).

**A check we can do ourselves.** Compare the Fourier spectrum of V(u, v) at barrier height for the
ideal and the SCM surface; then one reduced-grid run with ΔV blended to be periodic across the
edges, comparing reaction and specular probabilities.

**In your own runs:** treat per-snapshot diffraction and specular probabilities as unconverged in
nX, nY. **The C++ port:** reproduces the jump (it is the model as implemented) and reports its size
at a few heights in the set-up output.

## Q4. The lattice rescaling

**What the code does** (`EnableSurfExpan`, default on; `SCM_Wrapper.f90 SetSCMCRPParameters`,
`DistortedSurfaceWrapper.f90 GetCRPDistSurfPotential`, `SCM_Potential.f90 SCM_SurfInteract_3Body`).
With rescaling on, SPO-DVR's unit cell becomes a_T, the CRP potential is evaluated at (u/s, v/s),
and the molecule's centre in the subtracted term is scaled by 1/s laterally; Z, r, the bond vector
and the layer spacings are not scaled. With `EnableSurfExpan: false`, s = 1: the cell stays at a₀,
the CRP is evaluated as is, and ΔV compares the displaced atoms (on the a_T lattice) with the ideal
a₀ lattice directly, so it also contains the lattice mismatch, growing with distance from the
origin.

**The numbers.**
- a_T = 4.88967 bohr is a bulk lattice constant of 3.659 Å. The experimental Cu lattice constant is
  3.615 Å at room temperature and about 1.1–1.2 % larger at 925 K, about 3.656 Å [lit, estimate]:
  a_T agrees with it to about 0.1 %. So the MD box was probably set to the experimental 925 K
  lattice.
- a₀ = 4.916154 bohr (3.679 Å) is the CRP set's lattice: the SRP48 PES of Mondal et al. 2013, a
  paper on the effect of the thermal lattice expansion at 925 K [lit].
- The code divides the first by the second: s = 0.99461, a 0.54 % lateral **compression** of the
  CRP landscape.
- The rescaling alone, for a perfect slab at a_T: ΔV of +14 meV on average and at most +25 meV on
  the slice ((Z, r) = (2.25, 1.95) bohr, θ = 90°, 12 × 12 sites, 6 φ), about ±1 meV at Z = 4 bohr [T].
  This is the pair model's estimate of what compressing the lattice by 0.54 % does at fixed height
  and bond length (atoms closer, so more repulsion): a property of the model, not a coding
  inconsistency. It is probably what Spiering et al. 2018 call "corrected stretching" [inference;
  that paper was not consulted].
- The rescaling changes no CRP value at corresponding lattice positions: V_CRP(u/s, v/s) on the a_T
  cell is the a₀ cell's value at the same fractional position. It compresses the CRP landscape
  laterally by 0.54 % (the cell, the corrugation lengths and the diffraction channels follow), so
  the CRP's classical barrier heights are unchanged; barrier heights move only through ΔV.

**The questions.**
- Is `H2Cu111SRP_C3v` the ideal-lattice (DFT-relaxed) PES or the expanded-lattice (925 K) PES of
  Mondal et al. 2013? If it is the expanded one, rescaling it to a_T compresses an already expanded
  landscape back; if it is the ideal one, the rescaling maps the DFT lattice onto the experimental
  925 K lattice, and the 0.54 % "compression" only says that DFT's 0 K lattice is larger than
  experiment's 925 K lattice.
- Is the rescaling meant relative to the CRP's DFT lattice (as coded), or should the experimental
  relative expansion be applied?
- Should Z and the interlayer spacings scale too?
- Are the ideal positions in version-3 files always on the CRP lattice? The code relies on it
  (the subtracted term is consistent with V_CRP only then) and never checks it. The C++ port checks
  it: the ideal top layer on the CRP lattice, and layer 2 at the CRP's interlayer spacing under its
  hcp site (1/3, 1/3) (both true for the example [T]).

**Effect.** Systematic: every snapshot is rescaled the same way. The open point is the absolute
level: whether the CRP barrier is that of the 0 K or of the 925 K lattice (Mondal et al. 2013, and
Smits et al. 2024 on SRP48 and EAM lattice constants, may answer this; neither was consulted).
**In your own runs:** as coded. **The C++ port:** as coded.

## Q1. The units of the bond-length clamp

**What the code does.** `SCM_Types.f90 Setup_3Body_Rydberg_Par` clamps the bond length r to
[P17, P18] = [0.725, 1.847] before the eight pair parameters are formed from it. r is in bohr, as
everywhere in SPO-DVR, so the clamp acts in bohr.

**Why it looks wrong.**
- In bohr, the lower limit (0.725 bohr = 0.38 Å) acts only where the H₂ wave function is
  negligible.
- The upper limit, 1.847 bohr, lies just before the barrier (r ≈ 1.9–2.0 bohr): the r dependence
  stops where the bond starts to break, where it should matter most.
- Read as Å, the range is [1.37, 3.49] bohr: from the equilibrium bond length, through the barrier,
  to dissociation, the range one would choose for such a fit.
- **The strongest hint** (figure above) [T]: at the H–Cu distance of H chemisorbed in the fcc
  hollow, R ≈ 3.3 bohr (1.73 Å [lit]), f stays repulsive for every bond length as coded
  (f(3.3) ≈ +0.19 eV at t = 1.847). With the limits in Å, f(3.3) becomes attractive as the bond
  breaks (≈ −0.26 eV at r = 3.49 bohr), as an H–Cu bond should after dissociation.

**Two Å readings.** (b) only the limits are in Å: t = clamp(r, 0.725 Å, 1.847 Å) with r in bohr,
slopes per bohr; (c) the whole r dependence is in Å: t = clamp(0.529177·r, 0.725, 1.847) with r in
bohr, intercepts unchanged. ΔV(reading) − ΔV(as coded), over 12 × 12 sites and 6 φ at θ = 90° [T]:

| (Z, r) in bohr | (b) limits in Å: min … max (mean), eV | (c) r dependence in Å: min … max (mean), eV |
|---|---|---|
| (4.0, 1.4) | 0 | −0.055 … +0.019 (−0.018) |
| (2.5, 1.4) | 0 | −0.21 … +0.001 (−0.12) |
| (2.2, 1.95), near the barrier | −0.012 … +0.040 (+0.015) | −0.31 … +0.081 (−0.12) |
| (2.0, 2.3) | −0.16 … +0.20 (+0.039) | −0.27 … +0.20 (−0.058) |
| (1.75, 2.6) | −0.67 … +0.37 (−0.031) | −0.22 … +0.39 (+0.013) |
| (1.75, 3.5) | −2.24 … +0.95 (−0.089) | 0 (as coded and (c) both clamp to t = 1.847) |

The signs depend on the site, so the barrier shift is not simply the mean. With the full V
(CRP + ΔV), **the lowest V on the slice** is 0.488 eV as coded, **0.522 eV** with reading (b) (the
barrier rises by 0.034 eV: curves shift to higher energy) and **0.201 eV** with reading (c) (it falls
by 0.29 eV: curves shift to lower energy) [T, with the C++ CRP library]. Reading (c) changes ΔV by
−0.12 eV on average already at the H₂ bond length (Z = 2.5 bohr), which a check of the fit would
probably have caught; reading (b) is the more plausible alternative. Beyond the barrier, reading (b)
changes ΔV by −0.67 … +0.37 eV at (1.75, 2.6) bohr and −2.2 … +0.95 eV at (1.75, 3.5) bohr: the exit
channel, where recrossing is decided, so reaction probabilities may change more than the barrier
shift alone suggests.

**Effect.** Systematic, the same for every snapshot: a barrier shift of +0.03 eV (b) or −0.3 eV
(c), against chemical accuracy of 0.04 eV. Which reading produced the published
results is part of Q9. **Settled by:** the fit range of r in Spiering, Wijzenbroek & Somers 2018 (not
consulted here), or the authors' answer. **In your own runs:** as coded; allow a systematic
uncertainty of ≈ 0.03 eV if the limits are in Å (b), up to 0.3 eV if the whole r dependence is (c). **The C++ port:** as coded (bohr).

## Q11. Re-centring a snapshot, and its drift

**One snapshot is not representative.** Taking each of the example's 49 top-layer atoms in turn as
the cell origin (with the atoms within 22 bohr of the new origin, from periodic images of the 7 × 7
slab: 7·a_T for the displaced and 7·a₀ for the ideal positions), the lowest V on the slice is
0.585 eV above the gas phase on average, with a spread (sd) of 0.061 eV and a range of
0.470–0.707 eV; the ideal surface (the CRP term alone, whose values the rescaling does not change)
gives 0.625 eV, the example's own origin 0.488 eV, the second lowest of the 49 [P]. That the cell
average lies below the ideal value is not by itself a systematic barrier lowering: a minimum over
sites drops for any zero-mean disorder.

**The question.** The 49 cells tile the MD supercell, so averaging over them is the natural
per-snapshot average, and it comes for free from one MD run. Are they valid samples, or too
correlated to count as independent? (The sd is a spread, not a standard error.) The version-3 reader
has no option to re-centre; a user would have to translate the position file.

**Drift.** The snapshot is displaced on average against its ideal lattice, by different amounts per
layer [T] (bohr):

| Atoms | mean x | mean y |
|---|---|---|
| layer 1 (top) | −0.144 | +0.197 |
| layer 2 | −0.066 | +0.158 |
| layer 7 (bottom) | −0.009 | −0.087 |
| whole slab | −0.096 | +0.079 |

Part of the top layer's shift is drift of the whole slab (MD), the rest a layer-wise fluctuation. A
rigid translation must not change observables, but the pair model breaks that invariance. Change
of the cell-mean ΔV near the barrier ((Z, r) = (2.25, 1.95) bohr, θ = 90°, 12 × 12 sites, 6 φ) [T]:

| Mean removed | lateral only | 3D (incl. the mean height) |
|---|---|---|
| whole slab | +12 meV | +16 meV |
| top layer | +23 meV | +52 meV |

The top layer's mean height of −0.026 bohr alone is worth about 30 meV: the cell-mean ΔV changes by
about 1 eV per bohr of mean top-layer height. Removing the whole-slab drift is justified; removing
the top-layer mean also removes a thermal mode. The version-1 reader's refocus option
`SCMSlabRefocusMethod = 2` (hard-coded off) subtracts the top layer's mean position in all three
components (+52 meV here); version 3 removes nothing. Which is intended?

**In your own runs:** do not count re-centred cells as independent samples. **The C++ port:** the
example file only; re-centring and drift removal would be explicit options after an answer.

## Q3. Deep holes where the CRP clamps Z or r

**What the code does.** The CRP potential is defined only for Z ≥ Zmin and r ≤ Rmax (for SRP48,
Zmin = 0.472 bohr and Rmax = 4.35 bohr = 2.3 Å); outside, the CRP library uses the value at the
limit (it clamps Z and r). The SCM correction uses the true Z and r (`SCM_Wrapper.f90
GetSCMPotential` receives them unclamped).

**Why that makes holes.** Whenever an H atom comes within 1–2 bohr of Cu atoms, the pair function is
extremely repulsive (45 eV at 1 bohr), and ΔV is the difference of two such large numbers: the
repulsion from an atom's displaced position minus that from its ideal site. If the atom moved away
from the H atom, much more is subtracted than added, and ΔV is very negative. The CRP part, frozen at
its limit, does not compensate. This happens in two regions (V relative to the gas phase) [P, T]:

| Region | Where | Lowest V_SCM |
|---|---|---|
| below the CRP's Z range | (Z, r) = (0, 1.4) / (−0.5, 1.4) / (−0.5, 3.0) / (−1, 1.4) / (−1, 2.0) bohr | −29 / −99 / −255 / −2.8 / −122 eV |
| beyond the CRP's r range, tilted molecule | Z = 1.5–2.5 bohr, r = 6–9.8 bohr, θ = 45–68° (the lower H atom at the height of the top layer) | −357 eV (where V_CRP = +7.3 eV) |

On the ideal surface the lowest V in the first region is 1.2–3.1 eV above the gas phase, and in the
second region V_CRP lies between +0.6 and +774 eV [P, T]. For Z ≥ 1 bohr and r = 1–4.5 bohr there are
no spurious holes: the SCM minimum is only 55 meV below the ideal surface [P].

**What it changes.** Near the barrier, nothing. But the first region sits behind a thin wall (at
Z = 0.5 bohr the potential is only about 1.4 eV above the gas phase [P]), and the second lies in the
r range of dissociation, usually inside the r absorber. Amplitude trapped in a hole never reaches
the scattering analysis, so 1 − ΣR counts it as reaction, and the grid and time step are not
converged for a well of hundreds of eV, so the result there is inaccurate. Whether a run is affected depends on its grid: we have not
seen amplitude there in a run.

**Possible answers.** (a) Use the clamped Z and r in ΔV too, consistent with the CRP; (b) require
grids to stay within the CRP's ranges; (c) it does not matter for the published grids (which grids
were used?).

**In your own runs:** start the Z grid at or above Zmin, and set `ImposeVMin` at least 0.1 eV below
the ideal surface's lowest V on your grid (SCM legitimately lowers minima by up to 55 meV at
Z ≥ 1 bohr). **The C++ port:** reproduces the holes, warns when an SCM
grid reaches below Zmin or beyond Rmax, and reports the lowest V_SCM − V_ideal there.

## Q6, Q5, Q9. The data: production snapshots, parameters, provenance

**Why ask, and why last.** The repository holds one example snapshot and one `params.dat`. The C++
port needs nothing more: it is validated against the Fortran on the example. Reproducing a published
production run is not a goal, so the production data would only allow a comparison with published
curves (which could settle Q1 and Q4). Two of the questions stand on their own: the identical layer
parameters (Q5), and which code produced the published results (Q9).

**Q6. The production snapshots.** How many snapshots per T_s, at which T_s, with which weights
(equal weights for canonical sampling)? Is the example a 925 K snapshot from an embedded-atom-method
(EAM) MD run (the file tag `Multiseedv3_New` and the default `SCMType = 'EAMSCM'` suggest so; its
displacements fit 925 K, see the primer)? Could the production position files (e.g. the "Dist"
ensembles and the constrained one-, three- and five-atom "At1/3/5" sets of Smits, Litjens & Somers
2022), the production `params.dat`, the CRP data-set version and the `DVR.inp` be shared?

**Q5. Two layers, one parameter set.** Only Cu atoms of the top two layers enter ΔV (hard-coded:
`UsedLayers = 2` in `SCM_Types.f90`), each layer with its own line in `params.dat`. In the example
both lines are identical. A second-layer atom is partly screened by the top layer, so a fit might
give it weaker parameters. Are identical lines intended, and is the example `params.dat` the
published parameter set or a test file?

**Q9. Provenance.** Which code version and compiler options produced the published QD-SCM results?
Are reference outputs (probabilities of one snapshot with its input) available? One such output
would let us check the whole chain, and would settle Q1 by itself.

**Cost of an ensemble.** One 6D run per snapshot. A production-size grid takes hours of GPU time
per run (the C++ code's thesis benchmark grid: about 1.5–1.8 h for 10 000 steps on an RTX 4070 Ti),
so 100 snapshots are about a week and 1000 are months.

## Q12. HD

**What the code does.** SPO-DVR propagates the molecule's centre of mass and bond vector: the
kinetic operators use the total mass for X, Y, Z and the reduced mass for r (`InputNew.f90`,
`totmass`, `rdmass`). `HomoNuc` changes only the angular grid and the symmetrisation of the
potential (`angulargrid.f90 InitAngularGrid`, `setupdvr.f90`, `Potential.f90`). We found no
routine that moves the centre of mass to the bond midpoint before the CRP or SCM potential is
evaluated: SCM places the atoms at ±r/2 from the given point (`SCM_Wrapper.f90 ConvertCoMToXYZ`),
and so does the CRP library [code, read 2026-11-12; the other PES modules were not checked]. For HD
the centre of mass is r/6 from the bond midpoint (0.23 bohr at r = 1.4 bohr), so an HD run on these
potentials evaluates them with both atoms misplaced by r/6 along the bond.

**How large** [T, with the C++ CRP library], on the slice (Z, r) = (2.25, 1.95) bohr, comparing V
at the given point with V at the point shifted by ±r/6 along the bond: at points below the gas
phase + 1.5 eV, V changes by up to 0.42 eV on the ideal surface and 0.52 eV with SCM. The lowest V
on the slice changes by only +5 meV (+1 meV with SCM): it belongs to a parallel molecule, for which
the shift is lateral and moves the site, not the energy. The orientation dependence near the
barrier, and so rotational and steric effects, are what is affected.

**The question.** Is the shift applied somewhere we have not found (PES data, inputs), or should
every PES evaluation for a heteronuclear molecule shift by the mass ratio? This concerns all HD runs
on these potentials, not only SCM, and published HD results computed with this code path.

**In your own runs:** avoid HD until this is answered. **The C++ port:** it reproduces the Fortran
for the potentials already ported (HD is accepted there, as in the Fortran); for SCM, new code, it
rejects HD until the question is answered.

## Technical questions

No physics decision; short answers suffice.

- **Q7. File versions.** The C++ port reads position-file versions 3 (the example) and 2 (probably
  the constrained one-, three- and five-atom runs). Is version 1 still needed? Its reader has known
  bugs: an `INTEGER` bulk layer spacing (`BulkZDist`), `uCelly` used for `uCellx`, two hard-coded
  layers. The example's version-3 file carries the first bug's fingerprint: its ideal layers 4–7
  are exactly 4.000 bohr apart (√(2/3)·a₀ = 4.014 bohr), after three CRP spacings [T]. That is
  harmless for ΔV (layers 1–2 only), but which tool generated the file?
- **Q8. The asymptotic potential and the `.pes` cache.** V_asy(r) is evaluated at
  (0, 0, Z_asymptotic, r, π/2, 0) and includes ΔV there, except in analysis-only runs
  (`nread = 2`), where the SCM set-up is skipped and ΔV is zero. Should V_asy be the pure CRP? For
  the example, ΔV is exactly 0 at those points: every kept atom lies at least p8 + 0.5 bohr below
  Z_asymptotic = 11.9 bohr [P]; a snapshot (or a re-centring, Q11) that keeps an atom raised a
  further 0.5 bohr would break this (the example's highest atom, 1.24 bohr, would leave a margin of
  0.2 bohr), and then `nread = 2` runs would use other channel energies than the propagation. Also: the SCM
  `.pes` cache is never read in 5.24 (only set-up time is lost), because its file name is built
  before `ParameterVersion` is set and so contains version 0; its reader has an inverted return
  value (`ReadDistPESSetup`, cancelled by a second negation) and opens a unit that is never
  assigned (`TempParFileUnit` in `CheckPESVariables_SCM`); and the cache check compares only set-up
  scalars, not the contents of `params.dat` or the positions, so a working cache could return a
  stale potential. Is this known?
- **Q10. Symmetry.** The point group comes from the CRP set. The Fortran switches symmetry off for
  C3v; with a C6v set and `UseSymmetry` it would symmetrise the distorted V wrongly. Should symmetry
  be forced off with SCM? (The C++ port rejects `UseSymmetry` with SCM.)
- **Q13. The cut-off (a confirmation).** The cut-off sphere is centred at the cell corner (the
  origin), uses the displaced 3D distance, and has radius 16 bohr. Against the full 7 × 7 × 7 slab
  plus its periodic images, leaving out the rest changes ΔV by at most 2.4 meV (rms ≤ 0.2 meV) for
  r ≤ 2.6 bohr at θ = 90°, 3 meV at r = 4.5, 4 meV at r = 6 and 23 meV at r = 9.8 bohr (inside the
  r absorber) [P]. We take 16 bohr as adequate unless told otherwise.

## Literature

- Wijzenbroek & Somers, J. Chem. Phys. 137, 054703 (2012), doi:10.1063/1.4738956: the original SCM.
- Mondal et al., J. Phys. Chem. A 117, 8770 (2013), doi:10.1021/jp4042183: the thermal lattice
  expansion at 925 K, and the SRP48 PES of the CRP set `H2Cu111SRP_C3v` (its `CRPData.info`).
- Spiering, Wijzenbroek & Somers, J. Chem. Phys. 149, 234702 (2018), doi:10.1063/1.5058271:
  effective three-body terms and corrected stretching.
- Smits, Litjens & Somers, J. Chem. Phys. 156, 214706 (2022), doi:10.1063/5.0094985: the first
  quantum dynamics with SCM.
- Smits & Somers, J. Chem. Phys. 157, 134704 (2022), doi:10.1063/5.0112036, and J. Chem. Phys.
  158, 014704 (2023), doi:10.1063/5.0134817: H₂ on Cu(111) at 925 K.
- Smits et al., J. Chem. Phys. 161, 154702 (2024), doi:10.1063/5.0231559: effective Hartree
  potential against the sudden approximation; SRP48 and EAM lattice constants and thermal expansion.
