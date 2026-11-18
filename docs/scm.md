# SCM: the distorted-surface potential

The static corrugation model (SCM) of upstream SPO-DVR 5.24, SystemID `CRPDistrtSurface`: the 6D
CRP potential of an ideal surface plus a pairwise correction ΔV for the displaced metal atoms of
one frozen snapshot of a hot surface. One run uses one snapshot; the observable at the surface
temperature is the average over many runs, taken outside the program ([Ensembles](#ensembles-by-hand)).
The C++ port reproduces the Fortran: ΔV bitwise, gates 1–4 in the case `scm_cu111_srp48_small` on
the CPU and the GPU ([validation.md](validation.md)); its decisions (deviations and the Fortran
behaviour kept) are SCM D1–D15 in
[fortran_deviations.md](fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface).
Whether the model as coded is what its authors intended is open:
[scm_open_questions.md](scm_open_questions.md) explains the questions Q1–Q13. The repository holds
one example snapshot and no production data: enough for method tests, not for physics.

**Where the numbers come from.** [run]: the set-up report of a potential-only run of
[inputs/scm_example](../inputs/scm_example/README.md) (its input plus `"Potential": true`, run from
the repository root; 2026-11-15, build of 5e44ba8); [run −0.55]: the same with `ZStart` −0.55 bohr,
`nZ` 120, `nZSpecular` 210; a test name: the value that test asserts; [Q]: scm_open_questions.md
(measured with the C++ code on the samplings defined there).

## The model

```
V(X, Y, Z, r, θ, φ) = V_CRP(X/s, Y/s, Z, r, θ, φ)
                    + Σ_{H atoms a = 1, 2} Σ_{Cu atoms j} [ f_L(|h_a − R_j^distorted|) − f_L(|h′_a − R_j^ideal|) ]
```

- **Frozen surface** (the sudden approximation): H₂ crosses the barrier region fast compared with
  the surface vibrations, so each snapshot is held fixed during a run. Energy exchange with the
  lattice (recoil, phonons) is not described.
- **ΔV** moves every Cu atom j of the top two layers from its ideal site to its snapshot position:
  it removes that atom's H–Cu pair energy at the site and adds it at the displaced position. Version
  3 files keep the atoms whose displaced position lies within `SCMCutoff` (default 16 bohr) of the
  origin: 37 + 36 atoms in the example [run]. An undistorted slab without rescaling gives ΔV = +0
  exactly (`ScmCorrection.UndistortedSlabWithoutRescalingGivesExactlyPlusZero`).
- **The pair function** f_L (layer L's line of `params.dat`) is a modified Rydberg form switched
  smoothly to a constant between p7 and p8 (7.44 and 10.46 bohr: P7, P8 of the example file): an
  atom farther than p8 from both H atoms adds nothing. Its parameters are linear in the bond length
  r, clamped to [P17, P18] = [0.725, 1.847] bohr as coded; the published model clamps in Å (Q1).
- **Lattice rescaling** (`EnableSurfExpan`, default on): the grid's cell takes the slab's lattice
  constant a_T, V_CRP is evaluated at (X/s, Y/s) with s = a_T/a₀, and the subtracted term uses the
  molecule's centre scaled by 1/s (h′). Example: a_T 4.88967, a₀ 4.91615 bohr, s = 0.994612 [run],
  a 0.54 % lateral compression of the CRP, as published (Q4).
- **What it captures:** a raised atom adds local repulsion, a lowered or shifted one softens or
  moves the corrugation. **What it misses:** every non-pairwise effect of the distortion, atoms
  below layer 2, any change of the H–H interaction.

## The input

The SCM keys keep the Fortran names ([input_schema.md](input_schema.md#scm-keys-systemid-crpdistrtsurface)):

```json
"SystemID": "CRPDistrtSurface",
"CRPlibDataSet": "potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v",
"HomoNuc": true,
"DistortedPESMode": 1,
"SCMParFilepath": "potentials/H2Cu111SRP_C3v.dir/params.dat",
"PosFilepath": "potentials/H2Cu111SRP_C3v.dir/examplePosfile.dat",
"PosFileSlabNum": 1,
"EnableSurfExpan": true,
"SCMCutoff": 16.0
```

The data set must be C3v or C6v (SCM D14), `HomoNuc` true (SCM D13) and `UseSymmetry` false
(SCM D2); both file paths are required (SCM D12) and relative to the working directory. A complete
input: [inputs/scm_example/input.json](../inputs/scm_example/input.json); how to choose its grids:
[inputs/README.md](../inputs/README.md#12-scm-distorted-surface-systemid-crpdistrtsurface).

## The data files

- **Parameter file** (`SCMParFilepath`, e.g. `params.dat`): line 1 `VERSION 2` (upper case, from
  column 1; SCM D8), then one record of 18 numbers per layer, layer 1 then layer 2. P1–P8 are the
  intercepts of p1 … p8 (reference distance, decay, the four polynomial coefficients, switch start
  and end), P9–P16 their slopes per bohr of bond length, P17 < P18 the clamp range of r. Lengths in
  bohr, energies in hartree. The example's two lines are identical (Q5).
- **Position file** (`PosFilepath`): line 1 a tag and the structure (`FCC111`, the only one); line 2
  nx, ny, nz (atoms per row, rows per layer, layers ≥ 2), the number of snapshots, the version (3 or
  2; 1 is an error, SCM D5) and the atoms per snapshot of version 2; line 3 a_T in bohr. From line 4
  each snapshot has one record per atom (nx·ny·nz in version 3, the atoms per snapshot in version
  2): displaced x, y, z, ideal x, y, z (Cartesian, bohr, the ideal top layer at z = 0 with a site at
  the origin) and the layer (1 = top). The ideal atoms must lie on the CRP lattice (SCM D9). The
  example `examplePosfile.dat`: a 7 × 7 × 7 slab, one snapshot, version 3.

## Ensembles, by hand

SPO-DVR runs one snapshot; there is no ensemble tool yet. The example position file holds a single
snapshot, so this procedure needs a production file (Q9).

1. Write N inputs that differ only in `PosFileSlabNum` (1 … N of one position file) and run each
   with its own `--output-dir`. Paths are relative to the working directory: run all inputs from
   one working directory, or use absolute paths.
2. Check the SCM identity `scm` in every `analysis.json` ([output_files.md](output_files.md#results-tab-separated-tables)):
   `slab_number` runs through 1 … N; `crp_data_set`, `parameter_file`, `position_file`, `cutoff`,
   `rescaling`, `parameters` and `slab_lattice_constant` are equal, and so are `energies` and
   `analysis` (the same grid and channels); the snapshot statistics (`slab_atoms`,
   `rms_displacement`, `top_layer_shift`, `highest_atom`) differ, so two runs with equal statistics
   most likely used the same snapshot.
3. Average the linear columns row by row with equal weights (canonical sampling; the published
   runs average their 104 slabs, Q6): `P_scattered`, `P_reaction`, `P_v…` of
   `scattering_summary.tsv`, and `P` of `probabilities_rovib.tsv`, `probabilities_diffraction.tsv` and `s_matrix.tsv`. Give the standard
   error of the mean, s/√N with s the standard deviation over the N runs (N − 1 in its denominator).
4. Do not average `Re_S`, `Im_S` or `phase`, nor the per-scattered-molecule means `E_vib_mean`,
   `E_rot_mean`, `E_dif_mean` with equal weights: each is Σ P·E / P_scattered of one run i, so the
   ensemble value is Σᵢ P_scattered,i·Eᵢ / Σᵢ P_scattered,i.
5. **Which observables average across cells (Q2):** each run repeats one distorted cell
   periodically, so its diffraction comes from an artificial surface with jumps at the cell edges.
   SCM was introduced for reaction probabilities (their broadening with the surface temperature);
   whether averaged diffraction and specular probabilities mean anything is open.
6. **Re-centred cells (Q11):** shifting one snapshot onto another top atom gives more samples, but
   they are not independent: do not count them in N for the standard error.

## The SCM set-up report

The block "SCM DISTORTED SURFACE" follows the CRP data set in the set-up output
([output_files.md](output_files.md#results-tab-separated-tables)). Its lines, with the example's values [run]:

| Line | What it shows | Example |
|---|---|---|
| files, slab tag, structure, version, snapshots, snapshot used, cut-off | what was read | `Multiseedv3_New`, FCC111, version 3, 1 snapshot, 16 bohr |
| Atoms used; lattice constants; rescaling | the atoms ΔV sums over; a_T, a₀ and s | 37 + 36; 4.88967, 4.91615 bohr; s = 0.994612 |
| Displacement rms, largest; mean of layer 1; highest atom | displaced − ideal over the atoms used; the largest displaced z | 0.639, 1.211 bohr; (−0.174, 0.176, −0.083) bohr; 0.936 bohr |
| Ideal atoms off the CRP lattice (SCM D9) | the largest distance of an ideal atom from its CRP site: a warning above 1e-5, an error above 1e-3 bohr | 1.17e-6 bohr |
| Gas-phase points with V ≠ V_CRP (SCM D3) | the V_asy points (with switching also the switch's V_gas points) where ΔV changes V; any one stops the run (exit code 2) | 0 of 96 (Z = 11.9, 12.28 bohr) |
| Largest \|V − V_CRP\| on Z = z_asymptotic (SCM D3) | over the grid's (x, y, θ, φ) at r < `ROptStart`; a warning above 1 meV | 4.8e-5 eV (r = 4.0 bohr, θ = 25°) |
| Z planes below the CRP's Zmin (SCM D15) | the lowest V − V_CRP there at Rmin ≤ r ≤ Rmax; a warning | none (`ZStart` 0.5 ≥ Zmin 0.472 bohr); [run −0.55]: 7 planes, −335 eV at Z = −0.55, r = 2.0 bohr |
| r outside Rmin .. Rmax (SCM D15) | the lowest V − V_CRP where the CRP clamps r, every Z plane, the grid's first φ only; whether it lies below Rmin or inside the r absorber; no warning | −232 eV at Z = 1.55, r = 5.8 bohr, θ = 57.4°, inside the r absorber |
| Jump table (SCM D10) | the largest \|V(t·a, 0) − V(t·a, a)\| and \|V(0, t·a) − V(a, t·a)\| at θ = 90°, 12 points t per edge, 12 φ; "accessible": both sides below min V_asy + `EMaxWP` + 0.5 eV | (Z, r) = (2.25, 1.95) bohr: 1.10 eV, 0.61 eV accessible; (4, 1.4) bohr: 0.13 eV |

## Limits and open questions

Limits of the model (the first from the literature, the others inferred from the model): no
energy transfer to the lattice; pairwise and two layers only, so it is as good as its fit; not meant
for slow, trapping-mediated channels; lattice zero-point motion only if the MD sampling has it. "The
barrier slice" below is (Z, r) = (2.25, 1.95) bohr, near the SRP48 minimum barrier.
[scm_open_questions.md](scm_open_questions.md#the-questions-in-brief) explains the questions to the
authors; the C++ reproduces the Fortran until they answer.

| Question | Size in the example | Source | C++ |
|---|---|---|---|
| [Q1](scm_open_questions.md#q1-the-units-of-the-bond-length-clamp): the clamp, in Å in [Spiering et al. 2018](https://doi.org/10.1063/1.5058271), in bohr in the code (and, as we read them, in the published elbow cuts) | the lowest barrier of the origin cell: 0.501 eV as coded, 0.639 eV with the published clamp (ideal surface 0.625 eV) | [Q] | as coded (bohr) |
| [Q2](scm_open_questions.md#q2-the-jump-at-the-cell-edges): the jump at the cell edges | 0.700 eV on the barrier slice, both sides below the gas phase + 1.5 eV (Q2: 0.68 eV with 30° φ steps; the C++ takes 15°, SCM D10) | `ScmDiagnostics.ScmD10ExampleJumpNearTheBarrier` | reproduced, reported (D10) |
| [Q3](scm_open_questions.md#q3-deep-holes-where-the-crp-clamps-z-or-r): holes where the CRP clamps Z or r | V − V_CRP −335 eV below Zmin, −232 eV beyond Rmax (Q3 gives V − E_gas on its own sampling) | [run −0.55]; [run] | reproduced, reported (D15) |
| [Q12](scm_open_questions.md#q12-hd-on-the-crp-and-scm-potentials): HD | both atoms misplaced by r/6 along the bond: V changes by up to 0.52 eV on the barrier slice | [Q] | HD rejected (D13) |
| [Q7, Q8, Q10](scm_open_questions.md#technical-points-q7-q8-q10): technical | ΔV leaves V_asy unchanged (0 of 96 points) | [run] | version 1 rejected (D5), V_asy checked (D3), no `.pes` cache (D11), `UseSymmetry` rejected (D2) |
| [Q4, Q5, Q6, Q11, Q13](scm_open_questions.md#what-the-papers-answer-q4-q5-q6-q11-q13): answered in the papers | s = 0.994612 (the EAM's 925 K lattice over the CRP's 0 K one); identical layer parameters; one cell per slab, centred per trajectory; cut-off 16 bohr | [run]; [Q] | as coded |
| [Q9](scm_open_questions.md#data-we-would-like-to-have-q9): data | one snapshot, no production slabs | the files | the example only |

**In your own runs, until the authors answer:**
- Q1: results use the clamp in bohr, as the Fortran does. On the example snapshot the published clamp
  moves the lowest barrier by +0.14 eV; allow for an uncertainty of that order in barrier-sensitive
  results.
- Q2: treat per-snapshot diffraction and specular probabilities as unconverged in nX, nY.
- Q3: choose the Z grid and the r absorber as [inputs/README.md](../inputs/README.md#12-scm-distorted-surface-systemid-crpdistrtsurface)
  says (`ZStart` ≥ Zmin, or `ImposeVMin`; r beyond Rmax inside the r absorber).
- Q11: cells re-centred on another top-layer atom are correlated samples (item 6 of
  [Ensembles](#ensembles-by-hand)).

## Links and references

- [scm_open_questions.md](scm_open_questions.md): the open questions, with figures and the
  literature; [fortran_deviations.md](fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface):
  SCM D1–D15; [validation.md](validation.md): the gates; [upstream_sync.md](upstream_sync.md): the
  port status.
- [input_schema.md](input_schema.md#scm-keys-systemid-crpdistrtsurface), [output_files.md](output_files.md),
  [inputs/README.md](../inputs/README.md#12-scm-distorted-surface-systemid-crpdistrtsurface),
  [inputs/scm_example](../inputs/scm_example/README.md).
- Wijzenbroek & Somers, J. Chem. Phys. 137, 054703 (2012), doi:10.1063/1.4738956: the original SCM.
- Smits, Litjens & Somers, J. Chem. Phys. 156, 214706 (2022), doi:10.1063/5.0094985: the first
  quantum dynamics with SCM.
