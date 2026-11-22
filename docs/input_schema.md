# Input file (JSON)

`spodvr <input.json>` reads one JSON object. Keys are the Fortran `DVR.inp` field names
(`InputNew.f90`), **case-sensitive**, without renames. Unknown keys are an error, and so is a key
given twice in one object (the Fortran reader stops on a repeated field too). Integer keys
reject non-integer numbers; real keys accept any JSON number. Units are those of the Fortran
input: bohr, eV, amu, atomic units of time, degrees. `//` and `/* */` comments are accepted.
Example: `tests/reference/hardwall_small/input.json`.

## Required keys

These are the fields the Fortran reader pre-registers as compulsory, plus `TolSp`, which has no
Fortran default.

| Key | Type | Fortran | Meaning |
|---|---|---|---|
| `Title` | string | title | |
| `SystemID` | string | systemID | PES; implemented: `HardWallPES` (C2v), `H2Cu111DaiLight`, `H2Cu111Persson` (both C6v), `CRPlibrary` (a CRP data set, `CRPlibDataSet`; C3v, C4v or C6v from the data set), `CRPDistrtSurface` (a C3v or C6v CRP data set plus the SCM correction of a distorted slab, [SCM keys](#scm-keys-systemid-crpdistrtsurface)) (case-insensitive) |
| `NewCalc` | bool | newcalc | new calculation; `false`: restart from `RestartDir` (tStart = end of the earlier run) |
| `ZDoF` `RDoF` `XDoF` `YDoF` `ThetaDoF` `PhiDoF` | bool | | must all be true (6D only) |
| `Mass1` `Mass2` | real | amass1/2 | atom masses (amu) |
| `HomoNuc` | bool | lhomo | homonuclear molecule; must be true with `CRPDistrtSurface` (SCM D13) |
| `ZStart` `dZ` | real | zmin, delz | Z grid start, spacing (shared by both Z grids) |
| `nZ` `nZSpecular` | int | nz, nzsp | scattering / specular Z points (valid FFT sizes, nZSpecular ≥ nZ) |
| `RStart` `dR` `nR` | real, real, int | rmin, delr, nr | r grid; RStart > 0 (more than 1e-12 above 0, the Fortran's tolerance) |
| `nX` `nY` | int | nx, ny | x, y points (valid FFT sizes) |
| `JMax` | int | jmaxw | maximum J of the rotational basis |
| `ProjOp` | int | nproj | must be 2 (0 and 1 are rejected; 1 is a legacy option that is not ported) |
| `tStart` `tEnd` `dtStep` `nStep` | real×3, int | starttime, tottym, dtlan, nstep | propagation |
| `dtAna` `dtAve` | real | dtAnalysis, dtAverages | analysis / averages interval |
| `DynamZSpAna` | bool | dynamicalZSpAna | dynamical analysis surface on the specular grid |
| `Zana` | real | ZAnalysis | analysis line (on the Z grid) |
| `NrOfEs` | int | nenrgy | number of analysis energies, at least 2 (0 and 1 are input errors, analysis D1); −95 (Fortran `Constants_NEN95`): the collision energies k·π/tEnd (atomic units) between EMinWP + EParallel and EMaxWP + EParallel |
| `EMinWP` `EMaxWP` | real | | wave packet energy range (eV) |
| `Z0` | real | z0 | initial wave packet centre |
| `ZOptStart` `ZOptEnd` | real | zsplmi, zsplma | CAP in Z, scattering grid |
| `ROptStart` `ROptEnd` | real | rsplmi, rsplma | CAP in r |
| `ZSpOptStart` `ZSpOptEnd` | real | ssplmi, ssplma | CAP in Z, specular grid |
| `EParallel` `PhiParallel` | real | epar, phipar | parallel energy (eV), azimuth (deg) |
| `writeSMatrix` `writeProbs` `writeRDProbs` | bool | lwrit3, lwri25, lwdiro | output switches: `writeProbs` writes probabilities_rovib.tsv, `writeRDProbs` probabilities_diffraction.tsv; `writeSMatrix` writes fine_analysis/spectrum/s_matrix.tsv, the S-matrix of every state on the fine-analysis grid (Fortran: unit 3). The on-the-fly s_matrix.tsv is always written |
| `TolSp` | real | tolsp | specular-grid tolerance; also the maximum initial norm outside [ZAna, start of the specular absorber] (hard error) |

For each CAP, exactly one of `ZOptEnergy` / `ZOptFactor` (likewise `ROptEnergy` / `ROptFactor`,
`ZSpOptEnergy` / `ZSpOptFactor`) is required.
The energy is in eV, and the factor is the strength in hartree. An energy of 0 means no
absorption (absorber D1: the Fortran silently keeps a default strength instead).

## Optional keys (Fortran defaults)

| Key | Type | Default | Fortran |
|---|---|---|---|
| `MJMax` | int | JMax | mjmaxw. A negative value means JMax. **0 also gives JMax** for the basis (Fortran quirk) |
| `VIn` `JIn` `MJIn` | int | 0 | initial state |
| `onlyAnalysis` | bool | false | analysis-only run (nread = 2): the fine analysis of the BK file in `RestartDir` (format `RestartFormat`); no propagation and no 6D potential. `tEnd` sets the number of BK records used (tEnd / dtAna) |
| `alsoAnalysis` | bool | false | the fine analysis after the propagation, of the BK file the run has written (Fortran replaces the on-the-fly results by the fine ones; here both are written) |
| `Potential` | bool | false | potential-only run (nread = 6) |
| `JAna` `VAna` `DifOrdAna` | int | determined by the program | analysis limits |
| `ZAnaSp` | real | none | fixed analysis line on the specular grid |
| `doFluxAnalysis` | bool | false | doFluxAnalysis; needs the four keys below. The flux analysis (`FluxAnalyzor.f90`) is **not ported** ([flux_analysis_findings.md](flux_analysis_findings.md)): a propagating run (new or restart) with `true` stops with "doFluxAnalysis = true: the flux analysis is not implemented yet" (exit code 3) once the input checks have passed. Analysis-only and potential-only runs accept the keys and ignore them |
| `rFluxAnalysis` `EMinFluxAna` `EMaxFluxAna` `DeltaEFluxAna` | real | none (required with `doFluxAnalysis`) | rFluxAnalysis, EMinFluxAna, EMaxFluxAna, deltaEFluxAna. rFluxAnalysis (bohr): a point of the r grid, ≤ ROptStart. Energies (eV): EMinWP ≤ EMinFluxAna < EMaxFluxAna ≤ EMaxWP and 0 < DeltaEFluxAna ≤ EMaxFluxAna − EMinFluxAna. Two of these rules differ from the Fortran checks: EMaxFluxAna ≤ EMaxWP (flux Q1) and DeltaEFluxAna > 0 (flux Q2); see [fortran_input_bugs.md](fortran_input_bugs.md) |
| `ImposeVMin` `VMin` | bool, real | false, −5.5 eV | clamp the potential from below |
| `ImposeVMax` `VMax` | bool, real | false, 10.0 eV | clamp the potential from above |
| `UseSwitching` `ZminSwitch` `ZmaxSwitch` | bool, real, real | false, 6.0, 8.0 | switch V to the gas phase; zasymp becomes ZmaxSwitch |
| `write2DPES` `x2DPES` `y2DPES` `theta2DPES` `phi2DPES` | bool, real, real, real, real | false, 0.0, 0.0, 90.0, 0.0 | write2DPES: write V(r, Z) of the PES at the fixed x, y (bohr, along the lattice vectors), θ, φ (degrees) to `pes_cut_2d.dat` ([output_files](output_files.md#data-files-space-separated)) right after the PES is set up, then stop (exit code 0), in every run mode, before the 6D potential and the gas phase, as the Fortran. Every r grid point and every scattering-grid Z point (also beyond zasymp); V is the PES itself, without switching and without VMin/VMax. With `write2DPES`, 0 ≤ x2DPES ≤ alat, 0 ≤ y2DPES ≤ blat (the PES's lattice constants) and 0 ≤ θ, φ ≤ 360 (the Fortran's tolerant comparisons); the plot keys are refused (nothing writes them) |
| `UseSymmetry` | bool | false (usePotSymmetry) | compute V only at one (x, y) point per orbit of the PES's surface point group (C2v, C4v, C6v, C3v, C3) and map it to the others; operations that do not fit the grid (nφ, nX ≠ nY) are dropped and reported. Error if none is usable. The Fortran uses this for C6v only. An error with `CRPDistrtSurface` (SCM D2) |
| `nPhisOptimal` | bool | false | with `UseSymmetry`: nφ becomes a multiple of what the whole group needs (C2v 2, C4v 4, C6v and C3v 6, C3 3) |
| `CRPlibDataSet` | string | — | directory of the CRP data set (format: [crp_library.md](crp_library.md)), relative to the working directory; required with `SystemID` `CRPlibrary` and `CRPDistrtSurface`, an error with any other PES. The data set is read and checked at start-up; the setup output shows its ranges and the grid points outside them |
| `SymmetryTolerance` | real | 1e-5 eV | C++ only. With `UseSymmetry` the PES is first compared with its symmetry images at 1000 sample points (x, y, φ anywhere in the cell; θ, Z, r at grid points; images wrapped into the unit cell, VMin/VMax applied); a larger max \|ΔV\| is an error, and so is a non-finite V at a sample point or image (parameter `SystemID`, naming the point) |
| `ZOptLowStart` `ZOptLowEnd` | real | none: off | zsplmilow, zsplmalow (keys of upstream 5.24). An absorbing potential at the low end of the scattering Z grid: the high-end profile mirrored, strongest at Start, 0 at End. On only if both are given, with ZStart ≤ Start < End ≤ Zana. These two bounds are compared exactly, as the Fortran does (plain `<`, `>`). Start > End, Start = End and an incomplete set of the four keys are errors (low-Z D1–D3, [fortran_input_bugs.md](fortran_input_bugs.md#the-low-z-absorbers-of-upstream-524)). It may share the point Z = Zana with the high-end absorber (both are 0 there). The set-up warns when a grid point lies on different sides of Start here and in the Fortran ([inputs/README.md](../inputs/README.md#low-z-absorbers-and-the-grid-point-trap)) |
| `ZOptLowEnergy` `ZOptLowFactor` | real | none (exactly one when the absorber is on) | ZOptLowEnergy (eV, ≥ 0 within 1e-12; 0 means no absorption), vabszblow (strength in hartree) |
| `ZSpOptLowStart` `ZSpOptLowEnd` `ZSpOptLowEnergy` `ZSpOptLowFactor` | real | none: off | ssplmilow, ssplmalow, ZSpOptLowEnergy, vabssblow: the same on the specular grid, also with End ≤ Zana. Where it overlaps the high-end specular absorber (possible when Zana > ZSpOptStart), the low-end value replaces the high-end one, as in the Fortran |
| `DoDump` `FullSMatrix` | bool | false | accepted and ignored: the gas-phase states that `DoDump` writes in the Fortran are always written here (`gas_phase_states.dat`); `FullSMatrix` has no effect in the Fortran either ([fortran_deviations.md](fortran_deviations.md)) |

## SCM keys (SystemID `CRPDistrtSurface`)

The distorted-surface PES of upstream 5.24: the CRP potential of the data set `CRPlibDataSet`
plus the SCM correction ΔV of one snapshot of a distorted slab,
V = CRP(u/s, v/s, Z, r, θ, φ) + ΔV(x, y, Z, r, θ, φ), with (x, y) the Cartesian position of the
lattice coordinates (u, v) and ΔV converted to eV with the Fortran's precise factor (SCM D4). These
keys are an error with any other `SystemID` (SCM D1). With `CRPDistrtSurface`, `CRPlibDataSet`
must be a C3v or C6v data set (the slab is FCC111, SCM D14), `UseSymmetry` false (SCM D2) and
`HomoNuc` true (SCM D13). Example files: `potentials/H2Cu111SRP_C3v.dir/params.dat` and
`examplePosfile.dat` (for the SRP48 data set in that directory). The SCM decisions are listed in
[fortran_deviations.md](fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface).

| Key | Type | Default | Meaning |
|---|---|---|---|
| `DistortedPESMode` | int | 1 | DistPESMode: 1 = SCM, the only mode; any other value is an error, as in the Fortran (SCM D7) |
| `SCMParFilepath` | string | none: required (SCM D12) | ParamFile: the pair-function parameter file (`VERSION 2`, then 18 values per layer for layers 1 and 2), relative to the working directory |
| `PosFilepath` | string | none: required (SCM D12) | SlabPosFile: the slab position file (FCC111; versions 3 and 2, SCM D5), relative to the working directory. Its ideal atoms must lie on the CRP lattice: layer 1 on the sites, with one at the origin, at z = 0; with a C3v data set layer 2 under the site (1/3, 1/3), one interlayer spacing lower. More than 1e-3 bohr off is an error (SCM D9) |
| `PosFileSlabNum` | int | 1 | SurfNumb: the snapshot used, 1 … NumSlabConfigs of the file (0 or less is an error, SCM D6) |
| `EnableSurfExpan` | bool | true | EnabSCMExpan: rescale the slab to the CRP lattice. The lattice constants of the grid become the slab's a_T (line 3 of the position file; it must then be positive, SCM D6), the CRP is evaluated at (u/s, v/s) with s = a_T/a₀, and ΔV subtracts the ideal lattice in the CRP's frame. With false the grid has the CRP's lattice constant a₀ |
| `SCMCutoff` | real | 16.0 bohr | SCMCutoff: version 3 uses the atoms of layers 1–2 whose distorted position lies closer than this to the origin; must be positive |

## Checkpoint / restart keys (C++ only)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `RestartDir` | string | — | output directory of the run to continue; required if `NewCalc` is false. It must differ from `--output-dir`. For `onlyAnalysis`: the directory with the BK file to analyse and the checkpoint of the same run, whose analysis surface gives the specular analysis line that the initial state's S-matrix phases refer to (fortran_deviations.md, analysis D2) (required; it may equal `--output-dir`). The files the run reads from it (`RestartFormat`) are checked before the set-up: a missing one is an input error (exit code 2) |
| `RestartFormat` | string | `"native"` | `"native"`: checkpoint.bin + bk_coefficients.bin; `"fortran"`: fort.8 (or fort.9) + fort.22 (or fort.21) written by the Fortran code or by this program. `onlyAnalysis` reads the BK file and, of the checkpoint, only the analysis surface |
| `CheckpointFormat` | string | `"native"` | files written for a later restart: `"native"`, `"fortran"` (fort.9, fort.21) or `"both"`. bk_coefficients.bin is always written |

The fine analysis checks the BK file in the same way (the Fortran code ignores the header of
fort.22). Its records must lie every dtAna from t = 0.

A restart must use the same input as the earlier run except `NewCalc`, `tStart`, `tEnd`, `nStep`
(and output switches); the program checks the grids, masses, time steps, wave packet and
absorbers against the checkpoint and names any parameter that differs. This includes the number
of φ points, which `nPhisOptimal` with `UseSymmetry` can change, and the low-Z absorbers: a
checkpoint or BK file written with a low-Z absorber that the input does not have (or the other way
round) is refused with a message naming that absorber. Files written without low-Z absorbers have
the same contents as before these keys existed.
The checkpoint records the `SystemID` but not the CRP data set: a restart with SystemID
`CRPlibrary` must use the same `CRPlibDataSet` (not checked). With `CRPDistrtSurface` the native
files also record `PosFileSlabNum`, `SCMCutoff`, `EnableSurfExpan`, the slab's a_T, atom count and
position statistics (rms displacement, mean shift of layer 1, highest atom) and the 36 values of the
parameter file, and a restart or analysis-only run with others is refused naming them. Positions are
compared through these statistics, not atom by atom; the data set is not checked. The Fortran-format files (fort.8/fort.9, fort.21/fort.22)
hold no parameters, so a restart from them checks none of this.

## Plot data keys (C++ only)

These keys write extra data for the plot tool (`tools/plots/README.md`). They are output switches, not
restart parameters: a restart may set them differently from the run it continues.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `writeDensities` | bool | false | write `densities/` (output_files.md, "Series directories"): the probability per (Z, r) cell of the scattering and of the specular wave at tStart and every dtAve. Propagating runs only (`NewCalc` true, or a restart); refused with `Potential` or `onlyAnalysis`. It adds no computation to the time loop: the records are the densities the averages already compute, written once per averages row |
| `writePotential` | bool | false | write `potential/` (output_files.md, "potential/") once at set-up: the potential as the propagation uses it (switching and the VMin/VMax limits included) on the grid of Z < zasymp, its minimum over the DVR orientations at every grid point, and V at every orientation at the `PlotSites`. New, restart and potential-only runs (`Potential`: it is written before the run stops); refused with `onlyAnalysis`, which computes no 6D potential |
| `dtFrames` | real | none (no frames) | write `frames/` (output_files.md, "frames/"): at tStart and every dtFrames a.u., reductions of the whole wave packet ψ_full = ψ_scat + s2s·χ (the probability per (Z, r, x, y) cell and per (Z, rotational state), and χ beyond the scattering grid). Propagating runs only; refused with `Potential`, `onlyAnalysis` and `write2DPES`. dtFrames > 0, at least dtAna and a whole multiple of it (and dtAna a whole multiple of dtStep); tStart and tEnd whole multiples of dtFrames, so that the frames of a restart continue those of the run before. Each frame is one pass over the Z planes after that analysis step (on the GPU, a copy of each plane to the host); the PLOT DATA line projects its cost (about 0.03 % of the CPU propagation and 0.7 % of the GPU's at dtFrames = 25 dtAna on the thesis DZ grid), the end of the run prints the measured one |
| `FrameExtras` | list | none | optional datasets of `frames/` (output_files.md, "frames/"), each a list entry: `"angular"` (the probability per DVR node (θ, φ) and Z), `"diffraction"` (per parallel momentum k∥ + n b₁ + m b₂ and Z), `"vibration"` (per gas-phase state φ_vJ(r) and Z), `"rotation_moments"` (⟨P₂(cos θ)⟩ and ⟨J(J+1)/(2μr²)⟩ per Z), `"momentum"` (per k_r, and per k_Z of the specular wave). Needs `dtFrames`; known names only, each once. Each is computed in the frame's pass over the Z planes on the host (a transform per plane: about the cost of one propagation stage, more for `vibration`); the end of the run prints the measured seconds per frame of each |
| `NormBudget` | logical | false | write `norm_budget/` (output_files.md, "norm_budget/"): per analysis window, what each absorber removed, estimated from the cell probabilities at the window's start and end. Propagating runs only; refused with `Potential`, `onlyAnalysis` and `write2DPES`; dtAna must be a whole multiple of dtStep. It computes the block densities of ψ after every analysis step instead of every dtAve (the averages of those times use the same ones); the PLOT DATA line projects the cost (about 1.5 % of the CPU propagation and 1.9 % of the GPU's on the thesis DZ grid at dtAna = 8 dtStep) |
| `EnergyResolvedSites` | list | none | names of `PlotSites` at whose (x, y) grid points the run sums the wave packet over time at each energy of `EnergyResolvedEnergies` (output_files.md, "energy_resolved/"): s(Z, r, rot; E) = Σ_t ψ_scat(t) e^{i E_tot t} Δt_ana, the time-to-energy transform of the on-the-fly S matrix applied to ψ, and the same sum of the specular wave per energy; the plot tool makes the stationary scattering state ρ_E(Z, r) of it. Needs `EnergyResolvedEnergies` (and that key needs this one); propagating runs only (refused with `Potential`, `onlyAnalysis`, `write2DPES`); each name once; at most 256 (site, energy) pairs. A restart continues the sums of RestartDir's `energy_resolved/` and needs the same sites and energies. Cost: one copy of a site's (Z, r, rotation) column of ψ per site and analysis step (on the GPU a download) and a complex multiply-add per value and energy on the host; the PLOT DATA line projects it and warns above 10 % of an analysis step, the end of the run prints the measured time |
| `EnergyResolvedEnergies` | list | none | collision energies (eV) of the energy-resolved sums, each inside the packet's range EMinWP + EParallel .. EMaxWP + EParallel (as the analysis energies; E_tot = E + E(vin, jin), as the S matrix's); needs `EnergyResolvedSites` |
| `PlotSites` | list | none | named (x, y) grid points for the plot data, used by `writePotential` and `dtFrames` (the frames' manifest records them for the plots): a list of at most 16 objects, each with the fields `name` (a string, unique), `x` and `y` (fractions of the first and the second lattice vector, in [0, 1)), e.g. `[{"name": "top", "x": 0, "y": 0}, {"name": "bridge", "x": 0.5, "y": 0}]`. Each site must lie on a grid point (x·nX and y·nY whole numbers to 1e-9; the error names the nearest point). Refused when no plot key uses the sites; recorded in run_info.json (`plot_data`) |

The set-up prints a "PLOT DATA" section with the projected size of each plot-data directory. A run
whose plot data would not fit in the free space of the output directory's file system stops before the
propagation (exit code 2), naming the key whose directory is the largest.

## Not supported

These keys give a "not supported" error. Logical keys are accepted when set to `false`. The list
is: `SystemMix`, `SystemID2`, `SystemMixCoeff`, `UseRCutOff`, `RCutOff`, `PESMixingAlfa`.

## Validation

All checks of `InputNew.f90` ReadInput/ValidateInput are reproduced, with the same tolerant real
comparisons (absolute tolerance 1e-12). So are the grid checks of `Setup.f90` and `setupdvr.f90`
and the nφ checks of `transform.f90` (which are hard errors here). Every check has a unit test,
most of them in `tests/test_input.cpp`; the initial norm outside the allowed region (`TolSp`) in
`tests/test_initial_wave_packet.cpp`, the `RestartDir` checks in `tests/test_restart.cpp`.
Deviations are listed in [fortran_deviations.md](fortran_deviations.md).
Two Fortran flux analysis checks are wrong and are fixed here (flux Q1, flux Q2), and so are
the silent cases of the low-Z absorber keys (low-Z D1–D3), the default strength the Fortran
keeps for an absorber energy of 0 (absorber D1) and NrOfEs = 0 and 1 (analysis D1):
[fortran_input_bugs.md](fortran_input_bugs.md) describes them, with demonstration runs of the
Fortran program.
