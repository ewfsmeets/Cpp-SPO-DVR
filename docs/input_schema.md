# Input file (JSON)

`spodvr <input.json>` reads one JSON object. Keys are the Fortran `DVR.inp` field names
(`InputNew.f90`), **case-sensitive**, without renames. Unknown keys are an error. Integer keys
reject non-integer numbers; real keys accept any JSON number. Units are those of the Fortran
input: bohr, eV, amu, atomic units of time, degrees. `//` and `/* */` comments are accepted.
Example: `tests/reference/hardwall_small/input.json`.

## Required keys
These are the fields the Fortran reader pre-registers as compulsory, plus `TolSp`, which has no
Fortran default.

| Key | Type | Fortran | Meaning |
|---|---|---|---|
| `Title` | string | title | |
| `SystemID` | string | systemID | PES; implemented: `HardWallPES` (C2v), `H2Cu111DaiLight`, `H2Cu111Persson` (both C6v), `CRPlibrary` (a CRP data set, `CRPlibDataSet`; C3v, C4v or C6v from the data set) (case-insensitive) |
| `NewCalc` | bool | newcalc | new calculation; `false`: restart from `RestartDir` (tStart = end of the earlier run) |
| `ZDoF` `RDoF` `XDoF` `YDoF` `ThetaDoF` `PhiDoF` | bool | | must all be true (6D only) |
| `Mass1` `Mass2` | real | amass1/2 | atom masses (amu) |
| `HomoNuc` | bool | lhomo | homonuclear molecule |
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
| `NrOfEs` | int | nenrgy | number of analysis energies |
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
| `doFluxAnalysis` | bool | false | doFluxAnalysis; needs the four keys below. The flux analysis (`FluxAnalyzor.f90`) is **not implemented yet**: a propagating run (new or restart) with `true` stops with "doFluxAnalysis = true: the flux analysis is not implemented yet" (exit code 3) once the input checks have passed. Analysis-only and potential-only runs accept the keys and ignore them |
| `rFluxAnalysis` `EMinFluxAna` `EMaxFluxAna` `DeltaEFluxAna` | real | none (required with `doFluxAnalysis`) | rFluxAnalysis, EMinFluxAna, EMaxFluxAna, deltaEFluxAna. rFluxAnalysis (bohr): a point of the r grid, ≤ ROptStart. Energies (eV): EMinWP ≤ EMinFluxAna < EMaxFluxAna ≤ EMaxWP and 0 < DeltaEFluxAna ≤ EMaxFluxAna − EMinFluxAna. Two of these rules differ from the Fortran checks: EMaxFluxAna ≤ EMaxWP (flux Q1) and DeltaEFluxAna > 0 (flux Q2); see [fortran_input_bugs.md](fortran_input_bugs.md). Memory during the time loop: 32 bytes (φ and φ′) per flux energy, rotational function, (x, y) point and Z point of the stored 6D potential (the scattering-grid points with Z < zasymp − 1e-12, Fortran nzpot), plus up to 8 slots of 32 bytes (ψ and ∂ψ/∂r at r_c) per rotational function, (x, y) point and such Z point ([performance.md](performance.md#memory-and-numa)) |
| `ImposeVMin` `VMin` | bool, real | false, −5.5 eV | clamp the potential from below |
| `ImposeVMax` `VMax` | bool, real | false, 10.0 eV | clamp the potential from above |
| `UseSwitching` `ZminSwitch` `ZmaxSwitch` | bool, real, real | false, 6.0, 8.0 | switch V to the gas phase; zasymp becomes ZmaxSwitch |
| `UseSymmetry` | bool | false (usePotSymmetry) | compute V only at one (x, y) point per orbit of the PES's surface point group (C2v, C4v, C6v, C3v, C3) and map it to the others; operations that do not fit the grid (nφ, nX ≠ nY) are dropped and reported. Error if none is usable. The Fortran uses this for C6v only |
| `nPhisOptimal` | bool | false | with `UseSymmetry`: nφ becomes a multiple of what the whole group needs (C2v 2, C4v 4, C6v and C3v 6, C3 3) |
| `CRPlibDataSet` | string | — | directory of the CRP data set (format: [crp_library.md](crp_library.md)), relative to the working directory; required with `SystemID` `CRPlibrary`, an error with any other PES. The data set is read and checked at start-up; the setup output shows its ranges and the grid points outside them |
| `SymmetryTolerance` | real | 1e-5 eV | C++ only. With `UseSymmetry` the PES is first compared with its symmetry images at 1000 sample points (x, y, φ anywhere in the cell; θ, Z, r at grid points; images wrapped into the unit cell, VMin/VMax applied); a larger max \|ΔV\| is an error, and so is a non-finite V at a sample point or image (parameter `SystemID`, naming the point) |
| `ZOptLowStart` `ZOptLowEnd` | real | none: off | zsplmilow, zsplmalow (keys of upstream 5.24). An absorbing potential at the low end of the scattering Z grid: the high-end profile mirrored, strongest at Start, 0 at End. On only if both are given, with ZStart ≤ Start < End ≤ Zana. These two bounds are compared exactly, as the Fortran does (plain `<`, `>`). Start > End, Start = End and an incomplete set of the four keys are errors (low-Z D1–D3, [fortran_input_bugs.md](fortran_input_bugs.md#the-low-z-absorbers-of-upstream-524)). It may share the point Z = Zana with the high-end absorber (both are 0 there). The set-up warns when a grid point lies on different sides of Start here and in the Fortran ([inputs/README.md](../inputs/README.md#low-z-absorbers-and-the-grid-point-trap)) |
| `ZOptLowEnergy` `ZOptLowFactor` | real | none (exactly one when the absorber is on) | ZOptLowEnergy (eV, ≥ 0 within 1e-12; 0 means no absorption), vabszblow (strength in hartree) |
| `ZSpOptLowStart` `ZSpOptLowEnd` `ZSpOptLowEnergy` `ZSpOptLowFactor` | real | none: off | ssplmilow, ssplmalow, ZSpOptLowEnergy, vabssblow: the same on the specular grid, also with End ≤ Zana. Where it overlaps the high-end specular absorber (possible when Zana > ZSpOptStart), the low-end value replaces the high-end one, as in the Fortran |
| `DoDump` `FullSMatrix` | bool | false | |

## Checkpoint / restart keys (C++ only)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `RestartDir` | string | — | output directory of the run to continue; required if `NewCalc` is false. It must differ from `--output-dir`. For `onlyAnalysis`: the directory with the BK file to analyse (required; it may equal `--output-dir`) |
| `RestartFormat` | string | `"native"` | `"native"`: checkpoint.bin + bk_coefficients.bin; `"fortran"`: fort.8 (or fort.9) + fort.22 (or fort.21) written by the Fortran code or by this program. `onlyAnalysis` reads only the BK file |
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
`CRPlibrary` must use the same `CRPlibDataSet` (not checked).

## Not supported
These keys give a "not supported" error. Logical keys are accepted when set to `false`. The list
is: `SystemMix`, `SystemID2`, `SystemMixCoeff`, `UseRCutOff`,
`RCutOff`, `write2DPES`, `x2DPES`, `y2DPES`, `theta2DPES`, `phi2DPES`, `PESMixingAlfa`.

## Validation
All checks of `InputNew.f90` ReadInput/ValidateInput are reproduced, with the same tolerant real
comparisons (absolute tolerance 1e-12). So are the grid checks of `Setup.f90` and `setupdvr.f90`
and the nφ checks of `transform.f90` (which are hard errors here). Every check has a unit test in
`tests/test_input.cpp`. Deviations are listed in [fortran_deviations.md](fortran_deviations.md).
Two Fortran flux analysis checks are wrong and are fixed here (flux Q1, flux Q2), and so are
the silent cases of the low-Z absorber keys (low-Z D1–D3) and the default strength the Fortran
keeps for an absorber energy of 0 (absorber D1):
[fortran_input_bugs.md](fortran_input_bugs.md) describes them, with demonstration runs of the
Fortran program.
