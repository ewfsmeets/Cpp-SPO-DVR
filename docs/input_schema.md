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
| `SystemID` | string | systemID | PES; implemented: `HardWallPES` (case-insensitive) |
| `NewCalc` | bool | newcalc | new calculation; `false`: restart from `RestartDir` (tStart = end of the earlier run) |
| `ZDoF` `RDoF` `XDoF` `YDoF` `ThetaDoF` `PhiDoF` | bool | | must all be true (6D only) |
| `Mass1` `Mass2` | real | amass1/2 | atom masses (amu) |
| `HomoNuc` | bool | lhomo | homonuclear molecule |
| `ZStart` `dZ` | real | zmin, delz | Z grid start, spacing (shared by both Z grids) |
| `nZ` `nZSpecular` | int | nz, nzsp | scattering / specular Z points (valid FFT sizes, nZSpecular ≥ nZ) |
| `RStart` `dR` `nR` | real, real, int | rmin, delr, nr | r grid |
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
| `writeSMatrix` `writeProbs` `writeRDProbs` | bool | lwrit3, lwri25, lwdiro | output switches: `writeProbs` writes probabilities_rovib.tsv, `writeRDProbs` probabilities_diffraction.tsv; `writeSMatrix` has no effect yet (s_matrix.tsv is always written; in Fortran it controls the unit 3 output of the fine analysis) |
| `TolSp` | real | tolsp | specular-grid tolerance; also the maximum initial norm outside [ZAna, start of the specular absorber] (hard error) |

For each CAP, exactly one of `ZOptEnergy` / `ZOptFactor` (likewise `ROpt…`, `ZSpOpt…`) is required.
The energy is in eV, and the factor is the strength in hartree.

## Optional keys (Fortran defaults)

| Key | Type | Default | Fortran |
|---|---|---|---|
| `MJMax` | int | JMax | mjmaxw. A negative value means JMax. **0 also gives JMax** for the basis (Fortran quirk) |
| `VIn` `JIn` `MJIn` | int | 0 | initial state |
| `onlyAnalysis` `alsoAnalysis` | bool | false | run mode flags |
| `Potential` | bool | false | potential-only run (nread = 6) |
| `JAna` `VAna` `DifOrdAna` | int | determined by the program | analysis limits |
| `ZAnaSp` | real | none | fixed analysis line on the specular grid |
| `doFluxAnalysis` | bool | false | needs `rFluxAnalysis` `EMinFluxAna` `EMaxFluxAna` `DeltaEFluxAna` |
| `ImposeVMin` `VMin` | bool, real | false, −5.5 eV | clamp the potential from below |
| `ImposeVMax` `VMax` | bool, real | false, 10.0 eV | clamp the potential from above |
| `UseSwitching` `ZminSwitch` `ZmaxSwitch` | bool, real, real | false, 6.0, 8.0 | switch V to the gas phase; zasymp becomes ZmaxSwitch |
| `DoDump` `FullSMatrix` | bool | false | |

## Checkpoint / restart keys (C++ only)

| Key | Type | Default | Meaning |
|---|---|---|---|
| `RestartDir` | string | — | output directory of the run to continue; required if `NewCalc` is false. It must differ from `--output-dir` |
| `RestartFormat` | string | `"native"` | `"native"`: checkpoint.bin + bk_coefficients.bin; `"fortran"`: fort.8 (or fort.9) + fort.22 (or fort.21) written by the Fortran code or by this program |
| `CheckpointFormat` | string | `"native"` | files written for a later restart: `"native"`, `"fortran"` (fort.9, fort.21) or `"both"`. bk_coefficients.bin is always written |

A restart must use the same input as the earlier run except `NewCalc`, `tStart`, `tEnd`, `nStep`
(and output switches); the program checks the grids, masses, time steps, wave packet and
absorbers against the checkpoint and names any parameter that differs.

## Not supported
These keys give a "not supported" error. Logical keys are accepted when set to `false`. The list
is: `SystemMix`, `SystemID2`, `SystemMixCoeff`, `UseSymmetry`, `nPhisOptimal`, `UseRCutOff`,
`RCutOff`, `write2DPES`, `x2DPES`, `y2DPES`, `theta2DPES`, `phi2DPES`, `PESMixingAlfa`,
`CRPlibDataSet`.

## Validation
All checks of `InputNew.f90` ReadInput/ValidateInput are reproduced, with the same tolerant real
comparisons (absolute tolerance 1e-12). So are the grid checks of `Setup.f90` and `setupdvr.f90`
and the nφ checks of `transform.f90` (which are hard errors here). Every check has a unit test in
`tests/test_input.cpp`. Deviations are listed in CLAUDE.md.
