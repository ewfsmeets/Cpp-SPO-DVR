# Port status: the Fortran SPO-DVR features in the C++ program

The C++ program is a port of the Fortran SPO-DVR code. This page lists every feature of the Fortran
code (5.24 r1713) and its status in the C++, then what the C++ adds.

**Which Fortran version.** The port started from version 5.19 (2023). The gates compare with
upstream 5.24 at svn r1713 (2026), built with `-O3` ([reference data](validation.md#reference-data)).
The production-grid validation and the
[Fortran benchmark](performance.md#fortran-vs-c-full-dailight-run-8-threads-2026-11-05) used 5.19.
That makes no difference there: built with `-O3`, 5.24 is bitwise identical to 5.19 on every gate
case of 5.19 features, the benchmark case among them, and on the production grid's potential and
first 40 steps ([production grid](validation.md#production-grid)). r1713 adds a sign fix in the
fine-analysis S matrix; the screen labels differ.

**Where else to look:**

- each upstream change since 5.19 and its port status: [upstream_sync.md](upstream_sync.md);
- every known difference from the Fortran code: [fortran_deviations.md](fortran_deviations.md);
- the gates that compare the two programs: [validation.md](validation.md#gates).

**Status values:**

- **ported**: the C++ does what the Fortran does; a note in parentheses names a Fortran quirk that
  is kept, or a deliberate difference in a single check or edge case (with its ID, where it has one);
- **deviation**: ported, but different on purpose;
- **partly ported**: the cell says what is missing;
- **not ported**: the cell says what such an input does in the C++;
- **not planned**: decided not to port; the cell gives the reason.

"Input error": the program stops with a message naming the key (exit code 2). "Refused": the
program stops with a "not supported" message (exit code 3).

## Run modes

| Feature | What it does | Status | Details |
|---|---|---|---|
| New run (`NewCalc: true`; Fortran nread 0) | propagates the initial wave packet and analyses it on the fly | **ported** | [validation](validation.md#gates) |
| Restart (`NewCalc: false`; nread 1) | continues a run from its saved wave packet; the S matrix covers the whole run | **ported** | [restart](../README.md#continuing-a-run-restart) |
| Analysis-only run (`onlyAnalysis`; nread 2) | the fine analysis of the BK coefficients of an earlier run, without propagation | **ported** | [output_files](output_files.md#fine-analysis-onlyanalysis-alsoanalysis) |
| Fine analysis after the run (`alsoAnalysis`) | the same analysis right after the propagation | **deviation**: the fine results are written next to the on-the-fly ones (the Fortran replaces them) | [fortran_deviations: Fine analysis, kept](fortran_deviations.md#fine-analysis-kept) |
| Potential-only run (`Potential`; nread 6) | computes the potential and the gas-phase states, then stops | **partly ported**: V is computed and checked; `writePotential` saves it for plots (`potential/`), but the Fortran `.pes` cache (below) is not ported | [input_schema](input_schema.md#optional-keys-fortran-defaults), [.pes cache](fortran_deviations.md#not-ported-pes-cache) |
| 2D potential cut (`write2DPES`, `x2DPES`, `y2DPES`, `theta2DPES`, `phi2DPES`) | writes V(r, Z) at fixed x, y, θ, φ to a file, then stops | **ported** (`pes_cut_2d.dat`; the plot tool draws it) | [input_schema](input_schema.md#optional-keys-fortran-defaults), [output_files](output_files.md#data-files-space-separated) |

## Grids and wave packet

| Feature | What it does | Status | Details |
|---|---|---|---|
| Six dimensions (`ZDoF`, `RDoF`, `XDoF`, `YDoF`, `ThetaDoF`, `PhiDoF`) | the coordinates Z, r, x, y, θ, φ; both codes accept only all six | **ported** | [input_schema](input_schema.md#required-keys) |
| Molecule (`Mass1`, `Mass2`, `HomoNuc`) | the atom masses; homonuclear (one parity of j) or heteronuclear | **ported** | [input_schema](input_schema.md#required-keys) |
| Z grids (`ZStart`, `dZ`, `nZ`, `nZSpecular`) | FFT grids in the molecule–surface distance Z: the scattering grid and the longer specular grid | **ported**; a specular grid shorter than the scattering grid is an input error | [inputs/README](../inputs/README.md#1-terms) |
| r grid (`RStart`, `dR`, `nR`) | FFT grid in the bond length r | **ported** (`RStart ≤ 0` rejected, as in 5.24) | [upstream_sync](upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| Surface grid (`nX`, `nY`) | FFT grid over the surface unit cell, whose lattice and skew angle the PES defines | **ported** | [input_schema](input_schema.md#required-keys) |
| Rotational basis (`JMax`, `MJMax`) | spherical harmonics up to JMax, on a Gauss–Legendre grid in θ and an FFT grid in φ | **ported** (`MJMax ≤ 0` means JMax, as in the Fortran) | [fortran_deviations: MJMax <= 0](fortran_deviations.md#mjmax-zero) |
| Gas-phase states | the potential far from the surface, V_asy(r), and its rovibrational states (Fourier grid Hamiltonian) | **ported** | [validation](validation.md#gates) |
| Initial wave packet (`EMinWP`, `EMaxWP`, `Z0`, `VIn`, `JIn`, `MJIn`, `TolSp`) | a Gaussian in Z at Z0 that spans the energies EMinWP to EMaxWP, times the initial state (v, j, mJ) | **ported** (`TolSp` is required: the Fortran has no default) | [inputs/README](../inputs/README.md#3-the-wave-packet-is-set-by-its-energy-range) |
| Oblique incidence (`EParallel`, `PhiParallel`) | the parallel translational energy and its azimuth, by a momentum shift on the x, y grid | **ported** | [input_schema](input_schema.md#required-keys) |

## Propagation

| Feature | What it does | Status | Details |
|---|---|---|---|
| Split-operator propagation (`tStart`, `tEnd`, `dtStep`, `nStep`) | exp(−iT·dt/2) exp(−iV·dt) exp(−iT·dt/2) per step: the kinetic energy in Z, r, x and y by FFTs; the rotational energy, diagonal in the rotational basis, in half steps on either side of V, which acts on the angular grid | **ported** | [validation](validation.md#gates) |
| Two grids (`ProjOp: 2`) | the initial channel on the long specular grid, the rest of the wave on the scattering grid | **ported** | [input_schema](input_schema.md#required-keys) |
| Single grid (`ProjOp: 1`, `0`) | 1: moves the specular wave onto the scattering grid once little of it lies beyond that grid; 0: one grid from the start | **not planned**: 1 is a legacy option, and the Fortran rejects 0 itself; both are input errors | [fortran_deviations: ProjOp = 1](fortran_deviations.md#projop-1) |

## Potentials (`SystemID`)

Ported:

| Feature | What it does | Status | Details |
|---|---|---|---|
| `HardWallPES` | test potential: a 1/Z⁴ wall on an H₂ curve in r, without x, y, θ, φ dependence | **ported** | [validation](validation.md#gates) |
| `H2Cu111DaiLight` | H₂/Cu(111): the analytic LEPS PES of Dai and Zhang (J. Chem. Phys. 102, 6280 (1995); the DZ PES of the thesis), called Dai–Light in the Fortran and in the test and case names | **ported** | [validation](validation.md#gates) |
| `H2Cu111Persson` | H₂/Cu(111): analytic LEPS of Persson et al. (J. Chem. Phys. 110, 2240 (1999); the ER PES of the thesis, ER = Eley–Rideal) | **ported** (a lattice-period quirk kept) | [fortran_deviations: Persson PES quirk](fortran_deviations.md#persson-pes-quirk) |
| `CRPlibrary` (`CRPlibDataSet`) | any C3v, C4v or C6v data set of the Fortran CRP library (M. Wijzenbroek's implementation of the corrugation reducing procedure of Olsen et al., J. Chem. Phys. 116, 3841 (2002)) | **ported** ([deviations](fortran_deviations.md#crp-library-srccrplib-systemid-crplibrary)) | [crp_library](crp_library.md) |
| `CRPDistrtSurface` (`CRPlibDataSet` and the SCM keys; new in 5.20) | a CRP data set plus the static corrugation model (SCM) of one snapshot of a distorted slab | **ported** (position files v2/v3 only, SCM D5; [SCM deviations](fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface)) | [scm](scm.md), [keys](input_schema.md#scm-keys-systemid-crpdistrtsurface) |

Not ported: each of these SystemIDs is an input error. The four marked `PESMixingAlfa` take that
key; giving it is refused. The data files are in [`potentials/`](../potentials/), one directory
per system.

| SystemID | System and method | Data |
|---|---|---|
| `H2Cu111PW91` | H₂/Cu(111): CRP interpolation of PW91 DFT data | [H2Cu111PW91.dir](../potentials/H2Cu111PW91.dir) |
| `H2Cu111RPBE` | H₂/Cu(111): CRP interpolation of RPBE DFT data | [H2Cu111RPBE.dir](../potentials/H2Cu111RPBE.dir) |
| `H2Cu111PW91RPBE` (`PESMixingAlfa`) | H₂/Cu(111): α·PW91 + (1 − α)·RPBE of the two above (SRP) | [H2Cu111PW91RPBE.dir](../potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111SRPnMode4` (`PESMixingAlfa`) | 4-mode expansion of `H2Cu111PW91RPBE` about its barrier | [H2Cu111PW91RPBE.dir](../potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111SRPnMode5` (`PESMixingAlfa`) | 5-mode expansion of `H2Cu111PW91RPBE` about its barrier | [H2Cu111PW91RPBE.dir](../potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111Grow` | H₂/Cu(111): GROW fit (Collins) of DFT data | [H2Cu111Grow.dir](../potentials/H2Cu111Grow.dir) |
| `H2Cu110` | H₂/Cu(110): interpolation (A. Salin) of its data files | [H2Cu110.dir](../potentials/H2Cu110.dir) |
| `H2Cu100PES5` | H₂/Cu(100), PES number 5: interpolated ab initio data, corrugation-reducing and direct interpolation joined by a switching function | [PES5_INPUT_FILES.dir](../potentials/PES5_INPUT_FILES.dir) |
| `H2Pt111` | H₂/Pt(111): interpolated ab initio data, direct and corrugation-reducing interpolation joined by a switching function | [H2Pt111.dir](../potentials/H2Pt111.dir) |
| `H2Pt111nMode4` | 4-mode expansion of `H2Pt111` about its barrier | [H2Pt111.dir](../potentials/H2Pt111.dir) |
| `H2Pt111nMode5` | 5-mode expansion of `H2Pt111` about its barrier | [H2Pt111.dir](../potentials/H2Pt111.dir) |
| `HDPt111` | HD/Pt(111): `H2Pt111` in the HD centre-of-mass frame | [H2Pt111.dir](../potentials/H2Pt111.dir) |
| `H2Pd111` | H₂/Pd(111): interpolated ab initio data | [H2Pd111.dir](../potentials/H2Pd111.dir) |
| `H2Ru0001PW91` | H₂/Ru(0001): corrugation-reducing interpolation of PW91 data | [H2Ru0001.dir](../potentials/H2Ru0001.dir) |
| `H2Ru0001RPBE` | H₂/Ru(0001): corrugation-reducing interpolation of RPBE data | [H2Ru0001RPBE.dir](../potentials/H2Ru0001RPBE.dir) |
| `H2Ru0001PW91RPBE` (`PESMixingAlfa`) | H₂/Ru(0001): α·PW91 + (1 − α)·RPBE of the two above | [H2Ru0001PW91RPBE.dir](../potentials/H2Ru0001PW91RPBE.dir) |
| `H2CORu0001Grow` | H₂ on CO-covered Ru(0001): GROW fit | [H2CORu0001Grow.dir](../potentials/H2CORu0001Grow.dir) |
| `H2NiAl110` | H₂/NiAl(110): interpolated ab initio data, on a rectangular unit cell | [H2NiAl110.dir](../potentials/H2NiAl110.dir) |
| `H2TiAl100GrowPES` | H₂/TiAl(100): GROW fit | [H2TiAl100Grow.dir](../potentials/H2TiAl100Grow.dir) |
| `H2Al110` | H₂/Al(110), adapted from the Cu(110) module | [H2Al110.dir](../potentials/H2Al110.dir) |
| `H2W100` | H₂/W(100): corrugation-reducing interpolation (a "vdW" PES in its source) | [H2W100.dir](../potentials/H2W100.dir) |
| `H2W110` (new in 5.22) | H₂/W(110): corrugation-reducing interpolation (a "vdW" PES in its source) | [H2W110.dir](../potentials/H2W110.dir) |
| `HW110` (new in 5.23) | an H atom on W(110): a test potential built on `H2W110` | [H2W110.dir](../potentials/H2W110.dir) |

## Potential options

| Feature | What it does | Status | Details |
|---|---|---|---|
| Potential limits (`ImposeVMin`, `VMin`, `ImposeVMax`, `VMax`) | clamps V from below or above | **ported** | [input_schema](input_schema.md#optional-keys-fortran-defaults) |
| Switching to the gas phase (`UseSwitching`, `ZminSwitch`, `ZmaxSwitch`) | switches V over to V_asy(r) between ZminSwitch and ZmaxSwitch; the asymptotic region then starts at ZmaxSwitch | **ported** | [input_schema](input_schema.md#optional-keys-fortran-defaults) |
| Potential cache (`.pes` file) | saves V on disk and reads it back in a later run on the same grids | **not ported**: V is always computed; a `.pes` file is ignored | [fortran_deviations: .pes cache](fortran_deviations.md#not-ported-pes-cache), [upstream_sync](upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| PES morphing (`pesmorphparameters.dat`; new in 5.21) | scales or shifts expansion coefficients of V with Gaussians in Z and r | **not ported**: the file is ignored, without a warning | [upstream_sync](upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| PES mixing, r cut-off (`SystemMix`, `SystemID2`, `SystemMixCoeff`, `UseRCutOff`, `RCutOff`) | read by the Fortran, but without effect (no mixing code; the cut-off is forced off) | **not planned**: nothing to port; refused when switched on or given | [fortran_deviations: Not planned: PES mixing, r cut-off](fortran_deviations.md#not-planned-pes-mixing-and-r-cut-off) |

## Surface symmetry

| Feature | What it does | Status | Details |
|---|---|---|---|
| C6v symmetry (`UseSymmetry`) | computes and stores V only at the symmetry-unique (x, y) points | **ported** | [input_schema](input_schema.md#optional-keys-fortran-defaults) |
| C2v and C4v symmetry (`UseSymmetry`) | the same for these point groups; the Fortran switches it off with a warning | **deviation**: used | [fortran_deviations: Surface symmetry, deviations](fortran_deviations.md#surface-symmetry-deviations) |
| Optimal φ grid (`nPhisOptimal`) | raises nφ to a multiple that the whole point group can use | **ported** | [input_schema](input_schema.md#optional-keys-fortran-defaults) |

## Absorbers

| Feature | What it does | Status | Details |
|---|---|---|---|
| High-end absorbers (`ZOpt…`, `ROpt…`, `ZSpOpt…`: Start, End, Energy or Factor) | negative imaginary potentials, quadratic in the coordinate, at large Z on both Z grids and at large r | **ported** ([absorber D1](fortran_input_bugs.md#absorber-d1-an-absorber-energy-of-0-keeps-the-default-strength), [absorber D2](fortran_deviations.md#absorber-d2)) | [inputs/README](../inputs/README.md#7-absorbers) |
| Low-Z absorbers (`ZOptLow…`, `ZSpOptLow…`; new in 5.23) | the same, mirrored, at the low-Z end of both Z grids | **ported** ([low-Z D1–D4](fortran_input_bugs.md#the-low-z-absorbers-of-upstream-524)) | [inputs/README](../inputs/README.md#low-z-absorbers-and-the-grid-point-trap) |

## Analyses

| Feature | What it does | Status | Details |
|---|---|---|---|
| Analysis line (`Zana`, `dtAna`) | projects ψ at Z = Zana on the asymptotic channels every dtAna (Balint-Kurti coefficients) | **ported** | [validation](validation.md#gates) |
| Specular analysis line (`DynamZSpAna`, `ZAnaSp`) | where the initial channel is analysed on the specular grid: following the packet's tail inwards, or fixed | **ported** (Fortran quirks kept) | [inputs/README](../inputs/README.md#6-the-specular-analysis-line) |
| Analysed channels (`VAna`, `JAna`, `DifOrdAna`) | the vibrational, rotational and diffraction states analysed; by default the open ones | **ported** | [input_schema](input_schema.md#optional-keys-fortran-defaults) |
| Diffraction orders | the order of each diffraction channel, per surface point group | **deviation**: derived from the reciprocal lattice; the same orders | [fortran_deviations: Diffraction order](fortran_deviations.md#diffraction-order) |
| Analysis energies (`NrOfEs`) | the collision energies at which S is computed; −95 puts them on the grid of multiples of π/tEnd | **ported** (NrOfEs = 0 and 1 are input errors: [analysis D1](fortran_deviations.md#analysis-d1)) | [input_schema](input_schema.md#required-keys) |
| S matrix and probabilities | the energy-resolved S matrix from the BK coefficients during the run: state-to-state, diffraction and reaction probabilities | **ported** | [validation](validation.md#gates) |
| Mean final energies | the mean vibrational, rotational and diffraction energy of the scattered molecules | **partly ported**: the means per final v are not written | [fortran_deviations: Output not reproduced](fortran_deviations.md#output-not-reproduced) |
| Fine analysis | FFT of the BK time series (the records up to tEnd, cut to the largest valid FFT length that fits) onto the spectral grid, interpolated to the analysis energies | **deviation**: two Fortran interpolation bugs fixed | [fortran_fine_analysis_bugs](fortran_fine_analysis_bugs.md) |
| Optional tables (`writeProbs`, `writeRDProbs`, `writeSMatrix`) | the state-to-state probabilities per (v, j, mJ, n, m); the probabilities per diffraction order; the S matrix on the spectral grid | **ported** (as TSV files; the C++ always writes `s_matrix.tsv`, and `writeProbs` adds `probabilities_rovib.tsv`) | [output_files](output_files.md#results-tab-separated-tables) |
| Flux analysis (`doFluxAnalysis`, `rFluxAnalysis`, `EMinFluxAna`, `EMaxFluxAna`, `DeltaEFluxAna`) | the reaction probability from the flux through r = r_c, and its (x, y) density | **not ported**: a propagating run refuses it; its input checks are ported, two of them corrected ([flux Q1](fortran_deviations.md#flux-q1), [flux Q2](fortran_deviations.md#flux-q2)) | [flux_analysis_findings](flux_analysis_findings.md#1-status-2026-11-14) |
| Averages (`dtAve`) | the norm, ⟨Z⟩, ⟨r⟩ and the rest norm (between the two analysis lines) every dtAve | **ported** | [validation](validation.md#gates) |

## Restart files

| Feature | What it does | Status | Details |
|---|---|---|---|
| Fortran restart files (fort.9 → fort.8: wave packet; fort.21 → fort.22: BK coefficients) | what a run writes for a restart or a fine analysis, and what those read | **ported**: read and written (`RestartFormat`, `CheckpointFormat`) | [input_schema](input_schema.md#checkpoint--restart-keys-c-only) |
| Analysis state across a restart | the state of the specular analysis line carried into the restarted run | **deviation**: a Fortran bug fixed (`anlyTol` is not saved) | [fortran_deviations: Fortran restart bug, fixed](fortran_deviations.md#fortran-restart-bug) |

## Input and output

| Feature | What it does | Status | Details |
|---|---|---|---|
| Input file (`DVR.inp`, `InputNew.f90`; `Title`) | the run parameters as named fields with defaults | **deviation**: a JSON file with the Fortran's key names, case-sensitive | [input_schema](input_schema.md), [fortran_deviations: Input file](fortran_deviations.md#input-file) |
| Input checks (`InputNew.f90`, `Setup.f90`, `setupdvr.f90`, `transform.f90`) | stop on inconsistent parameters before the run | **ported**; stricter where the Fortran continues silently; exit code 2 (the Fortran's is 0) | [input_schema](input_schema.md#validation), [fortran_input_bugs](fortran_input_bugs.md) |
| Old input reader (`Input.f90`) | the earlier input format | **not planned**: dead code (the build uses `InputNew.f90`) | [fortran_deviations: Not ported by decision](fortran_deviations.md#not-ported-by-decision) |
| Debug output (`DoDump`; `FullSMatrix`) | `DoDump` writes the gas-phase states to `debug_asyp.dat`; `FullSMatrix` is only printed | **deviation**: both keys are ignored; the states are always written, to `gas_phase_states.dat` | [fortran_deviations: Output not reproduced](fortran_deviations.md#output-not-reproduced) |
| Screen output | input echo, grids, gas-phase states, averages, results per energy | **deviation**: a new layout with the same quantities, except those listed under "Output not reproduced" | [README: Running](../README.md#running), [fortran_deviations: Output not reproduced](fortran_deviations.md#output-not-reproduced) |
| Result files (fort.25, fil10, unit 3, `averages.out`) | S matrix, probabilities, spectrum and averages on disk | **deviation**: TSV and JSON files, and `averages.dat` (space-separated) for `averages.out` (additions below) | [output_files](output_files.md#correspondence-with-the-fortran-output) |
| CPU timings (`CpuTimings.f90`, fort.99) | the time spent per operator and per stage | **partly ported**: only the times of the potential, the propagation and the fine analysis | [fortran_deviations: Output not reproduced](fortran_deviations.md#output-not-reproduced) |
| Post-processing tools (upstream `tools.dir`) | read fil10 and fort.25, compare `.pes` files, site-specific reaction (flux density), n-H₂ averaging, plots, FFT sizes, wave-packet set-up, source-tree scripts, a PES-morphing example | **not ported**: the TSV files need no reader; the input checks FFT sizes; `tools/wave_packet_check.py` checks the wave packet; the wave-packet designer (`tools/wave_packet_designer/`) designs and checks the Z and r grids, the packet and the absorbers | [inputs/README](../inputs/README.md) |

## Parallelism and build

| Feature | What it does | Status | Details |
|---|---|---|---|
| OpenMP | threads share the grid points | **ported**: the results do not depend on the thread count | [performance](performance.md#threads-and-loop-scheduling) |
| MPI (`MPI_*.f90`, `USE_MPI`) | distributes the grid over several computers | **not planned**: OpenMP and the GPU instead | [fortran_deviations: Not ported by decision](fortran_deviations.md#not-ported-by-decision) |
| FFTs (`ffts.f90`, `fftw_2_scs.f90`) | the FFTs of the grids (batched 1D: complex, real↔complex, sine, cosine), through FFTW or oneMKL | **ported**: an FFT layer of its own (batched 1D and 2D; oneMKL's FFTW3 interface or FFTW) | [performance](performance.md#build-options-cmake) |
| Platform code (`*_hack.f90`; the ESSL and SCSL FFTs) | thread queries, timers and vendor FFT libraries per machine type | **not planned**: standard C++ and OpenMP instead | [fortran_deviations: Not ported by decision](fortran_deviations.md#not-ported-by-decision) |
| Build flags (upstream `GNU.mk` since 5.24) | `-Ofast -march=native`: fast-math and FMA | **not ported**: not relevant; strict floating point without FMA, and the reference is built with `-O3` | [upstream_sync](upstream_sync.md#port-status-of-the-upstream-changes-since-519) |

## Additions in the C++ program

| Addition | What it does | Details |
|---|---|---|
| GPU backend (`--device cuda`) | propagates on an NVIDIA GPU with the wave packet in GPU memory; agrees with the CPU code to rounding | [cuda](cuda.md) |
| C3v symmetry (C3 too, for a future PES); PES self-check (`SymmetryTolerance`) | symmetry-reduced V for point groups the Fortran has no code for; checks that the PES has the symmetry before using it | [input_schema](input_schema.md#optional-keys-fortran-defaults) |
| Native restart files (`checkpoint.bin`, `bk_coefficients.bin`; `RestartDir`) | restart files that record the run parameters; a restart checks them and names any that differ | [output_files](output_files.md#checkpoint-files-binary) |
| TSV and JSON results | tables with commented headers in full precision, and the run metadata in `analysis.json` | [output_files](output_files.md) |
| Loop scheduling (`--threads`, `--schedule`) | shares the uneven work per (Z, r) point over threads and NUMA nodes; the results do not depend on it | [performance](performance.md#threads-and-loop-scheduling) |
| Set-up reports | a memory estimate before V is computed; the ranges of a CRP data set; the SCM checks; the absorber grid-point warning | [performance](performance.md#memory-and-numa), [scm](scm.md#the-scm-set-up-report) |
| Build variants (`SPODVR_FAST_MATH`, PGO) | optional faster builds; the default build is strict floating point without FMA | [performance](performance.md#build-options-cmake) |
| Checking tools | compare the final wave packet of two builds or devices; check a wave packet; regenerate and compare the Fortran references | [CONTRIBUTING](CONTRIBUTING.md#testing), [fortran_reference](../tools/fortran_reference/README.md) |
| Wave-packet designer (`tools/wave_packet_designer/`) | designs and checks the Z and r part of an input before any run: grids, wave packet, absorbers, memory | [README](../tools/wave_packet_designer/README.md) |
| Plots of a run (`tools/plots/`), `run_info.json` and plot data | figures and a report of a finished run (restart segments joined), opened by a verdict on its completeness and convergence; every run writes its metadata for them; `writeDensities` records the probability per (Z, r) cell every dtAve (densities movie, carpets, norm per region); `writePotential` writes V on the run's grid at set-up (maps over the unit cell, elbows with the minimax barrier, V against the orientation at `PlotSites`) | [README](../tools/plots/README.md), [output_files](output_files.md#run-metadata), [series directories](output_files.md#series-directories-binary-plot-data), [plot data keys](input_schema.md#plot-data-keys-c-only) |
