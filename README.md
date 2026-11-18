```text
>>=============================================================================================<<
 _____                  _____  _   _ ______   ___     _____ ______  _____   ______  _   _ ______
/  __ \   _      _     /  __ \| | | ||  _  \ / _ \   /  ___|| ___ \|  _  |  |  _  \| | | || ___ \
| /  \/ _| |_  _| |_   | /  \/| | | || | | |/ /_\ \  \ `--. | |_/ /| | | |  | | | || | | || |_/ /
| |    |_   _||_   _|  | |    | | | || | | ||  _  |   `--. \|  __/ | | | |  | | | || | | ||    /
| \__/\  |_|    |_|    | \__/\| |_| || |/ / | | | |  /\__/ /| |    \ \_/ /  | |/ / \ \_/ /| |\ \
 \____/                 \____/ \___/ |___/  \_| |_/  \____/ \_|     \___/   |___/   \___/ \_| \_|
>>=============================================================================================<<
```

# Cpp-SPO-DVR

A modern C++20 rewrite of the Fortran SPO-DVR code: six-dimensional quantum wave-packet dynamics
of a diatomic molecule scattering from a surface (coordinates Z, r, x, y, θ, φ). It uses
split-operator propagation, FFT grids in Z, r, x and y, and a Gauss-Legendre/FFT DVR in θ and φ.

## About this project

This project is a full reimplementation of the original SPO-DVR Fortran 90 code, written by
Dr. Mark F. Somers during his PhD at Leiden University:

- M. F. Somers, *Six-Dimensional Quantum Dynamics of Dissociative Chemisorption of Hydrogen on
  Metal Surfaces*, PhD thesis, Leiden University (2004),
  [repository](https://scholarlypublications.universiteitleiden.nl/handle/1887/3250648).
- G. J. Kroes and M. F. Somers, *Six-dimensional dynamics of dissociative chemisorption of H₂ on
  metal surfaces*, J. Theor. Comput. Chem. 4, 493 (2005),
  [DOI](https://www.worldscientific.com/doi/abs/10.1142/S0219633605001647).

I am very grateful that I could use the original Fortran code during my own PhD. Making it faster
and porting it to CUDA has been a years-long endeavour: faster calculations open the door to more
avenues of approach with quantum dynamics calculations.

Over the years SPO-DVR has been my go-to code whenever I wanted to apply newly learned programming
skills to highly complex, performance-critical code. I tinkered with and refactored small parts of
it by hand, but it took Claude to pull all the ideas together into a working program and to
benchmark its accuracy. Without Claude this would have been too big for a hobby project.

## Status

This program is a port of the Fortran SPO-DVR code. The port started from version 5.19 (2023);
every regression gate compares with upstream 5.24 at svn r1713 (2026), built with `-O3`
([docs/validation.md](docs/validation.md#reference-data)). The production-grid validation and the
Fortran benchmark below used 5.19. Built with `-O3`, 5.24 is bitwise identical to it on every gate
case (the benchmark case among them) and on the production grid's potential and first 40 steps
([reference data](docs/validation.md#reference-data),
[production grid](docs/validation.md#production-grid)). The port status of each upstream change
since 5.19 is in [docs/upstream_sync.md](docs/upstream_sync.md), every known difference from the
Fortran code in [docs/fortran_deviations.md](docs/fortran_deviations.md), and the gates that
compare the two programs in [docs/validation.md](docs/validation.md#gates).

The table lists the features of the Fortran code (5.24 r1713) and their status in the C++:

- **ported**: the C++ does what the Fortran does; a note in parentheses names a Fortran quirk
  that is kept, or a deliberate difference in a single check or edge case (with its ID);
- **deviation**: ported, but different on purpose;
- **partly ported**: the cell says what is missing;
- **not ported**: the cell says what such an input does in the C++;
- **not planned**: decided not to port; the cell gives the reason.

"Input error": the program stops with a message naming the key (exit code 2). "Refused": the
program stops with a "not supported" message (exit code 3).

| Feature | What it does | Status | Details |
|---|---|---|---|
| **Run modes** | | | |
| New run (`NewCalc: true`; Fortran nread 0) | propagates the initial wave packet and analyses it on the fly | **ported** | [validation](docs/validation.md#gates) |
| Restart (`NewCalc: false`; nread 1) | continues a run from its saved wave packet; the S matrix covers the whole run | **ported** | [restart](#continuing-a-run-restart) |
| Analysis-only run (`onlyAnalysis`; nread 2) | the fine analysis of the BK coefficients of an earlier run, without propagation | **ported** | [output_files](docs/output_files.md#fine-analysis-onlyanalysis-alsoanalysis) |
| Fine analysis after the run (`alsoAnalysis`) | the same analysis right after the propagation | **deviation**: the fine results are written next to the on-the-fly ones (the Fortran replaces them) | [fortran_deviations: Fine analysis, kept](docs/fortran_deviations.md#fine-analysis-kept) |
| Potential-only run (`Potential`; nread 6) | computes the potential and the gas-phase states, then stops | **partly ported**: V is computed and checked but not saved (the `.pes` cache, below, is not ported) | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| 2D potential cut (`write2DPES`, `x2DPES`, `y2DPES`, `theta2DPES`, `phi2DPES`) | writes V(r, Z) at fixed x, y, θ, φ to a file, then stops | **not ported**: refused | [input_schema](docs/input_schema.md#not-supported) |
| **Grids and wave packet** | | | |
| Six dimensions (`ZDoF`, `RDoF`, `XDoF`, `YDoF`, `ThetaDoF`, `PhiDoF`) | the coordinates Z, r, x, y, θ, φ; both codes accept only all six | **ported** | [input_schema](docs/input_schema.md#required-keys) |
| Molecule (`Mass1`, `Mass2`, `HomoNuc`) | the atom masses; homonuclear (one parity of j) or heteronuclear | **ported** | [input_schema](docs/input_schema.md#required-keys) |
| Z grids (`ZStart`, `dZ`, `nZ`, `nZSpecular`) | FFT grids in the molecule–surface distance Z: the scattering grid and the longer specular grid | **ported**; a specular grid shorter than the scattering grid is an input error | [inputs/README](inputs/README.md#1-terms) |
| r grid (`RStart`, `dR`, `nR`) | FFT grid in the bond length r | **ported** (`RStart ≤ 0` rejected, as in 5.24) | [upstream_sync](docs/upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| Surface grid (`nX`, `nY`) | FFT grid over the surface unit cell, whose lattice and skew angle the PES defines | **ported** | [input_schema](docs/input_schema.md#required-keys) |
| Rotational basis (`JMax`, `MJMax`) | spherical harmonics up to JMax, on a Gauss–Legendre grid in θ and an FFT grid in φ | **ported** (`MJMax ≤ 0` means JMax, as in the Fortran) | [fortran_deviations: MJMax <= 0](docs/fortran_deviations.md#mjmax-zero) |
| Gas-phase states | the potential far from the surface, V_asy(r), and its rovibrational states (Fourier grid Hamiltonian) | **ported** | [validation](docs/validation.md#gates) |
| Initial wave packet (`EMinWP`, `EMaxWP`, `Z0`, `VIn`, `JIn`, `MJIn`, `TolSp`) | a Gaussian in Z at Z0 that spans the energies EMinWP to EMaxWP, times the initial state (v, j, mJ) | **ported** (`TolSp` is required: the Fortran has no default) | [inputs/README](inputs/README.md#3-the-wave-packet-is-set-by-its-energy-range) |
| Oblique incidence (`EParallel`, `PhiParallel`) | the parallel translational energy and its azimuth, by a momentum shift on the x, y grid | **ported** | [input_schema](docs/input_schema.md#required-keys) |
| **Propagation** | | | |
| Split-operator propagation (`tStart`, `tEnd`, `dtStep`, `nStep`) | exp(−iT·dt/2) exp(−iV·dt) exp(−iT·dt/2) per step: T by FFTs in Z, r, x, y and in the rotational basis, V on the angular grid | **ported** | [validation](docs/validation.md#gates) |
| Two grids (`ProjOp: 2`) | the initial channel on the long specular grid, the rest of the wave on the scattering grid | **ported** | [input_schema](docs/input_schema.md#required-keys) |
| Single grid (`ProjOp: 1`, `0`) | 1: moves the specular wave onto the scattering grid once little of it lies beyond that grid; 0: one grid from the start | **not planned**: 1 is a legacy option, and the Fortran rejects 0 itself; both are input errors | [fortran_deviations: ProjOp = 1](docs/fortran_deviations.md#projop-1) |
| **Potentials** (`SystemID`) | | | |
| `HardWallPES` | test potential: a 1/Z⁴ wall on an H₂ curve in r, without x, y, θ, φ dependence | **ported** | [validation](docs/validation.md#gates) |
| `H2Cu111DaiLight` | H₂/Cu(111): the analytic LEPS PES of Dai and Zhang (J. Chem. Phys. 102, 6280 (1995); the DZ PES of the thesis), called Dai–Light in the Fortran and in the test and case names | **ported** | [validation](docs/validation.md#gates) |
| `H2Cu111Persson` | H₂/Cu(111): analytic LEPS of Persson et al. (J. Chem. Phys. 110, 2240 (1999); the ER PES of the thesis, ER = Eley–Rideal) | **ported** (a lattice-period quirk kept) | [fortran_deviations: Persson PES quirk](docs/fortran_deviations.md#persson-pes-quirk) |
| `CRPlibrary` (`CRPlibDataSet`) | any C3v, C4v or C6v data set of the CRP library of M. Wijzenbroek | **ported** ([deviations](docs/fortran_deviations.md#crp-library-srccrplib-systemid-crplibrary)) | [crp_library](docs/crp_library.md) |
| `CRPDistrtSurface` (`CRPlibDataSet` and the SCM keys; new in 5.24) | a CRP data set plus the static corrugation model (SCM) of one snapshot of a distorted slab | **ported** (position files v2/v3 only, SCM D5; [SCM deviations](docs/fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface)) | [scm](docs/scm.md), [keys](docs/input_schema.md#scm-keys-systemid-crpdistrtsurface) |
| `H2Cu111PW91` | H₂/Cu(111): CRP fit of PW91 DFT data | **not ported**: input error | [data](potentials/H2Cu111PW91.dir) |
| `H2Cu111RPBE` | H₂/Cu(111): CRP fit of RPBE DFT data | **not ported**: input error | [data](potentials/H2Cu111RPBE.dir) |
| `H2Cu111PW91RPBE` (`PESMixingAlfa`) | H₂/Cu(111): α·PW91 + (1 − α)·RPBE of the two above (SRP) | **not ported**: input error (refused if `PESMixingAlfa` is given) | [data](potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111SRPnMode4` (`PESMixingAlfa`) | 4-mode expansion of `H2Cu111PW91RPBE` about its barrier | **not ported**: input error (refused if `PESMixingAlfa` is given) | [data](potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111SRPnMode5` (`PESMixingAlfa`) | 5-mode expansion of `H2Cu111PW91RPBE` about its barrier | **not ported**: input error (refused if `PESMixingAlfa` is given) | [data](potentials/H2Cu111PW91RPBE.dir) |
| `H2Cu111Grow` | H₂/Cu(111): GROW fit (Collins) of DFT data | **not ported**: input error | [data](potentials/H2Cu111Grow.dir) |
| `H2Cu110` | H₂/Cu(110), from its data files | **not ported**: input error | [data](potentials/H2Cu110.dir) |
| `H2Cu100PES5` | H₂/Cu(100): corrugation-reducing fit, PES number 5 | **not ported**: input error | [data](potentials/PES5_INPUT_FILES.dir) |
| `H2Pt111` | H₂/Pt(111): interpolated ab initio data (direct or corrugation-reducing) | **not ported**: input error | [data](potentials/H2Pt111.dir) |
| `H2Pt111nMode4` | 4-mode expansion of `H2Pt111` about its barrier | **not ported**: input error | [data](potentials/H2Pt111.dir) |
| `H2Pt111nMode5` | 5-mode expansion of `H2Pt111` about its barrier | **not ported**: input error | [data](potentials/H2Pt111.dir) |
| `HDPt111` | HD/Pt(111): `H2Pt111` in the HD centre-of-mass frame | **not ported**: input error | [data](potentials/H2Pt111.dir) |
| `H2Pd111` | H₂/Pd(111): interpolated ab initio data | **not ported**: input error | [data](potentials/H2Pd111.dir) |
| `H2Ru0001PW91` | H₂/Ru(0001): corrugation-reducing fit of PW91 data | **not ported**: input error | [data](potentials/H2Ru0001.dir) |
| `H2Ru0001RPBE` | H₂/Ru(0001): corrugation-reducing fit of RPBE data | **not ported**: input error | [data](potentials/H2Ru0001RPBE.dir) |
| `H2Ru0001PW91RPBE` (`PESMixingAlfa`) | H₂/Ru(0001): α·PW91 + (1 − α)·RPBE of the two above | **not ported**: input error (refused if `PESMixingAlfa` is given) | [data](potentials/H2Ru0001PW91RPBE.dir) |
| `H2CORu0001Grow` | H₂ on CO-covered Ru(0001): GROW fit | **not ported**: input error | [data](potentials/H2CORu0001Grow.dir) |
| `H2NiAl110` | H₂/NiAl(110): interpolated ab initio data, on a rectangular unit cell | **not ported**: input error | [data](potentials/H2NiAl110.dir) |
| `H2TiAl100GrowPES` | H₂/TiAl(100): GROW fit | **not ported**: input error | [data](potentials/H2TiAl100Grow.dir) |
| `H2Al110` | H₂/Al(110), adapted from the Cu(110) module | **not ported**: input error | [data](potentials/H2Al110.dir) |
| `H2W100` | H₂/W(100): corrugation-reducing fit with van der Waals | **not ported**: input error | [data](potentials/H2W100.dir) |
| `H2W110` (new in 5.24) | H₂/W(110): corrugation-reducing fit with van der Waals | **not ported**: input error | [data](potentials/H2W110.dir) |
| `HW110` (new in 5.24) | an H atom on W(110): a test potential built on `H2W110` | **not ported**: input error | [data](potentials/H2W110.dir) |
| **Potential options** | | | |
| Potential limits (`ImposeVMin`, `VMin`, `ImposeVMax`, `VMax`) | clamps V from below or above | **ported** | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| Switching to the gas phase (`UseSwitching`, `ZminSwitch`, `ZmaxSwitch`) | switches V over to V_asy(r) between ZminSwitch and ZmaxSwitch; the asymptotic region then starts at ZmaxSwitch | **ported** | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| Potential cache (`.pes` file) | saves V on disk and reads it back in a later run on the same grids | **not ported**: V is always computed; a `.pes` file is ignored | [upstream_sync](docs/upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| PES morphing (`pesmorphparameters.dat`; new in 5.24) | scales or shifts expansion coefficients of V with Gaussians in Z and r | **not ported**: the file is ignored, without a warning | [upstream_sync](docs/upstream_sync.md#port-status-of-the-upstream-changes-since-519) |
| PES mixing, r cut-off (`SystemMix`, `SystemID2`, `SystemMixCoeff`, `UseRCutOff`, `RCutOff`) | read by the Fortran, but without effect (no mixing code; the cut-off is forced off) | **not planned**: nothing to port; refused when switched on or given | [fortran_deviations: Not planned: PES mixing, r cut-off](docs/fortran_deviations.md#not-planned-pes-mixing-and-r-cut-off) |
| **Surface symmetry** | | | |
| C6v symmetry (`UseSymmetry`) | computes and stores V only at the symmetry-unique (x, y) points | **ported** (nX or nY = 1 is an input error; the Fortran warns) | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| C2v and C4v symmetry (`UseSymmetry`) | the same for these point groups; the Fortran switches it off with a warning | **deviation**: used | [fortran_deviations: Surface symmetry, deviations](docs/fortran_deviations.md#surface-symmetry-deviations) |
| Optimal φ grid (`nPhisOptimal`) | raises nφ to a multiple that the whole point group can use | **ported** | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| **Absorbers** | | | |
| High-end absorbers (`ZOpt…`, `ROpt…`, `ZSpOpt…`: Start, End, Energy or Factor) | negative imaginary potentials, quadratic in the coordinate, at large Z on both Z grids and at large r | **ported** ([absorber D1](docs/fortran_input_bugs.md#absorber-d1-an-absorber-energy-of-0-keeps-the-default-strength), [absorber D2](docs/fortran_deviations.md#absorber-d2)) | [inputs/README](inputs/README.md#7-absorbers) |
| Low-Z absorbers (`ZOptLow…`, `ZSpOptLow…`; new in 5.24) | the same, mirrored, at the low-Z end of both Z grids | **ported** ([low-Z D1–D4](docs/fortran_input_bugs.md#the-low-z-absorbers-of-upstream-524)) | [inputs/README](inputs/README.md#low-z-absorbers-and-the-grid-point-trap) |
| **Analyses** | | | |
| Analysis line (`Zana`, `dtAna`) | projects ψ at Z = Zana on the asymptotic channels every dtAna (Balint-Kurti coefficients) | **ported** | [validation](docs/validation.md#gates) |
| Specular analysis line (`DynamZSpAna`, `ZAnaSp`) | where the initial channel is analysed on the specular grid: following the packet's tail inwards, or fixed | **ported** (Fortran quirks kept) | [inputs/README](inputs/README.md#6-the-specular-analysis-line) |
| Analysed channels (`VAna`, `JAna`, `DifOrdAna`) | the vibrational, rotational and diffraction states analysed; by default the open ones | **ported** | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| Diffraction orders | the order of each diffraction channel, per surface point group | **deviation**: derived from the reciprocal lattice; the same orders | [fortran_deviations: Diffraction order](docs/fortran_deviations.md#diffraction-order) |
| Analysis energies (`NrOfEs`) | the collision energies at which S is computed; −95 puts them on the grid of multiples of π/tEnd | **ported** (NrOfEs = 0 and 1 are input errors: [analysis D1](docs/fortran_deviations.md#analysis-d1)) | [input_schema](docs/input_schema.md#required-keys) |
| S matrix and probabilities | the energy-resolved S matrix from the BK coefficients during the run: state-to-state, diffraction and reaction probabilities | **ported** | [validation](docs/validation.md#gates) |
| Mean final energies | the mean vibrational, rotational and diffraction energy of the scattered molecules | **partly ported**: the means per final v are not written | [fortran_deviations: Output not reproduced](docs/fortran_deviations.md#output-not-reproduced) |
| Fine analysis | FFT of the whole BK time series onto the spectral grid, interpolated to the analysis energies | **deviation**: two Fortran interpolation bugs fixed | [fortran_fine_analysis_bugs](docs/fortran_fine_analysis_bugs.md) |
| Optional tables (`writeProbs`, `writeRDProbs`, `writeSMatrix`) | probabilities per (v, j, mJ) and per diffraction order; the S matrix on the spectral grid | **ported** (as TSV files) | [output_files](docs/output_files.md#results-tab-separated-tables) |
| Flux analysis (`doFluxAnalysis`, `rFluxAnalysis`, `EMinFluxAna`, `EMaxFluxAna`, `DeltaEFluxAna`) | the reaction probability from the flux through r = r_c, and its (x, y) density | **not ported**: a propagating run refuses it; its input checks are ported, two of them corrected ([flux Q1](docs/fortran_deviations.md#flux-q1), [flux Q2](docs/fortran_deviations.md#flux-q2)) | [flux_analysis_findings](docs/flux_analysis_findings.md#1-status-2026-11-14) |
| Averages (`dtAve`) | the norm, ⟨Z⟩, ⟨r⟩ and the rest norm (between the two analysis lines) every dtAve | **ported** | [validation](docs/validation.md#gates) |
| **Restart files** | | | |
| Fortran restart files (fort.9 → fort.8: wave packet; fort.21 → fort.22: BK coefficients) | what a run writes for a restart or a fine analysis, and what those read | **ported**: read and written (`RestartFormat`, `CheckpointFormat`) | [input_schema](docs/input_schema.md#checkpoint--restart-keys-c-only) |
| Analysis state across a restart | the state of the specular analysis line carried into the restarted run | **deviation**: a Fortran bug fixed (`anlyTol` is not saved) | [fortran_deviations: Fortran restart bug, fixed](docs/fortran_deviations.md#fortran-restart-bug) |
| **Input and output** | | | |
| Input file (`DVR.inp`, `InputNew.f90`; `Title`) | the run parameters as named fields with defaults | **deviation**: a JSON file with the Fortran's key names, case-sensitive | [input_schema](docs/input_schema.md), [fortran_deviations: Input file](docs/fortran_deviations.md#input-file) |
| Input checks (`InputNew.f90`, `Setup.f90`, `setupdvr.f90`, `transform.f90`) | stop on inconsistent parameters before the run | **ported**; stricter where the Fortran continues silently; exit code 2 (the Fortran's is 0) | [input_schema](docs/input_schema.md#validation), [fortran_input_bugs](docs/fortran_input_bugs.md) |
| Old input reader (`Input.f90`) | the earlier input format | **not planned**: dead code (the build uses `InputNew.f90`) | [fortran_deviations: Not ported by decision](docs/fortran_deviations.md#not-ported-by-decision) |
| Debug output (`DoDump`; `FullSMatrix`) | `DoDump` writes the gas-phase states to `debug_asyp.dat`; `FullSMatrix` is only printed | **deviation**: both keys are ignored; the states are always written, to `gas_phase_states.dat` | [fortran_deviations: Output not reproduced](docs/fortran_deviations.md#output-not-reproduced) |
| Screen output | input echo, grids, gas-phase states, averages, results per energy | **deviation**: a new layout with the same quantities | [Running](#running) |
| Result files (fort.25, fil10, unit 3, `averages.out`) | S matrix, probabilities, spectrum and averages on disk | **deviation**: TSV and JSON files (additions below) | [output_files](docs/output_files.md#correspondence-with-the-fortran-output) |
| CPU timings (`CpuTimings.f90`, fort.99) | the time spent per operator and per stage | **partly ported**: only the times of the potential, the propagation and the fine analysis | [fortran_deviations: Output not reproduced](docs/fortran_deviations.md#output-not-reproduced) |
| Post-processing tools (upstream `tools.dir`) | read fil10 and fort.25, compare `.pes` files, site-specific reaction (flux density), n-H₂ averaging, plots, FFT sizes, wave-packet set-up | **not ported**: the TSV files need no reader; the input checks FFT sizes; `tools/wave_packet_check.py` checks the wave packet | [inputs/README](inputs/README.md) |
| **Parallelism and build** | | | |
| OpenMP | threads share the grid points | **ported**: the results do not depend on the thread count | [performance](docs/performance.md#threads-and-loop-scheduling) |
| MPI (`MPI_*.f90`, `USE_MPI`) | distributes the grid over several computers | **not planned**: OpenMP and the GPU instead | [fortran_deviations: Not ported by decision](docs/fortran_deviations.md#not-ported-by-decision) |
| FFTs (`ffts.f90`, `fftw_2_scs.f90`) | the 1D and 2D FFTs of the grids, through FFTW or oneMKL | **ported**: an FFT layer of its own (oneMKL's FFTW3 interface or FFTW) | [performance](docs/performance.md#build-options-cmake) |
| Platform code (`*_hack.f90`; the ESSL and SCSL FFTs) | thread queries, timers and vendor FFT libraries per machine type | **not planned**: standard C++ and OpenMP instead | [fortran_deviations: Not ported by decision](docs/fortran_deviations.md#not-ported-by-decision) |
| Build flags (upstream `GNU.mk` since 5.24) | `-Ofast -march=native`: fast-math and FMA | **not ported**: not relevant; strict floating point without FMA, and the reference is built with `-O3` | [upstream_sync](docs/upstream_sync.md#port-status-of-the-upstream-changes-since-519) |

### Additions in the C++ program

| Addition | What it does | Details |
|---|---|---|
| GPU backend (`--device cuda`) | propagates on an NVIDIA GPU with the wave packet in GPU memory; agrees with the CPU code to rounding | [cuda](docs/cuda.md) |
| C3v and C3 symmetry; PES self-check (`SymmetryTolerance`) | symmetry-reduced V for point groups the Fortran has no code for; checks that the PES has the symmetry before using it | [input_schema](docs/input_schema.md#optional-keys-fortran-defaults) |
| Native restart files (`checkpoint.bin`, `bk_coefficients.bin`; `RestartDir`) | restart files that record the run parameters; a restart checks them and names any that differ | [output_files](docs/output_files.md#checkpoint-files-binary) |
| TSV and JSON results | tables with commented headers in full precision, and the run metadata in `analysis.json` | [output_files](docs/output_files.md) |
| Loop scheduling (`--threads`, `--schedule`) | shares the uneven work per (Z, r) point over threads and NUMA nodes; the results do not depend on it | [performance](docs/performance.md#threads-and-loop-scheduling) |
| Set-up reports | a memory estimate before V is computed; the ranges of a CRP data set; the SCM checks; the absorber grid-point warning | [performance](docs/performance.md#memory-and-numa), [scm](docs/scm.md#the-scm-set-up-report) |
| Build variants (`SPODVR_FAST_MATH`, PGO) | optional faster builds; the default build is strict floating point without FMA | [performance](docs/performance.md#build-options-cmake) |
| Checking tools | compare the final wave packet of two builds or devices; check a wave packet; regenerate and compare the Fortran references | [Testing](#testing), [fortran_reference](tools/fortran_reference/README.md) |

## Dependencies

- **Build tools**: a C++20 compiler with OpenMP (tested with GCC 15, the Ubuntu 26.04 default), CMake ≥ 3.24,
  Ninja, pkg-config (CMake looks up FFTW3 through it).
- **Math libraries**, one of:
  - **Intel oneMKL** (recommended; tested with 2026.2, found through its `MKLConfig.cmake`):
    BLAS/LAPACK and the FFTs (through its FFTW3 interface);
  - FFTW3 + LAPACK (fallback when oneMKL is not installed).
- nlohmann/json and GoogleTest are downloaded by CMake (internet access needed on the first
  configure).
- **Optional, for the GPU backend**: an NVIDIA GPU with its driver, and the CUDA toolkit (nvcc,
  cudart, cuFFT, and at run time nvJitLink and NVRTC for the cuFFT callback). Recommended: the newest
  CUDA 13 the GPU driver supports (tested with 13.3,
  the first release that accepts GCC 15; nvcc uses the same compiler as the rest of the build). FP64 speed
  matters: data-centre GPUs (A100, H100) run FP64 at 1/2 of the FP32 rate, GeForce cards at 1/64.

On Debian/Ubuntu (also under WSL2):

**1. Build tools** (always needed):

```bash
sudo apt install build-essential cmake ninja-build pkg-config
```

Any C++20 GCC works; the project is tested with GCC 15 (Ubuntu 26.04's default). For another
compiler, configure a new build directory with e.g. `CXX=g++-14 cmake ...`; nvcc then uses it too.

**2a. Intel oneMKL** (recommended), from Intel's apt repository:

```bash
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg]" \
     "https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-mkl-devel
```

CMake finds oneMKL in `$MKLROOT` or in `/opt/intel/oneapi/mkl/latest` (no `setvars.sh` needed for
building). The program always uses the sequential, LP64 oneMKL.

**2b. Without oneMKL**, FFTW3 and LAPACK instead:

```bash
sudo apt install libfftw3-dev liblapack-dev
```

**3. CUDA** (optional, for `--device cuda`), CUDA 13 from NVIDIA's apt repository (repository
`wsl-ubuntu` under WSL2, `ubuntu2604` / `ubuntu2404` on a native Ubuntu 26.04 / 24.04):

```bash
nvidia-smi                                  # the GPU must be visible; note "CUDA Version" (newest CUDA the driver runs)
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install cuda-toolkit-13-3          # newest 13.x not above the driver's CUDA version; installs into /usr/local/cuda
```

- Do not install a toolkit newer than the driver supports: it compiles, but the kernels then fail
  at run time.
- **WSL2:** the GPU driver is the Windows NVIDIA driver (update it in Windows). Never install a
  Linux driver inside WSL; the `wsl-ubuntu` repository contains none.
- `/usr/local/cuda/bin` need not be on `PATH`: CMake finds nvcc there. It also sets the GPU
  architecture from `nvidia-smi` (override with `-DCMAKE_CUDA_ARCHITECTURES`, e.g.
  `"80-real;90-real"` for A100 and H100).
- Ubuntu's own packages (`sudo apt install nvidia-cuda-toolkit`) work as well, but are older
  (Ubuntu 24.04: CUDA 12.0 with g++-12 as nvcc's host compiler; untested since the kernels became C++20).
- nvcc compiles with the C++ compiler of the build. If the CUDA release rejects that GCC
  ("unsupported GNU version"; GCC 15 needs CUDA ≥ 13.3), set `CUDAHOSTCXX` or
  `-DCMAKE_CUDA_HOST_COMPILER` to a GCC it accepts.

## Building

```bash
git clone https://github.com/ewfsmeets/Cpp-SPO-DVR.git
cd Cpp-SPO-DVR
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

The configure step reports which libraries it uses, e.g.
`SPO-DVR: BLAS/LAPACK from oneMKL 2026.2 (sequential), FFT library MKL`.
How to run the tests: see [Testing](#testing).

With the GPU backend:

```bash
cmake -S . -B build-cuda -G Ninja -DCMAKE_BUILD_TYPE=Release -DSPODVR_USE_CUDA=ON
cmake --build build-cuda
```

The kernels are compiled for the GPU of the build machine (`native`); for another GPU set
`-DCMAKE_CUDA_ARCHITECTURES` (`89-real` RTX 40xx, `80-real` A100, `90-real` H100). The configure step prints the nvcc version, host compiler and architectures.

Useful CMake options:

| Option | Default | Meaning |
|---|---|---|
| `-DSPODVR_USE_MKL=OFF` | `ON` | ignore oneMKL; use the system LAPACK and FFTW |
| `-DSPODVR_FFT_LIBRARY=FFTW` | `AUTO` | use the real FFTW for the FFTs even when oneMKL is present (`AUTO` picks oneMKL when found) |
| `-DSPODVR_ARCH=native` | `x86-64-v3` | `-march` value (`x86-64-v3` = AVX2; `x86-64-v4` = AVX-512) |
| `-DSPODVR_PGO=GENERATE/USE` | `OFF` | profile-guided optimisation, two builds and a training run (see docs/performance.md) |
| `-DSPODVR_FAST_MATH=ON` | `OFF` | experimental: `-ffast-math` with NaN and inf kept (with CUDA also `--use_fast_math`, i.e. FMA, in the GPU kernels); changes results at the 1e-12 to 1e-11 level ([details and known test failures](docs/performance.md#the-fast-math-build-spodvr_fast_math)) |
| `-DSPODVR_USE_CUDA=ON` | `OFF` | build the GPU backend (`--device cuda`, NVIDIA GPU, CUDA toolkit; see [docs/cuda.md](docs/cuda.md)) |
| `-DSPODVR_BUILD_TESTS=OFF` | `ON` | skip the tests |

Performance notes, reproducibility across machines and the PGO workflow are in
[docs/performance.md](docs/performance.md).

### Build variants and `tools/build_all.sh`

The development builds live side by side in git-ignored directories.
`tools/build_all.sh [variant ...]` configures and builds them (all by default; the CUDA ones only
when nvcc is found). It does not use PGO; for that, see docs/performance.md:

| Directory | Options | Purpose |
|---|---|---|
| `build` | (defaults) | the reference build: oneMKL, no FMA, LTO |
| `build-fftw` | `SPODVR_FFT_LIBRARY=FFTW` | oneMKL LAPACK with the real FFTW |
| `build-nomkl` | `SPODVR_USE_MKL=OFF` | system LAPACK + FFTW: the fallback that must keep working |
| `build-poison` | `SPODVR_POISON_UNINITIALISED=ON` | NaN-fills new large arrays: finds reads before the first write |
| `build-fast` | `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fast CPU build (results differ at the 1e-12 to 1e-11 level; add PGO by hand for the fastest) |
| `build-cuda` | `SPODVR_USE_CUDA=ON` | GPU backend, no FMA (results comparable with the CPU build) |
| `build-cuda-fast` | `SPODVR_USE_CUDA=ON`, `SPODVR_FAST_MATH=ON`, `SPODVR_ARCH=native` | fastest GPU build (FMA in the kernels; no PGO) |

Native restart files (`checkpoint.bin`, `bk_coefficients.bin`) only work between builds with the
same `SPODVR_FAST_MATH` setting. A restart or `onlyAnalysis` run checks the stored set-up parameters
exactly, and the fast-math builds compute one of them (the effective `ZOptFactor`) 1 ulp
differently, so a default build refuses a fast-math build's native files and vice versa. Files in
the Fortran format (`CheckpointFormat` `"fortran"` or `"both"`, read with `RestartFormat`
`"fortran"`) carry no parameter list and can be used across builds.

### Formatting and linting

clang-format (C++) and ruff (Python: lint and format, settings in `ruff.toml`) run as git
pre-commit hooks; clang-tidy is run on demand. All versions are pinned in `.pre-commit-config.yaml`
and installed by pre-commit itself (no system packages needed):

```bash
pipx install pre-commit                                      # once per machine
pre-commit install                                           # once per clone
pre-commit run --all-files                                   # format everything, lint the Python files
pre-commit run clang-tidy --hook-stage manual --all-files    # lint (needs build/compile_commands.json)
```

## Testing

```bash
OMP_NUM_THREADS=8 ctest --test-dir build --output-on-failure                # everything (~9.5 min)
OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R Cuda --output-on-failure   # GPU tests (~1.5 min)
```

The times: 809 tests (803 pass, 6 long or reference-less tests skipped) and 24 GPU tests, measured
2026-11-15 at 8b3df8d with 8 threads on an AMD Ryzen 7 5800X and an RTX 4070 Ti.

- **Run ctest serially** (no `-j`): many tests already use all OpenMP threads. Set
  `OMP_NUM_THREADS` to the number of physical cores (not hyperthreads).
- **What is tested:**
  - unit tests, one program per module (`tests/test_<module>.cpp`);
  - regression tests against the original Fortran code (`tests/test_regression_*.cpp`, test names
    `Regression…`/`Gate…`): potential, gas-phase states, propagation (norms and averages),
    on-the-fly analysis (probabilities, BK coefficients) and fine analysis. The Fortran reference
    output is stored in `tests/reference/` (the CRP and SCM driver outputs in
    `tests/crplib/reference/` and `tests/scm/reference/`). The gates, their tolerances and cases:
    [docs/validation.md](docs/validation.md);
  - GPU tests (`tests/cuda/`, names `Cuda…`): only in builds with `-DSPODVR_USE_CUDA=ON`. They skip on
    a machine without a GPU.
- **Selecting tests** (ctest `-R`/`-E` take a regular expression on the test name;
  `ctest --test-dir build -N` lists them all):
  ```bash
  ctest --test-dir build -E Regression           # unit tests only (~2.5 min)
  ctest --test-dir build -R Gate3                # one acceptance gate
  ./build/tests/test_input --gtest_filter='Input.*'   # one test program directly (GoogleTest)
  ```
- **Long tests:** the full 7000-step regression runs (hard wall and Dai–Light; the latter alone takes
  about 10 minutes on 8 cores), the negative tests of the low-Z and the SCM gates (1000- and
  3000-step runs) and the CRP cluster-pruning check on ten million points
  (`CrpLib.ClusterPruningLong.TenMillionPoints`) are skipped unless `SPODVR_LONG_TESTS` is set:
  ```bash
  SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build -R 'Long($|\.)' --output-on-failure
  SPODVR_LONG_TESTS=1 OMP_NUM_THREADS=8 ctest --test-dir build-cuda -R 'Cuda.*Long($|\.)' --output-on-failure
  ```
  The regular expression selects the names that end in `Long` (or whose suite does); a plain
  `-R Long` also runs the short `ScmFiles.ScmD6LongPathsAreNotTruncated`.
- **Skipped regression tests:** a regression test skips (it does not fail) when its Fortran
  reference file is missing; the tests of the SCM cases `scm_cu111_srp48_small*` fail instead. The
  scripts in `tools/fortran_reference/` regenerate the files; they need a Fortran reference install
  ([tools/fortran_reference/README.md](tools/fortran_reference/README.md); how to get the Fortran
  code: [there](tools/fortran_reference/README.md#getting-the-fortran-code)). Running the tests
  needs only the committed reference files.
  `ctest -R ReferenceManifest` fails when a file listed in a case's `MANIFEST.tsv` is missing or
  edited, or a case has no `MANIFEST.tsv` ([MANIFEST.tsv](tools/fortran_reference/README.md#manifesttsv)).
- **Other build variants** ([Build variants](#build-variants-and-toolsbuild_allsh)): run the same
  ctest command in their build directory. Use `build-poison` after touching the uninitialised-memory
  code. The fast-math builds (`build-fast`, `build-cuda-fast`) are expected to fail a few tests that
  compare exact values; their gates pass. The list, with the reason for each:
  [docs/performance.md](docs/performance.md#the-fast-math-build-spodvr_fast_math).
- **Docs consistency:** `python3 -I -B tools/check_docs.py [--spodvr build/spodvr]` (under 1 s)
  checks the docs against the code: input keys, output files, links and anchors, decision IDs,
  local paths, command-line options (with `--spodvr`), gates, environment variables, and the
  upstream commit of `docs/upstream_sync.md` against the pin in
  `tools/fortran_reference/reference_config.sh` (and the install the runs use against it); it also
  lists numbers repeated across docs (report-only). ctest runs it as `Docs.Consistency`
  (`ctest --test-dir build -R Docs`). Accepted findings, each with a reason, are in
  `tools/check_docs_allowlist.tsv`; a row that no longer suppresses anything fails the check.
- **Numerically neutral changes** (refactors) must leave the results bitwise unchanged:
  `tools/check_bitwise_refactor.sh [base-ref]` (default `HEAD`). It builds both trees in a scratch
  directory. For changes that may alter rounding: `tools/compare_propagation.sh [base-ref]
  [tolerance]` (its comparison, `tools/compare/compare_averages.py`, is tested by
  `ctest -R ComparePropagation`). Both build through `tools/scratch_build.sh`; to configure without
  network access, set `FETCHCONTENT_SOURCE_DIR_NLOHMANN_JSON` to the nlohmann_json sources. The
  final wave packet of two builds, or of the CPU and the GPU: `spodvr_dump_wave` and `spodvr_compare_wave` ([docs/performance.md](docs/performance.md)).
- **Production-grid validation:** beyond the gates, the full thesis calculation
  (`inputs/thesis_fig5_3_dz`, 7000 steps) was run with the Fortran, the C++ CPU and the C++ GPU
  code and compared state by state: [docs/validation_thesis_fig5_3_dz](docs/validation_thesis_fig5_3_dz/README.md).
  With the default (strict) GPU build, all 111 reaction probabilities match the Fortran's 7 printed digits.

## Running

The program reads one JSON input file:

```bash
./build/spodvr inputs/hardwall_small/input.json --threads 8 --output-dir runs/hardwall_small
```

This example (`inputs/hardwall_small`, 7000 time steps) takes about 3 minutes with 8 threads on a
desktop machine (AMD Ryzen 7 5800X: the long gate `Gate3FullRunLong`, which runs this input, took
164 s on 2026-11-15). Its set-up prints two
warnings, "Rotational basis too small to describe all open channels" and "the diffraction basis is
too small to describe all open channels": JMax = 4 and the 4 × 4 surface grid are kept small for
speed. They are harmless here because the hard wall depends only on Z and r, so the wave packet
stays in the initial rotational state and the specular diffraction channel.

Command-line options:

```
spodvr <input.json> [--threads N] [--schedule dynamic|static] [--device cpu|cuda] [--output-dir DIR] [--no-color]
       [--help] [--version]
```

- `--device cpu|cuda`: propagate on the CPU (default) or on an NVIDIA GPU, where the wave packet
  stays in GPU memory (needs a build with `-DSPODVR_USE_CUDA=ON`; see [docs/cuda.md](docs/cuda.md)).
  The results agree with the CPU's to rounding (≈1e-13, the level of a different FFT library).
- `--threads N`: number of OpenMP threads (default: the OpenMP default, i.e. `OMP_NUM_THREADS` or
  all hardware threads). The results do not depend on it. Use the
  number of physical cores: hyperthreads do not speed up these memory- and FFT-bound loops.
- `--schedule dynamic|static`: how the (Z, r) points are shared out over the threads. `dynamic`
  balances the uneven work best (one NUMA node, i.e. one CPU socket); `static` keeps every thread on
  the memory it first wrote (several NUMA nodes). Default: chosen from the number of NUMA nodes. The
  results do not depend on it.
- `--output-dir DIR`: where the data files are written (default: the current directory).

Exit codes (`src/app/run.cpp`, `src/app/main.cpp`):

- 0: success (also `--help` and `--version`);
- 1: any other error, e.g. an output file that cannot be opened;
- 2: an input error (the message names the key) or a wrong command line;
- 3: not supported: a feature that is not ported ("Not supported: …");
- 4: a checkpoint or BK file error: a restart or analysis-only file that is missing, corrupt or
  written with other parameters, or a checkpoint that cannot be written.

Output:

- **Screen:** an input summary, the grids, a memory estimate, the surface symmetry (with
  `UseSymmetry`), the gas-phase states, the initial wave packet, the analysis basis, a table of averages during the propagation (on a terminal with a live progress line:
  step, time per step, estimated remaining time), and compact probability tables (one row per
  energy) at the end. Colour and the progress line are off when the output is not a terminal (a
  pipe, `tee`, a file) or with `--no-color`.
- **Data files** (in the output directory, numbers in full precision):
  - the results: tab-separated tables (`.tsv`, energies in eV) with `#` comment lines that explain
    every column, then one header line; they load with `pandas.read_csv(f, sep="\t", comment="#")`;
  - `analysis.json`: the run metadata;
  - the `.dat` files: space-separated columns after `#` comment lines (atomic units, unless a column
    says otherwise); they load with `numpy.loadtxt(f)`;
  - the binary restart files.

| File | Content |
|---|---|
| `scattering_summary.tsv` | per energy: total scattering and reaction probability, P(v), average final energies |
| `s_matrix.tsv` | per energy and open final state (v, j, mJ, n, m): S-matrix element, probability, phase |
| `probabilities_rovib.tsv` | per energy and (v, j, mJ), summed over diffraction |
| `probabilities_diffraction.tsv` | per energy, (v, j, mJ) and diffraction order |
| `analysis.json` | run metadata: basis sizes, initial state, energy grid, file list; the SCM identity of a `CRPDistrtSurface` run |
| `averages.dat` | time (a.u.), norm, ⟨Z⟩, ⟨r⟩ (bohr), rest norm during the propagation |
| `asymptotic_potential.dat`, `gas_phase_states.dat` | set-up data: V_asy(r) (hartree and eV); the gas-phase energies (hartree) and wave functions |
| `checkpoint.bin`, `bk_coefficients.bin` | final wave packet and Balint-Kurti coefficients, for a restart |

With `alsoAnalysis` (or an analysis-only run, `onlyAnalysis`) the fine-analysis results are written
to `fine_analysis/` in the output directory. The file formats are described in
[docs/output_files.md](docs/output_files.md).

### Continuing a run (restart)

Set `"NewCalc": false`, `"tStart"` to the end time of the earlier run, `"nStep"` to the number of
steps to add and `"tEnd"` to tStart + nStep·dtStep (the input check requires this), and
`"RestartDir"` to the earlier run's output directory; use a *different* `--output-dir`. All other
parameters must stay the same, apart from output switches (the program checks the grids, masses,
time steps, wave packet and absorbers against the checkpoint). The S-matrix of the restarted run
covers the whole propagation; a split run gives bitwise the same results as one long run. Runs can
also be continued from the Fortran code's files (`"RestartFormat": "fortran"`: fort.8/fort.9 and
fort.22/fort.21) and written for it (`"CheckpointFormat": "fortran"` or `"both"`); see
[docs/input_schema.md](docs/input_schema.md#checkpoint--restart-keys-c-only).

### Running a oneMKL build from another location

Inside the build tree the oneMKL library path is built into the executable. If you copy or install
the executable elsewhere, make the oneMKL libraries findable first:

```bash
source /opt/intel/oneapi/setvars.sh        # or: export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib
```

Optional: `export MKL_CBWR=AVX2` gives bitwise identical results across machines with different
vector extensions (see docs/performance.md).

### Inputs

`inputs/` holds ready-to-run example inputs:

| Directory | Content |
|---|---|
| `inputs/hardwall_small/` | The hard wall reference case (H₂ on a hard wall, Cu(211)-like grid, 0.2–0.6 eV). `input.json` is the C++ input. `DVR.inp` is the same calculation for the Fortran reference code. |
| `inputs/scm_example/` | A technical example of the SCM potential (`CRPDistrtSurface`): one snapshot on a small grid, not physics; see its README. |
| `inputs/thesis_fig5_3_dz/`, `inputs/thesis_fig5_3_er/` | H₂ + Cu(111) on the Dai–Zhang (DZ; SystemID `H2Cu111DaiLight`) and the Persson et al. ER PES (`H2Cu111Persson`), the two curves of Fig. 5.3 of M. F. Somers' PhD thesis, at the grid sizes of the thesis (except a finer Z grid, see their READMEs; JMax 28, ψ ≈ 6 GB; hours per run). Each directory has a README with the parameter choices. The DZ case is the GPU benchmark. |

- The JSON keys mirror the Fortran `DVR.inp` field names. All keys, defaults, units and checks
  are described in [docs/input_schema.md](docs/input_schema.md).
- How to choose the Z grids, the wave packet, the analysis lines and the absorbers for a new
  input: [inputs/README.md](inputs/README.md).
- Unknown keys and inconsistent parameters are reported as errors.
- The test copies of these inputs (with the Fortran reference output) live in
  `tests/reference/`, together with the corrugated test cases (`cu111_dailight_*`,
  `cu111_persson_small`, and `crp_*_synthetic_small`, which must be run from the repository root
  because their `CRPlibDataSet` path is relative to it), whose `input.json` files can be run in
  the same way.
- The hard wall and test grids are deliberately tiny so that they run quickly. Production
  calculations use much larger grids: the thesis inputs above need 7.9 GiB (DZ) and 9.7 GiB (ER) of
  memory on the CPU (their set-up memory estimates:
  [validation record](docs/validation_thesis_fig5_3_dz/data/logs/cleaned_cpu_fast.txt),
  [inputs/thesis_fig5_3_er](inputs/thesis_fig5_3_er/README.md)); finer grids and larger rotational
  bases need more.

## Benchmark

Full Dai–Light run (`tests/reference/cu111_dailight_small`): 7000 time steps on a 176 × 56 × 8 × 8
(Z, r, x, y) grid with 15 rotational states, plus the on-the-fly analysis. AMD Ryzen 7 5800X
(8 cores), 8 threads, the three runs one after another; measured 2026-11-05 (GCC 13, Fortran 5.19),
not re-run since:

| Binary | Wall time | CPU time | Threads busy on average |
|---|---|---|---|
| Fortran | 680 s | 5255 s | 7.7 |
| C++, current flags | 764 s | 5040 s | 6.6 |
| C++, fast-math + native + PGO | 675 s | 4489 s | 6.7 |

The C++ code needed less CPU time than the Fortran code, but kept the threads less busy, so it was
slower in wall time: load imbalance in the fused T_xy + V loop, which has since been scheduled
dynamically (−9 to −13 % on the thesis grids; docs/performance.md). The fast-math build (then plain
`-ffast-math`) passed the regression gates except the hard-wall V<sub>asy</sub> table (by 1.6e-14, with
GCC 15; since NaN and inf are kept it passes that too), and its results differ from the default build: by ~1e-12 on the gate grids, by 7.8e-12 on the thesis grid (GPU builds)
([validation report](docs/validation_thesis_fig5_3_dz/README.md#strict-build-gpu-run-added-2026-11-08)).
It is therefore not the default.
Details in [docs/performance.md](docs/performance.md).

GPU (`--device cuda`), DZ case of `inputs/thesis_fig5_3_dz` as of 3c1d91b (128 × 40 × 12 × 12 grid, 435 rotational
states, 4.8 GiB wave function; the input now has 160 Z points and a 6.0 GiB wave function), NVIDIA RTX 4070 Ti (a GeForce card: FP64 at 1/64 rate) against
8 threads of the Ryzen 7 5800X, per time step:

| Build | Time per step | Speed-up |
|---|---|---|
| CPU, default flags | 3.6 s | 1× |
| GPU, `build-cuda` (no FMA; agrees with the CPU to ≈1e-14 here) | 0.61 s | 5.9× |
| GPU, `build-cuda-fast` (FMA; S matrix of the thesis run 7.8e-12 from `build-cuda`) | 0.54–0.57 s | 6.3–6.7× |

At these speeds, the 10 000 steps that the DZ input had then (tEnd 25 000 a.u.; it now has 7000)
take about 1.5–1.8 h on the GPU instead of about 10 h. On the
current 160-point grid the validation record's 7000-step runs took 0.73 s per step on the GPU
(strict build) and 4.24 s on the CPU (fast-math build)
([timing table](docs/validation_thesis_fig5_3_dz/README.md#performance-as-run)). Details
and the limits of this card in [docs/cuda.md](docs/cuda.md).

## Repository layout

| Path | Content |
|---|---|
| `src/` | the program: `input/`, `grids/`, `potential/`, `symmetry/`, `gasphase/`, `propagation/`, `analysis/`, `io/`, `fft/`, `linalg/`, `core/`, `report/`, `app/`, `cuda/` (optional GPU backend), `crplib/` (the CRP library port), `fortran_text/` (the Fortran list-directed reader shared by `crplib/` and the SCM code) |
| `tests/` | unit tests (one file per module) and regression tests against Fortran reference output (`tests/reference/`, `tests/crplib/reference/`, `tests/scm/reference/`) |
| `tools/` | Fortran reference scripts, benchmark and comparison tools, the docs consistency checker |
| `docs/` | input schema, output files, validation against the Fortran (`validation.md`, `fortran_deviations.md`), the upstream changes and their port status (`upstream_sync.md`), the flux analysis findings (`flux_analysis_findings.md`), the SCM potential (`scm.md`, `scm_open_questions.md`), the CRP library (`crp_library.md`), the Fortran bugs found (`fortran_input_bugs.md`, `fortran_fine_analysis_bugs.md`), performance notes, GPU backend |
| `inputs/` | example inputs |
| `potentials/` | data files of the Fortran potentials, one directory per system (verbatim copy of `potential.dir` of the upstream Fortran repository) |

## License

This program is free software, licensed under the GNU General Public License, version 3. See
[LICENSE](LICENSE).
