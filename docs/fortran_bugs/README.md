# Fortran bugs

Bugs of the Fortran SPO-DVR code found while porting it and in later audits, each with the Fortran code, what goes
wrong, a demonstration with an unmodified Fortran build and what the C++ does. The short entries, one per ID, are in
[fortran_deviations.md](../fortran_deviations.md); the upstream fixes in [upstream_sync.md](../upstream_sync.md).

| Document | Contents | Status in the Fortran (r1727) |
| --- | --- | --- |
| [physics_q1_switching.md](physics_q1_switching.md) | physics Q1: with `UseSwitching` the asymptotic potential is the PES at ZmaxSwitch, not the gas phase | open |
| [physics_q2_restart.md](physics_q2_restart.md) | physics Q2: a restart is not checked against the run it continues (dtAna, tStart, analysis basis) | open |
| [physics_q3_morphing.md](physics_q3_morphing.md) | physics Q3: PES morphing changes V where nothing is morphed, compounds `cons`, runs out of bounds with `UseSymmetry` | open |
| [input_checks.md](input_checks.md) | input checks that are wrong: flux Q1/Q2, absorber D1, low-Z D1–D4, analysis D1, analysis D4, input Q1–Q3 | most fixed in r1714–r1719; analysis D4 and input Q1–Q3 open |
| [fine_analysis.md](fine_analysis.md) | the interpolation of the fine analysis (ATAN instead of ATAN2, signs, `s1_real`) | `s1_real` fixed in r1713; the rest open |
| [elastic_phase.md](elastic_phase.md) | the phase of the initial state's S-matrix elements referred to the wrong Z (analysis D2) | fixed on the fly in r1720; open in the fine analysis |
| [analysis_limits.md](analysis_limits.md) | the highest analysed j for odd JIn, and VIn against VAna | fixed in r1726/r1727 |

## The audit of 2026-11-25

An audit of the Fortran reference (r1727, git `3d00331`) for errors that make a computed quantity wrong: the active build
only (`InputNew.f90`, the FFTW path, the ported PES modules; not the MPI code, `*_hack.f90` or the flux analysis).
Four auditors read the code area by area; every candidate was then checked by an independent reviewer who saw the code
and a neutral question, not the auditor's conclusion, and run with the r1727 install where possible. Found:

| ID | Severity | File |
| --- | --- | --- |
| physics Q1 | medium | [physics_q1_switching.md](physics_q1_switching.md) |
| physics Q2 | medium | [physics_q2_restart.md](physics_q2_restart.md) |
| physics Q3 | medium (morphing only) | [physics_q3_morphing.md](physics_q3_morphing.md) |
| input Q1 | low (the run aborts) | [input_checks.md](input_checks.md#input-q1-an-r-absorber-from-the-first-r-point) |
| input Q2 | low | [input_checks.md](input_checks.md#input-q2-nphisoptimal-with-a-c3v-or-c3-potential-stops-the-run) |
| input Q3 | low | [input_checks.md](input_checks.md#input-q3-the-time-checks-reject-valid-decimal-steps) |

### Candidates that are not bugs

Read and reviewed in the same audit; recorded so that they need not be re-derived.

| Where | Candidate | Verdict |
| --- | --- | --- |
| `GasPhase.f90` lines 310–316, `Analyzor.f90` 868–871, 1535–1539 | The mean rotational energy uses E(0, j) − E(0, 0) for every v, so the mean vibrational and rotational energies do not add up to the mean internal energy for v ≥ 1 | quirk: a definition, not documented; the C++ the same (`src/gasphase/gas_phase.cpp`, `src/analysis/s_matrix.cpp`) |
| `GridEnergies.f90` lines 105–113 | Half-step kinetic factors as `sqrt` of the full-step ones: −exp(−iEdt/2) when E·dt mod 4π lies in (π, 3π) | quirk: only components with E·dt > π change sign, at both ends of each analysis step; the reference grids reach 0.48 rad |
| `Potential.f90` lines 1350–1361, `Absorption.f90` | `ZOptLow…` acts on the scattering grid only, so not on the initial channel | by design: each grid has its own absorbers (`ZSpOptLow…` for the specular grid) |
| `PESForCu111Persson.f90` line 57 | The Persson PES is not the gas phase at its asymptotic Z = 8 bohr: V_asy is +1.03 meV above the far field at r = 1.4 bohr (+0.89 to +1.71 meV for r = 0.8–3.05; demonstrated with a potential-only run and a cut to Z = 38.9 bohr) | quirk: a property of the PES model, consistent within a run |
| `CRPlibrary.dir/Spline.f90` DPL, DPLP | The linear extrapolation beyond the end knots does not continue the clamped spline's slope (a kink of m·e/3) | quirk: reached only by an atom cut whose z grid is narrower than cut 1 (not SRP48) |
| `CRP_C{3v,6v,4v}_MoleculeSurface.f90` | `Rmin`, `Rmax`, `Zmin` initialised in the declaration (implicitly SAVEd) | unreachable: one data set per process |
| CRP data sets | Not checked: the atom cuts' z grids (C3v, C6v), the cut order in CRPData, the Z and r coverage of `extrapol.dat` and `VH2.dat` | the SRP48 data set satisfies all; the C++ the same |
| `MyOperators.f90` `.SEQ.` | Characters 32 codes apart compare equal | harmless: no two names of the active build collide |
