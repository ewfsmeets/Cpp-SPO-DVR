# physics Q1: with UseSwitching the asymptotic potential is not the gas phase

A bug of the Fortran SPO-DVR code (upstream r1727, git `3d00331`, and earlier versions as stated below), found in an
audit of the active build for errors that make a computed quantity wrong. The candidate was read by an auditor, checked
by an independent reviewer who saw the code but not the auditor's conclusion, and run where possible. Line numbers refer to
`sources.dir/SOURCES.dir/` at `3d00331`. The short entry is in [fortran_deviations.md](../fortran_deviations.md) under
the same ID; the other bugs of the audit are listed in [README.md](README.md).

**Status:** open in the Fortran (5.19 and r1727); not yet reported upstream (2026-11-25). The C++ does the same (it
reproduces the Fortran); when `ZmaxSwitch` is absent the two codes differ (below).

## The demonstration runs

- **Executable:** the reference install `r1727-3d00331` (upstream git `3d00331` built with `O3FLAGS = -O3`;
  [tools/fortran_reference/README.md](../../tools/fortran_reference/README.md#the-reference-install)); for the size
  of the effect, also a scratch build of the same commit with the same flags and the one-line fix proposed below
  (`H2SurfaceSystem.f90` line 426).
- **Inputs:** the `DVR.inp` of a reference case in `tests/reference/` with the keys of each case changed or appended
  (and its data files, as `tools/fortran_reference/reference_config.sh` `stage_inputs` stages them).
- **How:** each run in a new scratch directory through `reference_config.sh` (`run_clean`: a clean environment, 4
  OpenMP threads). The C++ runs used the program at git `48e5bdd`, 4 threads.
- A Fortran input error stops with exit status 0 ([input_checks.md](input_checks.md#the-demonstration-runs)).

## The Fortran code

`UseSwitching` switches the 6D potential smoothly to the gas phase between ZminSwitch and ZmaxSwitch, so that the
asymptotic region can start at ZmaxSwitch instead of at the PES's own asymptotic Z. `InputNew.f90 SetupInput` moves
the asymptotic Z (lines 66–71):

```fortran
    ! If the potential is switched, the new asymptotic value of Z is given by 'ZmaxSwitch'.
    ! First, store the REAL, i.e. unswitched, asymptotic value.
    zasymp_unswitched = zasymp
    Result = GetLogicalField( 'UseSwitching', lswitch )
    if ( lswitch ) Result = GetRealField( 'ZmaxSwitch', zasymp )
```

The 6D potential is stored for Z < zasymp only (`Setup.f90` lines 343–349: `nzpot` is the last point below
zasymp). The switch blends it towards the gas phase, the PES at the *unswitched* asymptotic Z
(`Potential.f90 switchPotentialToGasPhase`, lines 903–908 and 926–935):

```fortran
    ! Setup the asymptotic (gas phase) potential. Use the 'zasymp_unswitched' value because this
    ! tells us for what Z value the potential module will give us the gas phase potential. ...
    do ir = 1, nr
      asypot(ir) = getVH2Surface( 0.0, 0.0, zasymp_unswitched, rgrid(ir), Constants_PI/2.0, 0.0, nzpot+1, ir ) / chaev
    enddo
    ...
          vdvr(:,:,izr) = fswitch(iz) * (vdvr(:,:,izr) - asypot(ir)) + asypot(ir)
```

The asymptotic potential V_asy(r), however, is evaluated at the *switched* zasymp (`H2SurfaceSystem.f90` line 426,
`Potential.f90 setupTheAsymptoticPotential` line 302):

```fortran
      CALL setupTheAsymptoticPotential( getVH2Surface, zasymp )
...
      vasy(ir) = getVH2Surface( 0.0, 0.0, Z, rgrid(ir), Constants_PI/2.0, 0.0, nzpot+1, ir ) / chaev
```

V_asy is the potential of every Z ≥ ZmaxSwitch on both grids (`GasPhase.f90 SetupAsympWF` line 325 makes its
factor `XVGasPhase`, applied by `PerformEXPGasPhaseOnSpecGrid` and `PerformEXPGasPhaseOnScatGrid`, lines 353–408), and the
asymptotic vibrational-rotational states, their energies and so every channel energy of the analysis are computed in
it (`Main.f90` line 320, `call SetupAsympWF( vasy )`).

## What goes wrong

- Below ZmaxSwitch the wave packet feels a potential that tends to the gas phase; from ZmaxSwitch on it feels the
  raw PES at Z = ZmaxSwitch (at x = y = 0, θ = 90°, φ = 0). Unless the PES is already the gas phase there, V steps at
  ZmaxSwitch, and every asymptotic state is an eigenstate of the wrong potential.
- The switch exists for the case where it is not: a PES whose interaction reaches beyond the grid. The printout of
  the Fortran says so (`InputNew.f90` lines 1043–1045: "the Z-value at which the potential first reaches the gas
  phase, is the new asymptotic Z-value").
- The step lies at ZmaxSwitch, which in the reference cases is also the analysis line `Zana`.
- Size: the step is V_PES(ZmaxSwitch) − V_gas. For the Dai–Light PES (its own asymptotic Z is 15 bohr) switched
  between 6 and 8 bohr, as in the reference cases, it is +8.8 to +10.9 meV for r = 0.8 to 2.0 bohr (below). The
  vibrational levels lie about 10 meV too high, which shifts the energy axis of every result. For the SRP48 CRP data
  set switched at 11.3–11.9 bohr (its asymptotic Z is 12.28 bohr) the step is about 0.06 meV (read from the data
set's `extrapol.dat`; not run).

## Demonstration

**Case (a), the step.** `cu111_dailight_small` with `Potential : LOGICAL : TRUE;` (a potential-only run), and the same
case with `write2DPES : LOGICAL : TRUE;` and `nZ` 256, so that the cut V(r, Z) at x = y = 0, θ = 90°, φ = 0 reaches
Z = 18.4 bohr, where the Dai–Light PES is the gas phase. Both exit with status 0. The printed V_asy is the cut at
Z = 8 to every digit; the stored 6D potential at the last point below ZmaxSwitch (Z = 7.92, read from the `.pes` file
the run writes; switching factor 0.0039) is the gas phase to 0.12 meV at every x, y, θ, φ:

| r (bohr) | V_asy printed (eV) | gas phase, cut at Z = 16 (eV) | V_asy − gas (meV) | max \|V(Z = 7.92) − gas\| (meV) |
| --- | --- | --- | --- | --- |
| 0.80 | −1.283978 | −1.292766 | +8.788 | 0.115 |
| 1.10 | −4.117126 | −4.126367 | +9.241 | 0.114 |
| 1.40 | −4.735868 | −4.745668 | +9.800 | 0.112 |
| 1.70 | −4.401053 | −4.411451 | +10.398 | 0.108 |
| 2.00 | −3.729128 | −3.740037 | +10.909 | 0.100 |

**Case (b), the results.** The full `cu111_dailight_small` run (7000 steps) with the fix (V_asy evaluated at
`zasymp_unswitched`), compared with the reference output of the unmodified r1727 (`tests/reference/cu111_dailight_small`,
the same input):

- The stored 6D potential (`.pes`) of the two runs is identical byte for byte: the fix changes V_asy only.
- The lowest gas-phase level E(v = 0, j = 0) moves from −4.467141 eV to −4.477030 eV, 9.889 meV lower; every channel
  energy and so every collision energy of the analysis moves with it.
- The reaction probability (1 − the sum of the probabilities) at the same collision energies:

  | E_trans (eV) | P_reaction r1727 | P_reaction with the fix | sum of the probabilities r1727 | with the fix |
  | --- | --- | --- | --- | --- |
  | 0.20 | 6.58e-03 | 2.42e-04 | 0.993420 | 0.999758 |
  | 0.22 | −8.46e-03 | 7.09e-04 | 1.008464 | 0.999291 |
  | 0.30 | −6.56e-04 | 5.93e-04 | 1.000656 | 0.999406 |
  | 0.40 | −3.44e-04 | −9.60e-05 | 1.000344 | 1.000096 |
  | 0.50 | 5.23e-04 | 8.43e-04 | 0.999477 | 0.999157 |
  | 0.55 | 5.19e-03 | 3.83e-03 | 0.994805 | 0.996165 |
  | 0.59 | 3.11e-02 | 2.14e-02 | 0.968867 | 0.978642 |
  | 0.60 | 4.49e-02 | 3.15e-02 | 0.955090 | 0.968521 |

- Two effects show. Where P_reaction rises steeply (0.55–0.60 eV) the fixed curve is the r1727 curve shifted by
  about 10 meV (r1727 at 0.59 eV ≈ the fix at 0.60 eV), so P_reaction there is 30 % too high. At low energies, where
  the fixed run's P_reaction is 1e-4 to 1.5e-3, the r1727 sums swing between 0.9934 and 1.0085 (P_reaction negative at
  16 of the 41 energies), the fixed run's stay within 0.0015 of 1 below 0.5 eV: the step at the analysis line disturbs the
  flux the analysis measures there. Largest |ΔP_reaction| over the 41 energies: 1.34e-2.
- This is a short test grid (`Zana` = ZmaxSwitch = 8 bohr, the Dai–Light PES 9.8 meV from the gas phase there); the
  size scales with V_PES(ZmaxSwitch) − V_gas of the PES and the switch used.

**Without `ZmaxSwitch`.** `ZmaxSwitch` is not pre-registered (`InputNew.f90` line 151 is commented out) and
`inputfields.f90 GetRealField` leaves its argument unchanged when the key is absent, so with `UseSwitching` and no
`ZmaxSwitch` the Fortran keeps the PES's asymptotic Z while the switch uses the default 8.0 (`globaldata.f90`). That
combination is consistent (V_asy is the gas phase). `cu111_dailight_small` without its `ZmaxSwitch` line stops with
`ERROR: ZAna < ZAsymptotic.` (the asymptotic Z printed as 0.15000000E+02).

## The C++ behaviour

- The same as the Fortran: `switch_to_gas_phase` blends towards the PES at `z_asymptotic_unswitched`
  (`src/potential/potential_grid.cpp` line 37), `compute_asymptotic_potential` evaluates V_asy at `z_asymptotic`
  = ZmaxSwitch (line 141). It keeps the Fortran's behaviour until upstream fixes it (the project's decision,
  2026-11-25): gates 1–4 compare the two codes on the eleven reference cases that use `UseSwitching`, so the fix goes
  into both codes at once and those references are regenerated from the fixed upstream version.
- Without `ZmaxSwitch` the C++ differs: `build_config` sets `z_asymptotic` to `z_max_switch`, whose default is 8.0
  (`src/input/process.cpp`, `src/input/config.hpp`), so the case above runs (exit code 0, "Z asymptotic 8.000000 bohr")
  with the step of case (a). A deviation, recorded under physics Q1 in [fortran_deviations.md](../fortran_deviations.md).

## Proposed fix

Evaluate V_asy where the switched potential arrives, the PES's own asymptotic Z
(`H2SurfaceSystem.f90` line 426):

```fortran
      CALL setupTheAsymptoticPotential( getVH2Surface, zasymp_unswitched )
```

Without switching `zasymp_unswitched = zasymp`, so nothing else changes. `nzpot`, the checks of Zana, Z0 and
ZSpOptStart, and the switch keep the switched zasymp.
