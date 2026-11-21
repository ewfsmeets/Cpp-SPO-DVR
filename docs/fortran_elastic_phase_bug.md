# Bug: the phase of the specular S-matrix elements is referred to the wrong Z

**Status:** open. The bug is in the Fortran (5.19, 5.24 and upstream master) and in the C++ port, which reproduces
it. It was found in a code review on 2026-11-20 and verified on the hard-wall case against the exact answer.

**In one paragraph.** With two grids (`ProjOp = 2`), the S-matrix elements of every channel are extracted from
Balint-Kurti (BK) coefficients taken at the analysis line `Zana`. The exception is the initial rotational-diffraction
state (the initial j, mJ and diffraction order (0, 0), for every v). That state lives on the specular grid and is
measured at a different line, Z_sp. The S-matrix formula then removes the free outgoing phase e^{ikZ} at `Zana` for
every channel, this one included. So the reported elements of that state carry an extra phase:

```text
S_reported = S_true · exp(+i k_f (Z_sp − Zana))        k_f: the channel's own final normal momentum
```

|S|, and therefore every probability, is unaffected. The phase (Re S, Im S, `phase`) of the specular elements is
wrong by k_f (Z_sp − Zana), which is several radians in every reference case, and it changes from run to run with
Z_sp.

## Contents

- [Background: how S is extracted](#background-how-s-is-extracted)
- [The bug in the code](#the-bug-in-the-code)
- [What it causes](#what-it-causes)
- [Evidence](#evidence)
- [History](#history)
- [Proposed fix: Fortran](#proposed-fix-fortran)
- [Proposed fix: C++](#proposed-fix-c)
- [Open points](#open-points)

## Background: how S is extracted

- **The BK method.** The wave packet starts at Z0 as (2α/π)^¼ exp(−α(Z−Z0)² − i k_c Z), a superposition of incoming
  waves with weights η(k). Take a stationary scattering state with outgoing part S_f(k) e^{+i k_f Z} in channel f,
  and project the packet onto channel f at a line Z_L. The time-to-energy transform of that coefficient c_f(t) is:

  ```text
  ∫ c_f(t) e^{iEt} dt  ∝  η(k) · S_f(E) · e^{+i k_f Z_L}
  ```

  To obtain S_f, the code divides out |η|, the phase of η (the `(rkz − sk0)·z0` term) and e^{+i k_f Z_L}
  (the `−skz·ZAnalysis` term). The `−π` makes S = −R for a pure reflection R. The last step is only right if the
  coefficient was taken at Z_L = `Zana`.
- **The two grids.** With `ProjOp = 2` the initial rotational-diffraction state is kept on its own reduced (Z, r)
  grid, the specular grid, which extends further out in Z. In every potential step the code adds it into the
  scattering grid for the interaction and projects that channel back onto the specular grid afterwards
  (`TwoGrids.f90` `AddSpecToScat` / `SubSpecFromScat`; C++ `src/propagation/potential_operator.cpp:115-135`).
  The scattering grid's component of that channel is therefore zero after every potential step, and the reflected
  elastic and vibrationally inelastic specular waves travel out on the specular grid. That this makes the
  scattering-grid component zero to rounding everywhere, also beyond the potential region, is inferred from the code
  and not measured.
- **The specular analysis line.** The incoming packet also lives on the specular grid. So the specular channel cannot
  be analysed at a line the incoming packet has yet to pass. The line Z_sp is placed:
  - with a fixed line (`ZAnaSp`, or the start of the specular absorber `ZSpOpt`), there;
  - with the dynamic line (the default), it moves inward behind the incoming packet and is fixed when the reflected
    part comes back (`FindNewAnalysisPoint`; C++ `AnalysisSurface::update`). It ends anywhere between `Zana` and the
    start of `ZSpOpt`, depending on `TolSp`, the grid, the time step and the dynamics.

  From the step at which the line is fixed, the BK coefficient of the initial state is *replaced* by the projection
  of the specular wave at Z_sp (a documented kept quirk, [fortran_deviations.md](fortran_deviations.md)).

## The bug in the code

The coefficient is taken at Z_sp, but the phase is removed at `Zana`.

### Fortran

- 5.19 `Analyzor.f90`; upstream master (dce0ecc) line numbers in brackets.
- `BalintKurtiCoeffs`, `Analyzor.f90` lines 380–389 [393–401]:

  ```fortran
  if ( FixedZAnalysisSpec .AND. (nproj == 2) ) then
     do iv = 1 , nviba
       bcoef = (0.0,0.0)
       do ir = 1 , nr
         bcoef = bcoef + SpecWave(ir,izanlySpecChan) * eivea(ir,jin,iv)   ! measured at zgrid(izanlySpecChan) = Z_sp
       enddo
       balcf(inchna,iv) = bcoef * delr
     enddo
  endif
  ```

- `smatrix_calc`, `Analyzor.f90` line 635 [647]. On-the-fly analysis, `nread` 0 and 1. This line runs for every
  channel, `is == inchna` included:

  ```fortran
  arg = - skz * ZAnalysis + ( rkz - sk0 ) * z0 - pi
  ```

- `probaf`, `Analyzor.f90` line 1446 [1397]: the fine analysis (`nread = 2`) does the same with `fskz(idat,is,iv)`.

### C++

On `main`, unchanged since the port:

- `src/analysis/balint_kurti.cpp:74-82`: `if (surface.fixed())`, `bk(a.initial_state(), v)` = the projection of
  `specular(ir, surface.index())`.
- `src/analysis/s_matrix.cpp:106` (`SMatrixAccumulator::result_at`):
  `const double arg = -kz * p.z_analysis + (d.kz - k0) * p.z0 - pi;` for every channel.
- `src/analysis/fine_analysis.cpp:181` (`fine_spectrum`): the same line.

The CPU and GPU backends share this analysis code (the backend only supplies the specular wave and the `Zana`
plane), so both are affected the same way.

**The size of the error.** The extra factor is e^{+i k_f (Z_sp − Zana)}, where k_f = √(2M(E − ε_f)) is the
channel's own final normal momentum. For the elastic channel this is `sk0`, not the packet's centre momentum `rkz`
and not the total momentum: the oblique-incidence run below tells these apart.

## What it causes

- **Wrong:**
  - Re S, Im S and `phase` of every element of the initial rotational-diffraction state: the elastic element and the
    vibrationally inelastic specular elements (v ≠ v_in, same j, mJ and diffraction order).
  - Where they are written: Fortran stdout and fort.25 (on-the-fly and fine analysis); C++ `s_matrix.tsv`,
    `analysis.json` and the fine analysis's spectrum files.
- **Size:** k_f (Z_sp − Zana) radians.

  | Case (Fortran reference) | `Zana` (bohr) | Z_sp (bohr) | Offset per unit k_f | At 0.4 eV normal incidence (k ≈ 10.4 a.u.) |
  | --- | --- | --- | --- | --- |
  | `hardwall_small` | 8.0 | 8.64 | 0.64 rad | 6.65 rad |
  | `cu111_dailight_small` and its `_short` siblings | 8.0 | 10.88 | 2.88 rad | ≈ 30 rad |
  | `scm_cu111_srp48_small` | 12.05 | 14.30 | 2.25 rad | ≈ 23 rad |

  Z_sp is the "specular channel fixed at Z" line of each case's `fortran_stdout.txt`. Every propagating reference
  case is affected; only a line that happens to fix exactly at `Zana` is not.
- **Not reproducible between runs.** With the dynamic line, Z_sp depends on `TolSp` and the run. On the hard wall,
  Z_sp = 8.64 with `TolSp` 1e-7, 10.80 with `TolSp` 1e-11 and 14.96 with `ZAnaSp` 15.04. Two runs of the same physics
  therefore report elastic phases that differ by an arbitrary amount (up to π after wrapping), while all other
  channels keep their phase.
- **Derived quantities:**
  - A time delay from the phase, τ = d arg S / dE, gains a spurious M (Z_sp − Zana) / k_f: ≈ 226 a.u. (5.5 fs) on
    the hard wall and ≈ 1000 a.u. (25 fs) on Dai-Light, both at 0.4 eV.
  - Interference between the specular and other channels, and comparisons of Re S / Im S with other methods (e.g.
    time-independent close coupling, or a wave-packet code that analyses the specular channel at another line), are
    off by the same phase.
- **Not affected:**
  - |S|, all probabilities, their sums and averages, the reaction probability;
  - the BK coefficients themselves (`bk_coefficients.bin`, fort.21): they are correct coefficients at their own lines;
  - gates 1–4 (gate 4 compares probabilities and BK coefficients, no phases);
  - the plot tool's time delays (`tools/plots/compute/time_series.py` `time_delays`). They are computed from the BK
    coefficients with each channel's own line Z_f, so they are already right.
- **Gate 5 compares phases.** It compares fort.25 Re S / Im S / phase with the C++, and both carry the same wrong
  reference plane, so it passes. After a fix it must keep comparing with the Fortran's convention (see the C++ fix
  below).

## Evidence

### 1. The Fortran's own reference output, against the exact answer

The hard wall is V(r) + V_Z(Z) with no
coupling, so the elastic element is exactly the 1D reflection amplitude of V_Z: S = −R, |R| = 1, and below 0.6 eV
no other channel is open (ΣP(other) ≤ 2e-9 in the runs).

- The script below computes R by Numerov integration of V_Z (`HardWallPES.f90`). It compares arg R with the elastic
  S that the Fortran printed in the committed reference fort.25 (`tests/reference/hardwall_small/fortran_fort25.txt`,
  5.24 r1713 `-O3`; Z_sp = 8.64 from `fortran_stdout.txt`).
- Run from the repository root:
  `.venv/bin/python -I -B elastic_phase_check.py tests/reference/hardwall_small/fortran_fort25.txt`.
- Output, 41 energies from 0.20 to 0.60 eV (excerpt):

  ```text
  E_z/eV  k*(Zsp-Zana)/rad  arg S - arg(-R)  arg S - arg(-R e^{ik(Zsp-Zana)})
  0.200    4.703            -1.5898          -0.0101
  0.400    6.652            +0.3650          -0.0035
  0.600    8.147            +1.8589          -0.0045
  max |arg S - arg(-R)|                   = 1.8589 rad
  max |arg S - arg(-R e^(ik(Zsp - Zana)))| = 0.0101 rad
  ```

- The plain difference equals k (Z_sp − Zana) wrapped to (−π, π], e.g. 4.703 − 2π = −1.580 at 0.2 eV. After that
  factor is removed, 0.010 rad is left, which is the time-step error of the propagation (next point).

<details>
<summary>elastic_phase_check.py</summary>

```python
"""Hard wall: the Fortran's elastic S phase (fort.25) against the exact 1D reflection amplitude.
Usage: python3 elastic_phase_check.py tests/reference/hardwall_small/fortran_fort25.txt"""

import re
import sys

import numpy as np

CHAEV, MASS = 27.2114, 2 * 1.0078 * 1822.89  # legacy Ha->eV; H2 mass in m_e
Z_ANA, Z_SP = 8.0, 8.64  # DVR.inp Zana; "ZAnalysis on specular channel fixed at Z = 0.864E+01"


def reflection(k, h=2.0e-4):
    """R of psi ~ e^{-ikZ} + R e^{ikZ} for V_Z of the hard wall (HardWallPES.f90), by Numerov."""
    z = np.arange(0.30, 7.2, h)
    vz = 1.0 / (z**4 + 0.001)
    switch = np.where(z < 0.1778, 1.0, 0.5 + 0.5 * np.cos(np.pi * (z - 0.1778) / (7.0 - 0.1778)))
    f = 2 * MASS * np.where(z < 7.0, switch * vz, 0.0) / CHAEV - k * k
    w, psi = 1 - h * h * f / 12, np.zeros(z.size)
    psi[1] = 1e-30
    for i in range(1, z.size - 1):
        psi[i + 1] = ((12 - 10 * w[i]) * psi[i] - w[i - 1] * psi[i - 1]) / w[i + 1]
        if abs(psi[i + 1]) > 1e200:
            psi *= 1e-200
    i1, i2 = z.size - 400, z.size - 1  # both beyond Z = 7, where V_Z = 0
    a = np.array([[np.exp(-1j * k * z[i]), np.exp(1j * k * z[i])] for i in (i1, i2)])
    c_in, c_out = np.linalg.solve(a, np.array([psi[i1], psi[i2]], dtype=complex))
    return c_out / c_in


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


text = open(sys.argv[1]).read().replace("D", "E")
worst = {"plain": 0.0, "shifted": 0.0}
print("E_z/eV  k*(Zsp-Zana)/rad  arg S - arg(-R)  arg S - arg(-R e^{ik(Zsp-Zana)})")
for block in text.split("total energy =")[1:]:
    ez = float(re.search(r"trans energy in z =\s*(\S+)", block).group(1))
    m = re.search(r"v j mj dif.order are:\s+0\s+0\s+0\s+0\s*\n.*\n\s*0\s+0\s+(\S+)\s+(\S+)", block)
    s = complex(float(m.group(1)), float(m.group(2)))  # elastic Re S, Im S
    k = np.sqrt(2 * MASS * ez / CHAEV)
    r = reflection(k)
    plain = wrap(np.angle(s) - np.angle(-r))
    shifted = wrap(np.angle(s) - np.angle(-r * np.exp(1j * k * (Z_SP - Z_ANA))))
    worst = {"plain": max(worst["plain"], abs(plain)), "shifted": max(worst["shifted"], abs(shifted))}
    print(f"{ez:.3f}   {k * (Z_SP - Z_ANA):6.3f}            {plain:+.4f}          {shifted:+.4f}")
print(f"max |arg S - arg(-R)|                   = {worst['plain']:.4f} rad")
print(f"max |arg S - arg(-R e^(ik(Zsp - Zana)))| = {worst['shifted']:.4f} rad")
```

</details>

### 2. C++ runs of the hard wall

- Set-up: variants of `hardwall_small`'s input; s_matrix.tsv against the same exact R; 0.2–0.6 eV.
- The residual after the predicted factor does not depend on where the line ends up, while the plain difference
  follows Z_sp:

  | Run | Z_sp | max \|arg S − arg(−R)\| | max residual after e^{ik(Z_sp − Zana)} |
  | --- | --- | --- | --- |
  | dt 4 / 2 / 1 / 0.5 | 8.64 | 1.85–1.86 rad | 2.8e-2 / 1.0e-2 / 5.6e-3 / 4.4e-3 rad (the time-step error, O(dt²)) |
  | `TolSp` 1e-11, longer specular absorber | 10.80 | 3.07 rad | 6.7e-3 rad |
  | fixed line `ZAnaSp` 15.04 | 14.96 | 3.14 rad | 8.9e-2 rad (rest of the packet not yet across the line at tEnd) |
  | oblique incidence, E∥ = 0.05 eV | 8.64 | 1.86 rad | 6.4e-3 rad (k_f = the normal momentum; the total k would leave ≈ 0.5 rad) |
  | j = 2, mJ = 1 initial state | 8.64 | 1.86 rad | 6.3e-3 rad |

### 3. The derivation

- An independent reviewer re-derived the S formula from the packet's momentum expansion. They confirmed the factor,
  the sign and k_f, found no other Z-dependent phase in the chain, and found no special handling of the initial
  state anywhere in the S-matrix code.
- The momentum shift acts in x, y only. The specular wave stores χ(Z, r) with exp(i k∥·R) factored out. The two Z
  grids share their origin.

### 4. A possible observation in the literature (inferred, not verified)

- In §4.3 of Somers' thesis, the split-operator calculation analyses the specular channel at Z∞ = 12.90 a0 and the
  other channels at 8.55 a0.
- Its Table 4.8 compares the result with a method that analyses at 8.55 a0. Only the specular element differs, and
  only in phase: |S| 0.1101 vs 0.1097, phase gap 2.96 rad. The thesis attributes the gap to the split-operator
  propagator.
- k (Z_sp − Zana) at 0.303 eV (k = 9.046) is 1.65 rad for ΔZ = 4.35 bohr, and 3.007 rad for ΔZ = 4.50 (one 0.15 bohr
  grid step further out). The latter is within 0.04 rad of the table.
- That thesis code is not this source, so this is a reading, not a proof.

## History

- The comment in `SetupAnalysis` (`Analyzor.f90` lines 170–175, 5.19) records a change of 14 Aug 2021. Before it, a fixed
  specular line was set to `izanly` (the `Zana` point), where the phase reference is right, but the line then sat in
  front of the initial packet:

  ```fortran
  ! this is where the stupid bug was. it put the ZSpAna in front of
  ! the initial WP where Zana is on the scat. grid [but the smatrix
  ! code does not correct for this incomming WP spec channel].
  ! So MFSomers commented out the boo boo on 14 aug 2021.
  ! izanlySpecChan = izanly
  ```

- The fix moved the line off `Zana` but left `smatrix_calc` referring the phase to `ZAnalysis`.
- The dynamic line can fix anywhere in [`Zana`, start of `ZSpOpt`], so it has the same offset. When the dynamic line
  was introduced was not checked.
- A commented-out block after `smatrix_calc`'s `arg` line adds the eq. 26a term of J. Chem. Phys. 103, 1216 (1995) for a line
  at `ZAnalysis` (`corrTerm = exp(-2i sk0 ZAnalysis)`): the code once handled the specular channel at `ZAnalysis`.

## Proposed fix: Fortran

In `smatrix_calc` (`nread` 0/1), remove the free phase at the line where the coefficient was taken. The initial
rotational-diffraction state is `is == inchna` for every `iv`, as in `BalintKurtiCoeffs`:

```fortran
               skz = sqrt(skzsq)

               ! The initial rot-dif state (all iv) is measured on the specular grid at
               ! zgrid(izanlySpecChan) once that line is fixed (BalintKurtiCoeffs):
               ! remove the free outgoing phase exp(i skz Z) there, not at ZAnalysis.
               zRef = ZAnalysis
               if ( nproj == 2 .and. FixedZAnalysisSpec .and. is == inchna ) zRef = zgrid(izanlySpecChan)

               arg = - skz * zRef + ( rkz - sk0 ) * z0 - pi
```

- `zRef` is `REAL` through the routine's `implicit real(8) (a-h,o-z)`, i.e. double, as `ZAnalysis` and `zgrid` are.
- `smatrix_calc` runs at the end of the run (and of a restarted run), when `FixedZAnalysisSpec` and `izanlySpecChan`
  hold their final values. A line that never fixed leaves the coefficient at `Zana`, so `zRef = ZAnalysis` is right
  there too.
- `probaf` (`nread = 2`, line 1446 [1397]) needs the same change with `fskz(idat,is,iv)`. It also needs Z_sp in an
  analysis-only run: check whether `nread = 2` reads `izanlySpecChan` (fort.8); otherwise take it from an input key.
- Coefficients from before the line is fixed are taken from the scattering grid at `Zana`. That component is zero by
  construction (see [Background](#background-how-s-is-extracted)), so one reference plane for the whole sum is exact
  up to rounding. That is inferred from the code: a check is to print `balcf(inchna,:)` before fixing.
- Test: the script above, applied to the fort.25 of the fixed code, must give ≈ 0.01 rad in the plain column
  (max |arg S − arg(−R)|) instead of 1.86 rad.

## Proposed fix: C++

The fix is a deliberate deviation from the Fortran. It changes only the phases of the initial state's elements.

1. **S matrix (on-the-fly analysis).** `SMatrixAccumulator::results` / `result_at` (`src/analysis/s_matrix.cpp`)
   take the line where the specular coefficients were measured, `std::optional<double> z_specular_line`. The caller
   (`src/app/run.cpp:360`, `propagation.scattering_results()`) passes `z_specular[surface.index()]` when
   `surface.fixed()`, else `std::nullopt`. In the channel loop:

   ```cpp
   // The initial rot-dif state is measured on the specular grid at Z_sp once that line is fixed
   // (balint_kurti.cpp); remove the free outgoing phase there. The Fortran uses Zana for every channel.
   const bool on_specular_line = z_specular_line && is == a.initial_state();
   const double z_measured = on_specular_line ? *z_specular_line : p.z_analysis;
   const double arg = -kz * z_measured + (d.kz - k0) * p.z0 - pi;
   ```

   A restart restores the surface state from the checkpoint (native, or the Fortran files), so a split run gets the
   same Z_sp.
2. **Fine analysis.** `fine_spectrum` (`src/analysis/fine_analysis.cpp:181`) gets the same parameter. An
   analysis-only run must know Z_sp of the run that wrote the BK file:
   - the native files: read the analysis surface (index, fixed) from the header of `checkpoint.bin` in `RestartDir`,
     which the plot tool already reads for the specular line. Read only the header, not ψ (6 GiB on the DZ grid);
   - the Fortran files: the surface record of fort.9 (`read_fortran_wave_function` already returns it);
   - if neither is there, stop with an error naming the missing file (the input rules: never continue silently).
     Alternatively, record Z_sp in `bk_coefficients.bin`, which changes the file format.
3. **Gate 5 keeps comparing with the Fortran.** Its test calls `fine_spectrum` with `std::nullopt` (every channel
   referred to `Zana`, the Fortran convention). That needs no production mode and no input key, and the gate's
   tolerances stay as they are. The alternative, correcting the Fortran's printed Re/Im by the factor, would mix the
   rounding of Re and Im and loosen the tolerance. No other gate compares phases.
4. **Tests:**
   - Unit (`tests/test_analysis.cpp`): synthetic BK sums with Z_sp ≠ `Zana`.
     - The initial state's S is the `nullopt` result times exactly e^{−ik_f (Z_sp − Zana)}, for every v.
     - Every other state is unchanged.
     - |S| is unchanged to rounding.
   - Physics (a new deviation needs a test): the hard-wall elastic phase from the existing `hardwall_small` gate run
     against arg(−R) from a small Numerov solver in the test, tolerance ≈ 1.5e-2 rad. The fixed code measured 1.0e-2
     on the Fortran's fort.25 and on C++ runs at the gate's time step, so the margin is small: set the tolerance from
     the measured worst case.
   - Negative test: the `nullopt` phase must fail. Assert the number of failing energies: all 41 except those where
     k (Z_sp − Zana) lies within the tolerance of a multiple of 2π (near 0.37 eV on this case).
5. **Docs:**
   - [fortran_deviations.md](fortran_deviations.md): a new entry with the next free `analysis` ID, citing this file
     and the test.
   - [output_files.md](output_files.md): the `phase` column of `s_matrix.tsv` refers every element to the line where
     it was measured.
   - [validation.md](validation.md): gate 5 uses the Fortran convention.
   - The docstring of `tools/plots/compute/time_series.py` `time_delays` ("s_matrix.tsv refers every phase to
     Zana").
6. **Checks:**
   - Not a bitwise refactor. `tools/check_bitwise_refactor.sh` must show differences only in the S values and phases
     of the initial state's rows (`s_matrix.tsv`, `analysis.json`, the fine-analysis spectrum files and the printed
     tables). P may move in the last bits.
   - The full suite, CPU and CUDA: the GPU shares the analysis code, so no CUDA change is needed.

Effort: small (two formula lines, one parameter, a header read for the fine analysis, tests and docs). Risk: low;
phases only, and gate 5 keeps its Fortran comparison.

## Open points

- **Fix, or only document?** A fix is a deviation from the Fortran. Documenting alone (in output_files.md) leaves
  run-dependent elastic phases in the output.
- **Report it upstream?** The Fortran is affected in 5.19, 5.24 and master; the patch above is a few lines.
- **The scattering-grid part before fixing.** Its being zero is inferred from the projection in every potential
  step and not measured. If it were not, apply the factor per coefficient during the time accumulation, only to
  coefficients taken after fixing (a separate sum for the initial state).
- **The thesis Table 4.8 reading** (Evidence point 4) is a hypothesis.
