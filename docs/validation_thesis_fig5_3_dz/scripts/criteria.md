# Criteria (pre-registered)

Written on 2026-11-07, before any comparison of the three runs (only the run logs' layout had
been looked at, to build the inventory). Sources: docs/validation.md (gate tolerances),
`tests/analysis_reference.hpp` (gate 4 code), `tests/test_regression_pes.cpp` (gate 1/2 code),
docs/fortran_deviations.md, docs/performance.md and docs/cuda.md (fast-math measurements).
Later changes go to "Criteria changes" at the end, dated, with the reason.

## Criteria

### Verdict words (fixed list)
| Word | Meaning |
|---|---|
| **identical** | bitwise equal (or equal integers / text after index conversion) |
| **rounding** | full-precision data that differ, within the full-precision criterion below |
| **printed** | within the printed precision of the Fortran value (criterion below) |
| **explained** | outside its criterion, with a quantitative explanation (predicted size, or the difference reproduced) |
| **FINDING** | outside its criterion and not explained |
| **shared code** | CPU↔GPU rows where both runs execute the same host code (set-up, potential, analysis); the equality is expected, not a test of the GPU kernels |
| **n/a** | not comparable (reason given) |

### Integers, flags, decisions, channels
- Exact after index conversion (C++ 0-based, Fortran 1-based; counts equal).
- Open/closed channels: exact. Decisions at a threshold: report the margin where an output
  gives it (e.g. the specular analysis line: the C++ log gives Z and t, the Fortran only Z; the
  Fortran step is pinned through fort.21).

### Printed Fortran values (C++ value at full precision against the Fortran text)
- |Δ| ≤ 1.01 · ½ unit of the last printed digit, floor 1e-20 (`printed_text_tolerance`).
- P_reaction: |Δ| ≤ max(1.01 · ½ printed unit, 1e-11) absolute (gate 4 floor; P_reaction = 1 − ΣP
  is known only to the absolute accuracy of the sum).
- Printed values whose half unit is below the predicted absolute rounding of that quantity
  (see "Prediction") may exceed the printed criterion; they count as **explained** only if
  |Δ| stays within that predicted size, which for a probability P = |S|² is
  |ΔP| ≤ 2 √P · δS + δS², with δS the largest |ΔS| measured at full precision for the same code
  pair (item 18b), or, before item 18b exists, the predicted δS = 1e-11.
- Values that are reported as a count of "compared / within" always include the worst ratio
  |Δ| / tolerance and where it occurs.

### BK coefficients (all analysis steps, all channels, both vibrational levels)
- Global: max |Δc| ≤ 1e-10 · max(max|c|, 1e-4) (gate 4).
- Per channel: relative error |Δc| / |c| reported for |c| above the floor 1e-6 · max|c|
  (a small channel could be wrong and still pass the global criterion). No pass/fail threshold
  is set for it in advance; it is checked against the prediction below.
- Time resolution: the error growth law of max |Δc(t)| / max|c| is fitted (√t random rounding,
  ∝ t systematic bias, a jump = a logic difference, e.g. a different analysis-line step).

### Potential (identical physical points)
- |ΔV| ≤ 1e-13 · max(|V|, 1 eV) (the `.pes` is in hartree: floor 1 eV = 0.036749… hartree via
  `legacy` 27.2114). Also reported separately for the region V < 2 eV.
- Different symmetry representatives (C++: first point of an orbit in memory order; Fortran:
  triangle 1) are reported separately, bounded by the PES's own asymmetry computed over the
  whole grid (max |V(q) − V(g q)|, C++ PES without symmetry).
- V_asy(r): ≤ 1e-14 · max(|V|, 1 eV) (gate 1 criterion; both sides at full precision).

### Full-precision C++ data (CPU↔GPU tables, replay, `onlyAnalysis` pairs)
- P, S (complex, |ΔS|) and P_reaction: ≤ 1e-10 absolute (never relative for P_reaction).
- Phases only where |S|² > 1e-8, difference wrapped into (−π, π]; Fortran quirk: phase 0 when
  Re S = 0 (reproduced by the C++, so compared as printed).
- Averages (norm, ⟨Z⟩, ⟨r⟩, rest norm) CPU↔GPU: ≤ 1e-10 relative (rest norm absolute 1e-10).
- ψ: ‖Δψ‖₂ / ‖ψ_final‖₂ (no global-phase alignment: a phase drift is diagnostic). Reported, no
  threshold beyond the prediction.

### Thesis (item 17), defined before any run is drawn on the digitised curve
- Digitisation: the dotted curve is extracted from the PDF's vector drawing (pymupdf
  `get_drawings`), axes calibrated on the tick marks. u_dig = the largest calibration residual
  of the ticks (in P and E) plus half the drawn line width, both converted to data units.
- Compared only on 0.40–0.95 eV (the runs' energy range), after confirming that the thesis
  x-axis is the collision (normal incidence translational) energy E_coll used by the runs.
- **Agrees**: |P_run(E) − P_thesis(E)| ≤ u_dig + 0.01 at every digitised point in range (0.01 =
  allowance for the different set-up: ΔZ 0.15 vs 0.12, N_Z 128 vs 160, T 25 000 vs 17 500 a.u.,
  unknown packet and absorbers; no convergence study is done here). Reported: max |ΔP|, where,
  and the fraction of points within. A difference with the thesis is not a port error when the
  three codes coincide; it would be marked "open" (set-up), not FINDING.

### Every comparison reports
n, median, p95, max, and where the max occurs (E, channel, t, grid point).

## Prediction (written before comparing)

Known sources of C++↔Fortran differences in this run, with the sizes documented for the small
gate grids (context, not expectations for this grid):

| Source | Where | Documented size |
|---|---|---|
| quad-precision sites (`-fdefault-real-8` promotes `DFLOAT`/`d0`) | grids, momenta, kinetic energies, CAPs, analysis energies | ~1 ulp per value; ≈ 1e-12 phase after 7000 steps |
| libmvec vs scalar libm | θ grid acos, Legendre exp, cos/sin of exp(−iV dt) | a few ulp per value |
| oneMKL vs FFTW (FFT algorithms) | every FFT | ≤ 1e-14 per stage; CPU vs GPU ≈ 1e-13 after 1000 steps |
| fast-math (host) / FMA (GPU kernels) | both C++ runs | norm 2e-12 relative, P 1.5e-12 absolute after 7000 steps (cu111_dailight_small) |
| LAPACK (C++) vs Jacobi (Fortran) | gas-phase states | gate 2: energies ≤ 1e-12 relative, ψ ≤ 1e-10 (loose bound; measured far smaller) |
| BK C++ vs Fortran, strict build | cu111_dailight_small, 7000 steps | 6.4e-13 of max\|c\| |
| P_reaction C++ vs Fortran | cu111_dailight_small | 3.6e-13 / 4.1e-12 strict, 2.0e-12 / 5.3e-12 fast-math |

This grid has the same number of steps (7000) but JMax 28 (435 rotational functions instead of
15), 12 × 12 (x, y) points and longer sums in the angular transform, so per-step rounding is
somewhat larger (sums of n terms: ~√n ulp).

Predicted for this run:
1. Set-up integers, decisions, channel lists: identical. Printed set-up reals: within the
   printed precision. Grids: ≤ 2 ulp (quad-precision sites). θ nodes ≤ 1e-15 absolute.
2. Potential at identical points: ≤ 1e-13 · max(|V|, 1 eV) (gate: ≤ 3.6e-14); V_asy ≤ 1e-14.
   CPU↔GPU potential: identical (same host code and flags; V is computed on the host).
3. Propagation, C++↔Fortran at full precision: BK max |Δc| / max|c| ≈ 1e-12, at most 1e-11;
   growth ∝ √t to ∝ t, no jump. |ΔS| ≲ 1e-11, |ΔP| ≲ 1e-11 absolute, |ΔP_reaction| ≲ 1e-11.
4. CPU↔GPU (both fast-math, different rounding: host fast-math vs FMA kernels, oneMKL vs cuFFT):
   the same order as item 3, ≈ 1e-12 (the in-run rounding floor). Claim to test: C++↔Fortran
   is of the order of CPU↔GPU, i.e. within a factor ~10.
5. Printed comparisons: everything within ½ printed unit, except printed values whose half unit
   is below the absolute rounding of item 3 (probabilities below ~1e-8 printed with 7 digits,
   P_reaction near rounding boundaries): those are expected to exceed the printed criterion by
   at most |ΔP| ≤ 2 √P · δS + δS² (δS ≈ 1e-11).
6. All three are far below the print resolution (≥ 5e-8 · P) and the visible artefacts of the
   calculation (norm left at T ≈ 1e-4, a sub-threshold reaction floor ~1e-4, asymmetry of
   equivalent channels): (b) ≈ (a) ≪ (c) ≪ (d).
7. Analysis line: fixed at the same analysis step in all three runs (C++ log: t = 2220 a.u.;
   the Fortran step is to be pinned through fort.21).

A C++↔Fortran difference above 1e-10 (absolute in P or S, relative to max|c| in BK) would exceed
this prediction by more than an order of magnitude and would be investigated as a possible
logic difference; above the fast-math/FMA level it may motivate a strict-build run (D2).

## Criteria changes
- 2026-11-07, while writing `setup_items.py` (before looking at its results): three set-up cases
  the list above did not cover.
  1. Set-up scalars that both programs hold in double precision (fort.21 header values, `.pes`
     grid records, C++ `%.16e` log values): **identical**, or **rounding** if
     |Δ| ≤ 1e-13 · |x| (the potential criterion; the predicted size is a few ulp).
  2. Rovibrational energies at full precision: ≤ 1e-12 relative (the gate 2 criterion,
     `tests/test_regression_pes.cpp`, LAPACK vs Jacobi).
  3. Where the C++ side is itself printed with fewer than 16 digits (e.g. `%.6f` in the log), the
     tolerance is the sum of both half units (each 1.01 × ½ unit of its last digit).
- 2026-11-07, after the first run of `item14_16_analysis.py` (so after seeing results; this changes
  how an exceedance is *explained*, not any pass/fail threshold). The predicted bound
  2 √P · δS + δS² did not cover P_reaction (it is 1 − ΣP) nor single S components, and it needs
  a δS. Replaced by a measured decomposition through the replay (item 18b): a printed value
  outside its half unit is **explained** when (a) replay(fort.21), i.e. the Fortran's own BK
  coefficients through the C++ analysis, reproduces the printed Fortran value within its half
  unit, and (b) |C++ − replay(fort.21)| ≤ 1e-10 (the full-precision criterion). Where (a) fails
  (20 S components, the analysis side), the first-order bound of the analysis difference is
  computed per entry from known deviations: |factor| dt Σ_t |c_t| (δE t + ulp(|E| t) + 4ε) +
  |S| (Z_ana δk + δk/(2k)), with δE the eviroi difference (LAPACK vs Jacobi, item 6) and δk the
  channel-momentum difference from the channel energies (item 9); the value is explained when
  |replay − F| ≤ half unit + bound.
- 2026-11-07, after review round 2: the pre-registered "fitted growth law" of ladder 2 (√t / ∝ t /
  jump) is not identifiable from these data: the signal window spans a factor 2.2 in t, an
  unweighted linear fit of the amplitude part gives R² 0.66, the |c|²-weighted fit a different
  slope. Reported instead: the window quartile medians of the amplitude part, both fits, and the
  absence of a jump at the fixing step (a discrete decision would show as one).
- 2026-11-08, strict-build GPU run (the D2 proposal; `strict_gpu.py`). Written while the run was
  still propagating, before any of its results were read. Same input, `build-cuda`
  (`SPODVR_FAST_MATH` OFF: no `-ffast-math` on the host, no FMA in the kernels), into
  `runs/thesis_fig5_3_dz_cuda/`. The question is where the uniform amplitude part κ and the
  specular rotationally inelastic excess come from (README, composition of the difference).
  Estimators are those of items 13c and 18b, with the strict run in the GPU's place; replays use
  the same `build-fast` helper as before, so replay vs replay differences are propagation only.
  1. Validity: the run finished (exit status 0, 7000 steps), the binary is built from the same
     source as the fast runs, the analysis line is fixed at the same step (t = 2220 a.u.).
  2. κ, the median over energies of the uniform amplitude part of replay(strict) − replay(fort.21)
     (fast runs: CPU 6.0e-12, GPU 6.6e-12). **Fast math is the cause** if |κ| ≤ 1.5e-12 (a drop by
     a factor 4, to within 3× the GPU−CPU spread of the fast runs); **not fast math** if
     |κ| ≥ 4.5e-12 (≥ 75% of the CPU value); in between: partly. The BK-level amplitude part
     (window median, fast CPU 6.7e-12) is reported as the second estimator.
  3. Specular rotationally inelastic family, the j = 10 norm-weighted residual after the uniform
     factor (fast runs 8.4e-9): **fast math involved** if ≤ 8.4e-10 (a factor 10); **not fast
     math** if ≥ 4.2e-9 (≥ 50%); in between: partly.
  4. Strict vs fast GPU, i.e. the size of the fast-math effect: predicted ≈ 1e-12 relative (the
     documented fast-math level), as the norm-weighted relative S difference.
  5. End to end, the criteria of items 13a and 14: averages and P_reaction within half a printed
     unit of the Fortran (P_reaction floor 1e-11); |ΔP_reaction| > 1e-6 would be investigated.
  6. Replay self-check with a helper linked against `build-cuda`: predicted **bitwise** (without
     fast math the rounding does not depend on inlining). If not, the offset is reported.
- 2026-11-08, after the first run of `strict_gpu.py` (after seeing its results): the validity check
  "same input" first required the INPUT echo of the strict run's log to equal the fast GPU run's
  verbatim. The echo also prints derived values at full precision (E_trans, α, k_Z, θ nodes and
  weights, the PES symmetry self-check), which differ by fast-math rounding (largest 1.3e-13, a θ
  quadrature weight). The input identity now rests on `input.json` (unchanged since before the fast
  runs, clean working tree, named by both commands) and on an echo whose text agrees apart from
  numbers; the numeric differences are reported. The replay of the strict BK file with the
  `build-fast` helper goes through a header copy (`strict_bk_copy.py`): the reader's exact check of
  the stored restart parameters rejects the strict file over one 1-ulp parameter.
