# Review log

Reviewers are fresh subagents without edit, build or run rights (grep/sed/short Python on the raw
files allowed). Each round gets the document, the raw run paths, docs/validation.md and
docs/fortran_deviations.md, not the scripts. Triage: accept / reject, with the reason; numeric
claims of reviewers are checked before acting; a rejected high-severity point needs logged data.

## Round 1 (2026-11-07): professor (P), PhD student (S), adversarial verifier (V)

Reviewers' numeric claims were re-derived before acting (`item13_propagation.py`
"specular_inelastic_family", "amplitude_vs_time"; `item18b_replay.py` "S_relative_normweighted",
residual breakdown; a per-operation symmetry check). The verifier's four headline re-derivations
(1.884e-12; 109/111 per run; 6.07e-12 / 7.30e-12 / 1.24e-12; step 111 in all three files) match
the document.

| # | From | Sev. | Point | Decision | Reason / action |
|---|---|---|---|---|---|
| 1 | V | high | The C++↔F difference is not a uniform amplitude scale: after removing α, ~60% of ‖ΔS‖ remains; the specular rotationally inelastic channels (v=0, (0,0), mJ=0, J=2–10) differ far more than GPU↔CPU | accept | Verified: residual 57% of ‖ΔS‖ (4.2e-12 relative, spread over all channels); in BK time domain 98.5% of ‖Δc‖² before the fixing step and 2.2% over the run sit in that family, \|Δc\| 5e-14…3.3e-13 (1.7e-11 rel. at J=2, 7.9e-10 at J=10), GPU↔CPU ~8e-13 rel. Document rewritten: three components (uniform 6e-12, diffuse 4e-12, channel-local ~1e-13), all within criteria, none attributed: **open**. Hypothesis stated as such (angular round trip of the specular component; θ nodes differ ≤ 6 ulp, Fortran libmvec Legendre). |
| 2 | P, S, V | high | Prediction (≲1e-11, within ~10× of CPU↔GPU) was missed (1.4e-11, 1.6e-11; 13×); "upper end of the prediction" understates it | accept | Wording: "missed"; a prediction-vs-outcome table in the README. Strict-build GPU run proposed to the user (D2). |
| 3 | P, S, V | high/med | GPU↔CPU is not "the rounding floor": they share all tables, V and the operator decomposition | accept | Rung 1 renamed "two executions of the same C++ code"; scope says it cannot show table-level deviations. |
| 4 | S | high | Nothing shows the comparisons could detect a failure; negative control (deliberate perturbation) | partly | A perturbed-table run is a new run (D6: not allowed) → rejected. Added "What these comparisons resolve, and what they would miss" (resolution per comparison; deviations of a few ulp that were detected and traced: 3/11-ulp k_Z/α, 1.2e-15 hartree eviroi; the analysis-line step; the channel-local excess of point 1) and the verifier's smallest-undetected-error analysis. |
| 5 | P, S | med/high | "per-step norm bias ~8e-16" is α/7875 restated and mislabelled; t^0.69 fitted over a factor 2.2 in t cannot tell √t from t; GPU↔CPU t^1.76 unexplained | accept | Replaced by the measured amplitude part α(t) of the BK difference per step: 6.3e-12 → 9.8e-12 across the signal window, slope 4.9e-16 per a.u. (≈1.1e-15 per kinetic block in amplitude, 2.2e-15 in norm), intercept 2.2e-12. Growth-law exponents removed from the text and Figure 3. |
| 6 | P, V | med | "absolute floor that GPU↔CPU lacks … shared BK projection" is not what the data show | accept | Verified: before t = 2220 C++↔F is 62× GPU↔CPU (median), carried by the specular inelastic family; at 2240–5000 a.u. both pairs are ~2–3e-14. Attribution to the BK projection removed. |
| 7 | V | med | "explained" rule written after seeing results; its condition (b) (≤1e-10) can hardly fail | accept | README says so with the counts (2 + 17 + 56 + 46 + 502 + 380 values, CPU) and that (b) is weak; the decisive part is (a), the analysis-only reproduction. |
| 8 | V, S | med | Verdict's "1.4e-11 C++ vs Fortran at full precision" is propagation only (C++ analysis of both BK files) | accept | Verdict and Figure 1 caption say so; the analysis is checked separately at printed precision. |
| 9 | V | low | "2α·ΣP reproduces ΔΣP" is a first-order identity for the least-squares α | accept | Removed. |
| 10 | V | med | "180/180" averages counts the time column | accept | 144 values (36 times × 4 columns). |
| 11 | V | low | Item 19 mirror: Fortran "2.3e-16" is from printed values; in the BK data both codes break the mirror at ~1e-13 | accept | Row labelled "printed values"; verdict wording "the same in all codes". |
| 12 | P | low | Item 19 verdict "identical" while precisions differ | accept | "same in all codes, within each code's precision". |
| 13 | P | med | Symmetry breaking only "consistent with" the grid; every operation except (n,m)↔(m,n) is broken by 3–5%; orbits inside the grid break as much as those crossing it → dynamics (T_xy on the parallelogram k-grid), not labelling; show with nX = nY = 18 | partly | Verified per operation: median 2.5–3.9% (p95 46–59%), the diagonal mirror 2.6e-14. Wording changed to the dynamics; the larger-grid run is a new run (D6) → rejected, listed as a follow-up. Added: a bug common to both codes cannot be seen by a port comparison. |
| 14 | P | med | Sub-threshold floor not quantified (replay truncated BK files) | reject | Physics of the input, identical in all codes; not needed for the port verdict. Listed as a follow-up. |
| 15 | P, V | low/med | Thesis: 0.01 allowance undeveloped, tolerates ±8 meV; the −4.6 meV offset is real; item 17 tests the reconstructed input, not the port | accept | Item 17 reframed: input reconstruction at the 1e-2 level, not a port test; a GPU run on the thesis grid is a follow-up (D6). |
| 16 | P | — | Explain the amplitude excess (α vs t̄(E), per-operator unitarity budget, strict GPU run) | partly | α(t) measured (point 5); FFT scales quantified (≈6% of the per-block rate); unitarity budget of the Fortran operators needs Fortran-side tables we cannot run (D6) → follow-up; strict GPU run proposed (D2). |
| 17 | S | med | Undefined terms (BK, c, replay, specular grid, analysis line, rung, signal window, gate 4, Phase 0, D2, CAP, eviroi, kinetic block), DZ vs Dai–Light, thesis not cited, footnote marks read as squares, α overloaded, figure order 1-4-3-2 | accept | Terms box in the README; DZ/Dai–Light explained; thesis cited; markers (p)/(f); amplitude excess called κ; figures renumbered in reading order. |
| 18 | S | med | First screen lacks the method; Fig. 1 crosses look like failures; GPU printed mismatches not plotted; speed in the correctness verdict | accept | One method sentence; crosses for both runs labelled "outside ½ unit, explained"; speed moved below the verdict. |
| 19 | V | low | "within 1.3e-11 of a rounding boundary": the GPU case is 1.34e-11 | accept | "within 1.4e-11". |
| 20 | P | high | Headline should say both binaries are the experimental fast-math builds that fail one project gate | accept | Verdict names the builds; the gate failure stays in "What ran" and limitations. |

## Round 2 (2026-11-07): professor (P), PhD student (S), adversarial verifier (V), on the revised document

The verifier re-derived every new number of round 1 (α(t) fit, 98.52% / 2.22%, per-j table, κ,
57.2% residual, per-operation breaking, 144 averages, 218/222, explained counts) and found them
correct, with the corrections below. New analyses added in response: |c|²-weighted and
through-origin fits with R² and window quartiles (`item13_propagation.py`), the amplitude part of
ψ(T) on both grids, a global-bias test of the analysis and κ for GPU↔CPU (`item18b_replay.py`),
the per-energy print margin (`fig2_ladder.py`), and a negative control (`negative_control.py`).

| # | From | Sev. | Point | Decision | Reason / action |
|---|---|---|---|---|---|
| 1 | P, V | high/med | The per-block bias (2.2e-15) and "FFT scales 6%" rest on an unweighted fit driven by the weak tail; a constant offset fits as well | accept | Verified: unweighted slope 4.9e-16 (R² 0.66), \|c\|²-weighted 2.1e-16 with intercept 4.2e-12, window quartile medians 6.5, 6.2, 6.9, 9.5e-12. Text now: an offset of ~6–7e-12 established early, plus a rise in the weak tail; the per-step reading is model-dependent and dropped; the FFT scales predict 4.9e-13 by T (< 10% of the offset). Amplitude part of ψ(T): specular grid 7.7e-12 (of 7.8e-12 total), scattering grid 1.4e-11. |
| 2 | S, P | high | Detection power asserted, not shown; an analysis-side negative control is in scope | accept | `negative_control.py`: a 1e-10 relative error injected into one channel of a fort.21 copy stands 4.7× above the all-channel p95 of the per-channel residual at all 111 energies (P_reaction moves by only 7.2e-12); 1e-11 rises 8× above that channel's own baseline but not above p95. Global analysis bias from the print-rounding residuals: \|S\|² −1.4e-9 ± 1.3e-9, phase 6.6e-10 ± 5.9e-10 rad, residual std 0.577 half units (pure rounding: 0.577). |
| 3 | S | high | The per-channel residual did flag a family (22–766× GPU↔CPU), waived as "open"; put its size in the verdict or test the θ-node hypothesis | partly | Its relative size is now in the verdict and the composition table (j = 2…10, 1.0e-11…8.4e-9 relative in S). The θ-node test needs the Fortran's FBR⇄DVR matrix (Legendre values with libmvec) or a code change: not done, listed as follow-up; the item stays **open**. |
| 4 | S | high | The verdict covers only the experimental fast-math binaries, one of which fails a gate | accept | Verdict names this; strict-build GPU run proposed to the user (D2). |
| 5 | S | med | Prediction table grades "yes" where the rule came afterwards / the central estimate was missed | accept | "partly" and "within the bound, above the central estimate". |
| 6 | S | med | Item 13d verdict cell changed metric (pre-registered ‖Δψ‖/‖ψ_T‖) | accept | Both numbers shown. |
| 7 | S | med | "5×10³ above rung 2" compares medians; at 0.40–0.53 eV the margin is ~4× | accept | Minimum margin 3.7× (0.40 eV) stated. |
| 8 | S | low | Ratio 12× vs ×13; per-channel 9.5e-13 vs P(v,j,mJ) 1.1e-11 | accept | Metrics named; P(v, j, mJ) sums channels. |
| 9 | P | med | "Primary evidence \|ΔS\| ≤ 1.7e-12" is small only because no \|S\|² > 0.079; quote relative | accept | Large channels (\|S\|² > 1e-2): \|ΔS\|/\|S\| 7.0e-12 (median), ≤ 7.7e-12. |
| 10 | P, V | med/low | Mixed units in the three-part composition; j = 10 missing | accept | All relative (S domain); j = 10 included. |
| 11 | P | med | Limitation (ii) understates the blind spot: anything below the 1e-10 criterion passes | accept | Rewritten with the negative-control numbers. |
| 12 | V | med | Smallest analysis error missed: per S element up to ~5e-7 (fort.25 has 6 digits), not 5e-8 | accept | Limitation (i) corrected; the global bias test bounds a channel-wide analysis bias to ~1e-9. |
| 13 | V | med | GPU↔CPU also has a uniform κ (median 5.3e-13, slope 2.0e-16) — not a C++↔F signature; "stays near 5.5e-13" is wrong | accept | GPU↔CPU column in the composition table; text corrected. |
| 14 | V | low | Item 16: 20 analysis-side components are in the replay; only 14 of the 882 explained C++ values fail condition (a) | accept | Corrected. |
| 15 | V, P | low | ×17 at the fixing is on ‖Δc‖, Figure 3 plots the max metric (×165); worst averages ratio is relative to the tolerance | accept | Both stated. |
| 16 | P | low | Item 4: "GPU potential bitwise" is inferred from the relinked helper | accept | "the helper's V from both builds". |
| 17 | P | low | The pre-registered growth-law fit was dropped without a dated entry | accept | criteria.md entry and a prediction row. |
| 18 | S | med | First screen: verdict reads hedged (15 numbers); Fig. 1 crosses look like near-failures, 4 crosses but 2 visible; terms (quad-precision sites, libmvec, Jacobi, onlyAnalysis, fort.25, …), α/β overloaded, j vs J | accept | Verdict = claim + three evidence bullets + one open item; caption explains the crosses (print rounding included, CPU and GPU overlap); Terms box extended; κ(t) for the amplitude part; j throughout. |

## Round 3 (2026-11-07, final; cap reached): professor (P), PhD student (S), adversarial verifier (V)

The verifier re-derived the round-2 numbers (composition table both columns, κ(t) quarters, fits,
specular ψ(T) amplitude, global bias, 218/222, 20 / 14 of 882, 98.5% / 2.2%, step ratios, gate
margin) and found them correct. The professor would accept the conclusion "with a minor revision".
New script outputs: per-element statistics of the specular inelastic family, the per-channel
residual by channel size, the de-duplicated bias, κ(E) vs 1/√E, the GPU print margin
(`item18b_replay.py`, `item14_16_analysis.py`, `fig2_ladder.py`).

| # | From | Sev. | Point | Decision | Reason / action |
|---|---|---|---|---|---|
| 1 | P, S, V | high | Headline "~1e-11" holds for P (absolute) and populated channels; per channel the specular inelastic family reaches 8.9e-10 (j = 6, 8, P > 1e-8) and per element up to 1.6e-6 relative (j = 10); state it in the conclusion, with units | accept | Verified (j = 8: median 3.6e-10, max 8.0e-7; j = 10: median 9.7e-8, max 1.6e-6). Verdict rewritten: propagation (absolute in P; relative per populated channel), analysis (print resolution), and the family's size in the conclusion sentence. |
| 2 | V | high | The pre-registered "within ~10× of CPU↔GPU" fails per channel by 31–303× (S) / 22–766× (BK) | accept | Prediction-table row added, marked "no". |
| 3 | S | high | The family is graded "rounding, within the criteria" only because no relative per-channel threshold was set; a FINDING was impossible | accept | Verdict rows 13c, 16: "rounding (criterion) · open: channel-local excess"; the text says the absolute criterion could not flag it and no relative threshold was pre-registered. |
| 4 | P | med | Selection-rule argument: the specular wave is (x, y)-uniform, φ-independent, j = 0, so any angular round-trip defect lands exactly in v=0, (0,0), m_J=0, even j, with \|Δc\| independent of \|c_j\|; testable offline with the θ nodes of both codes | partly | Argument added (it matches the observation); "not yet attributed; signature consistent with the angular transform acting on the specular component, testable offline". The offline test needs the Fortran's own weights and Legendre values (libmvec), not in its output; listed as follow-up, not done. |
| 5 | P | med | "GPU↔CPU shows a uniform part too, so this is the form rounding takes" runs backwards: κ(GPU−F) 6.6e-12 ≈ κ(CPU−F) 6.0e-12 + κ(GPU−CPU) 0.5e-12, so the common 6e-12 lies in what both C++ executions share (host tables, fast-math host code, operator decomposition) or on the Fortran side | accept | Rewritten. This also motivates the strict-build run (it changes the host tables). |
| 6 | V | med | Negative control follows from linearity; its "all-channel p95" covers only \|S\|² > 1e-8; the detection threshold depends on channel size | accept | Rewritten as a sensitivity statement by channel size: natural residual max 3.2e-11 (\|S\|² 1e-4…1e-2), 2.1e-10 (1e-6…1e-4), 8.9e-10 (1e-8…1e-6). |
| 7 | V | med | 10 of the 46 explained \|S\|² values are the family with real propagation differences ~1e-6 relative, not "small values or rounding boundaries" | accept | Stated. |
| 8 | V | low/med | Bias uncertainty ignores symmetry duplicates (−1.5e-9 ± 1.9e-9 de-duplicated): "≲ 5e-9", not "~1e-9" | accept | Both numbers; "≲ 5e-9". |
| 9 | V, P | low | The 20 analysis-side components: 3 miss by 2.0–2.5 half units, all are small components (≤ 3.7e-3 \|S\|), the bound is 10–1000× loose; the set-up-correction test (20 → 20) was not reported | accept | Wording: a loose bound, a weak explanation; the null result of the correction test reported. |
| 10 | S | med | Prediction table omits the missed "printed set-up reals" (\|ψ(k)\| tails, ~340 half units) | accept | Row added ("partly": 2 values, explained). |
| 11 | S | med | Rung 3 (print resolution) is that of P_reaction only; the state-resolved values have ~1000 printed values outside ½ unit | accept | Stated in the ladder. |
| 12 | S | low | Worst-ratio locations (criteria.md requires where) missing in items 14–16 tables | accept | "at" column added. |
| 13 | S | med | Terms: fast-math vs strict, FMA, DVR⇄FBR, rest norm, iLow/iHigh, enymo, "Phase 0", α in a table header; Terms box comes after Figure 1; why no strict run | accept | Terms box moved above the verdict and extended; header κ; "Phase 0" replaced; the strict run is the owner's decision (D2), said so. |
| 14 | V | low | GPU print margin 3.2× (CPU 3.7×); 57% is a norm share (33% of ‖ΔS‖²); 62× vs 58× pre-fix ratio depends on the subset | accept | Both margins; "57% of the norm (33% of ‖ΔS‖²)"; ratio described as subset-dependent (median over the pre-fix records). |
| 15 | V | — | Smallest error still missed: a uniform amplitude error ~5e-12 on the incoming path (κ(E) correlates with 1/√E, r = 0.96); channel-local errors up to ~1e-9 at \|S\|² ~ 1e-8; up to ~1e-6 in the family | accept | Added to "What these comparisons would miss". |
| 16 | S | low | BK agreement listed under "Set-up"; thesis 67/67 in the performance line | accept | Moved. |
