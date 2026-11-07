# Status

One line per item and figure: status, script, verdict. Tiers from the task (A must, B, C).

| Item | Tier | Status | Script | Verdict |
|---|---|---|---|---|
| Phase 0: runs finished | A | done | (shell) | all three exit 0, final files present |
| Phase 0: provenance | A | done | `provenance.py` | builds at abc47f0 (src unchanged to HEAD e0c231b); binaries older than runs |
| Phase 0: regression gates | A | done | (ctest) | build-fast 20/21 run pass, build-cuda-fast 23/24 run pass; both fail RegressionHardWall.Gate1AsymptoticPotentialTable by 1.64e-14 vs 1e-14 (fast-math; strict `build` passes it): explained |
| 1 input parameters | A | done | `setup_items.py` | identical/printed; input.json ≡ DVR.inp except Title (by design); echo checks all pass; n/a: r cut-off, mixing, MjIn message |
| 2 grids | A | done | `setup_items.py`, `item04_potential.py` | identical (grids, φ, counts); θ nodes ≤ 6 ulp (rounding, libmvec acos) |
| 3 MISC SETUP lines | A | done | `setup_items.py` | identical/printed; nmoder n/a |
| 4 potential (full grid) | A | done | `potential_check.cpp`, `item04_potential.py` | rounding: same point max 2.8e-14, used 7.5e-14 (≤ 1e-13); CPU↔GPU identical (shared code); helper self-check bitwise |
| 5 V_asy(r) | A | done | `setup_items.py` | rounding: max 5.1e-16 (≤ 1e-14), 28/40 identical |
| 6 rovibrational states | A | done | `setup_items.py` | printed: 348 values, worst 0.99 half unit; no sign flips; E(0,0) full 7e-15 rel. (rounding) |
| 7 diffraction basis | A | done | `setup_items.py` | identical (incl. warning, maxnum derived) |
| 8 rotation channels | A | done | `setup_items.py` | identical |
| 9 rotation-diffraction channels | A | done | `setup_items.py` | identical labels of all 9504 states; channel energies rounding (max 1.3e-14 rel.) |
| 10 optical potentials | A | done | `setup_items.py` | printed |
| 11 initial wave packet | A | done | `setup_items.py`, `item11_psi_k.py` | printed/identical/rounding; \|ψ(k)\| tails explained (reproduced: Fortran α, k_Z quad-precision site + 1e-10 cancellation) |
| 12 analysis set-up | A | done | `setup_items.py` | identical/printed; energy grid 111/111 within half unit |
| 13a averages | A | done | `item13_propagation.py` | printed: 144/144 (CPU and GPU; 36 times × 4 columns), worst 0.95 half unit; GPU↔CPU norm 5.0e-12 rel., ⟨Z⟩ 5.6e-14, ⟨r⟩ 1.9e-14, rest norm 1.1e-13 abs (rounding) |
| 13b analysis line fixing | A | done | `item13_propagation.py` | identical: Z = 15, record 111 (t = 2220 a.u.) in all three; jump 1e-19 → 4.4e-4 pins the Fortran step |
| 13c BK coefficients | A | done | `item13_propagation.py` | rounding: max \|Δc\|/max\|c\| 6.1e-12 (C++−F), 7.3e-12 (GPU−F), 1.2e-12 (GPU−C++) ≤ 1e-10; amplitude part α(t) 6.3e-12 → 9.8e-12 in the window; channel-local excess in the specular inelastic channels (98.5% of ‖Δc‖² before the fix): **open** (not attributed) |
| 13d final ψ | B | done | `item13_propagation.py --psi` | rounding: scattering grid rel. 2.0e-10 (C++−F), 3.4e-11 (GPU−C++); relative to the initial norm 1.2e-13 / 1.9e-14; specular 8e-12 rel. |
| 14 probabilities | A | done | `item14_16_analysis.py` | printed + explained (replay): P_reaction 109+2/111, P(v,j,mJ) 14635+17/14652, ΣP, P(v), ⟨E⟩ printed; full precision C++−F ≤ 1.4e-11 (uniform κ 6e-12 + diffuse 4e-12 + channel-local): rounding by criterion, prediction narrowly missed, composition **open**; GPU↔CPU ≤ 1.9e-12 |
| 15 diffraction probabilities | A | done | `item14_16_analysis.py` | printed + explained: 82432+56/82488 (CPU), 0 outside; GPU↔CPU 6.4e-13 |
| 16 S-matrices | A | done | `item14_16_analysis.py`, `item18b_replay.py` | printed + explained: \|S\|² 636327+46/636373, Re/Im S 0 outside, phases 327734/327734; full precision \|ΔS\| ≤ 1.7e-12; specular inelastic family up to 1.6e-6 relative per element: **open**; GPU↔CPU \|ΔS\| 5.5e-13 |
| 17 thesis Fig. 5.3 | A | done | `thesis_extract.py`, `item17_thesis.py` | agrees (criterion: u_P + slope·u_E + 0.01): 67/67 points in 0.40–0.95 eV for all three; max \|ΔP\| 0.014 at 0.632 eV; runs +0.007 on average (−4.6 meV shift: rms 0.003) |
| 18a onlyAnalysis cross-check | A | done | `item18a_only_analysis.py` | rounding: all pairs ≤ 1e-10; C++−F P_reaction median 1.03e-11 (max 1.36e-11), S ≤ 1.7e-12; GPU−C++ 7.9e-13 median: same picture as the replay |
| 18b replay | B | done | `replay.cpp`, `item18b_replay.py`, `negative_control.py` | analysis-only: every printed probability and phase reproduced, all but ≤ 20 S elements (loose bound); global bias ≲ 5e-9; propagation-only C++−F \|ΔS\| ≤ 1.7e-12, P_reaction ≤ 1.4e-11; composition κ 6e-12 (shared by both C++ executions) + diffuse 4e-12 + specular inelastic family: **open**; self-check not bitwise (≤ 2.0e-14, fast-math LTO) |
| 19 physics checks | B | done | `item19_physics.py` | same in all codes: open-channel sets equal (margin 1.3e-4 eV); 0 ≤ P, ΣP ≤ 0.99975; P(mJ) = P(−mJ) to 4.7e-12; (n,m)↔(m,n) mirror exact to rounding, all other C6v operations broken by median 2.5–3.9% identically (dynamics on the parallelogram (x, y) momentum grid); sub-threshold P_r 2.5–5.2e-4; norm left 1.2e-4 |
| Figure 1 reaction probability | A | done | `fig1_reaction.py` | rendered and inspected (167 KB) |
| Figure 4 execution time (was 2) | A | done | `fig4_time.py` | rendered and inspected; table `tables/timing.md` |
| Figure 3 BK differences in time | B | done | `fig3_bk_growth.py` | rendered and inspected; bottom panel: amplitude part κ(t) and its window-quarter medians (after review round 2) |
| Figure 2 error ladder (was 4) | B | done | `fig2_ladder.py` | rendered and inspected |
| README.md / stages.md / appendix.md | A | done | `tables.py`, `fill_tables.py` | revised after three review rounds; verdict with units; open items stated |
| logs in data/ | A | done | `make_all.sh --extract` | raw `.log.gz` (gzip -9 -n), cleaned logs |
| review round 1 | A | done | `review_log.md` | 20 points: 16 accepted, 3 partly, 1 rejected; new finding (verifier): channel-local excess, document rewritten |
| review round 2 | A | done | `review_log.md` | 18 points: 17 accepted, 1 partly, 0 rejected; added negative control, global analysis-bias test, ψ(T) amplitude, fit variants; per-step reading dropped |
| adversarial verifier | B | done (rounds 1, 2) | `review_log.md` | headline numbers re-derived and matching |
| review round 3 | B | done (cap) | `review_log.md` | 16 points: 15 accepted, 1 partly; professor: accept with minor revision; open follow-ups: offline angular-matrix test, strict-build GPU run (D2) |
| number check | B | done | `check_numbers.py` | 0 unmatched after round 3; weak check (coincidental matches possible), headline numbers also checked by hand |
| Markdown check | A | done | `check_markdown.py` | 0 problems (links, anchors, tables, details, ignored files, sizes, whitespace) |
| make_all.sh --heavy | C | cut (tier C) | `make_all.sh` | default and `--extract` modes run clean; every `--heavy` step ran individually during the work, the script as a whole was not rerun |
| sha256 manifest | C | done | `manifest.py` | 66 files |
| mermaid diagram | C | cut (tier C) | — | not drawn: time spent on the review findings instead |
| byte-identical PNG reruns | C | done | `make_all.sh --extract` | all 5 PNGs byte-identical after a full extract rerun |
