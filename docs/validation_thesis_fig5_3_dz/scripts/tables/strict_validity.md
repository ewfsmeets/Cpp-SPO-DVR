| Check | Result |
|---|---|
| Run | exit status 0, 7000 steps, wall 1.43 h (1.12× the fast GPU run) |
| Build | `build-cuda/spodvr` (2026-11-06 15:57), FAST_MATH OFF, ARCH x86-64-v3; source commits after the build: none |
| Input | `inputs/thesis_fig5_3_dz/input.json` unchanged since before the fast runs and clean: yes |
| Input echo vs the fast GPU run | 84 lines: 65 identical, 17 with derived values differing (≤ 1.3e-13 rel., 1 above 1e-13), 2 PES self-check diagnostics, 0 text differences |
| BK file header vs the fast GPU run | 34/35 restart parameters identical; ZOptFactor (effective): 1.4e-16 rel. |
| Analysis line | fixed at Z = 15 a₀, t = 2220 a.u. |
| V_asy and gas-phase states vs the fast GPU run | max \|Δ\| / max \|value\|: 6.5e-17 and 2.3e-14 |
| Averages vs Fortran (item 13a rule) | 144/144 within ½ printed unit, worst 0.95 of it |
| P_reaction vs Fortran printed (item 14 rule) | 111/111 within ½ unit (floor 1e-11), 0 explained, 0 outside |
| Replay self-check, helper linked against `build-cuda` | bitwise |
