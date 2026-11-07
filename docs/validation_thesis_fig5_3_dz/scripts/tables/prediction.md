| Pre-registered prediction | Outcome | Met? |
|---|---|---|
| set-up integers, decisions, channel lists identical | identical (items 1–3, 7–9, 12) | yes |
| potential at identical points ≤ 1e-13 | 2.8e-14 | yes |
| C++↔Fortran BK ≈ 1e-12, at most 1e-11 of max\|c\| | 6.1e-12 | within the bound, above the central estimate |
| growth ∝ √t to ∝ t, no jump | no jump at the fixing step; the growth law is not identifiable (window spans a factor 2.2 in t) | partly |
| C++↔Fortran \|ΔP_reaction\| ≲ 1e-11 | 1.4e-11 (CPU), 1.6e-11 (GPU) | **no**, narrowly |
| C++↔Fortran within ~10× of CPU↔GPU | 12× (‖ΔS‖/‖S‖), 13× (P_reaction medians) | **no**, narrowly |
| … also per channel | specular inelastic channels 31–303× (S), 22–766× (BK) | **no** |
| printed set-up reals within printed precision | all but the two \|ψ(k)\| tails (337 and 346 half units, explained, item 11) | partly |
| printed comparisons within ½ unit except small values | 1003 values outside (CPU), all explained by a rule adopted after seeing them | partly |
| analysis line fixed at the same step | same step, t = 2220 a.u. | yes |
