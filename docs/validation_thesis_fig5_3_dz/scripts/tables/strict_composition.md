| Propagation only (replay vs replay, relative to \|S\|) | C++ CPU fast − F | C++ GPU fast − F | C++ GPU strict − F | GPU strict − GPU fast | Criterion (2026-11-08) |
|---|---|---|---|---|---|
| κ, uniform amplitude part (median over E) | 6.0e-12 | 6.6e-12 | 1.5e-13 | -6.4e-12 | fast math if \|κ\| ≤ 1.5e-12, not if ≥ 4.5e-12 |
| κ range over E | 5.4e-12 … 6.9e-12 | 6.0e-12 … 7.8e-12 | -3.8e-14 … 4.1e-13 | -7.4e-12 … -5.8e-12 |  |
| BK amplitude part (window median) | 6.7e-12 | 7.3e-12 | 2.0e-13 | -7.1e-12 | second estimator |
| residual outside the family (norm-weighted) | 4.2e-12 | 4.2e-12 | 2.5e-13 | 4.2e-12 |  |
| family j = 10, residual (norm-weighted) | 8.4e-9 | 8.4e-9 | 9.5e-10 | 9.3e-9 | fast math involved if ≤ 8.4e-10, not if ≥ 4.2e-9 |
| family j = 10, per element (max) | 1.6e-6 | 1.6e-6 | 1.8e-7 | 1.8e-6 |  |
| family j = 8, per element (max) | 8.0e-7 | 8.0e-7 | 3.1e-8 | 7.7e-7 |  |
| S, norm-weighted relative difference | 7.5e-12 | 7.9e-12 | 3.5e-13 | 7.8e-12 | strict − fast predicted ≈ 1e-12 |
| \|ΔS\| max | 1.7e-12 | 2.2e-12 | 2.5e-13 | 2.0e-12 |  |
| \|ΔP_reaction\| max | 1.4e-11 | 1.6e-11 | 8.2e-13 | 1.5e-11 |  |
