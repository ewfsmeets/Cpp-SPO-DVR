# Output files

All files are written to `--output-dir`. Numbers that are compared with the Fortran code are
written in full precision (`%.16e`). Energies are in eV, lengths in bohr, times in atomic units.

## Results (tab-separated tables)
Every table starts with `#` comment lines that explain the columns and units, followed by one
header line and one row per record ("long" format: one quantity per row, labels in columns).
```python
import pandas as pd
s = pd.read_csv("s_matrix.tsv", sep="\t", comment="#")
s[(s.v == 0) & (s.j == 2)].groupby("E_coll").P.sum()   # P(v=0, j=2) summed over mJ, n, m
```
gnuplot: `set datafile separator "\t"; plot "scattering_summary.tsv" using "E_coll":"P_reaction"`.

| File | One row per | Columns |
|---|---|---|
| `scattering_summary.tsv` | analysis energy | `E_coll` collision energy, `E_total` = E_coll + E(vin, jin), `E_z` = E_coll − initial parallel energy, `P_scattered` (Σ\|S\|²), `P_reaction` (1 − P_scattered), `P_v0`, `P_v1`, … (summed over rotation and diffraction), `E_vib_mean` `E_rot_mean` `E_dif_mean` (average final energies of the scattered molecules) |
| `s_matrix.tsv` | energy and **open** final channel | `E_coll`, `E_total`, `v`, `j`, `mJ`, `n`, `m` (diffraction channel n b₁ + m b₂), `order`, `Re_S`, `Im_S`, `P` = \|S\|², `phase` = arg S (rad; 0 when Re S = 0, as in Fortran) |
| `probabilities_rovib.tsv` | energy and (v, j, mJ) | `E_coll`, `v`, `j`, `mJ`, `P` (summed over diffraction). Written if `writeProbs` |
| `probabilities_diffraction.tsv` | energy, (v, j, mJ) and diffraction order | `E_coll`, `v`, `j`, `mJ`, `order`, `P` (summed over the (n, m) of that order). Written if `writeRDProbs` |
| `averages.dat` | averages time | t, norm, ⟨Z⟩, ⟨r⟩, rest norm (space-separated, Fortran averages.out) |

Closed channels (energetically inaccessible at that energy) are omitted from `s_matrix.tsv`;
in the summed tables they contribute 0. The diffraction order depends on the surface point group
(|n| + |m| for C2v and C4v).

`analysis.json` holds the run metadata: program version, title, SystemID, units, initial state,
basis sizes (`n_vibrational_levels` = nviba, `j_max` = jana, `max_diffraction_order` = mdifa,
`n_rotational_states` = nana, `n_diffraction_channels` = ndifa, `n_states` = nstata), the
analysis surface, the energy range and the list of result files.

Set-up files: `asymptotic_potential.dat` (vasy(r)) and `gas_phase_states.dat` (rovibrational
energies and wave functions).

## Correspondence with the Fortran output
| Fortran | C++ |
|---|---|
| stdout per-energy blocks (`smatrix_combine`) | screen summary tables; `scattering_summary.tsv`, `probabilities_rovib.tsv`, `probabilities_diffraction.tsv` |
| fort.25 (state-to-state S-matrix, `writeProbs`) | `s_matrix.tsv` |
| fil10 / unit 3 (fine analysis) | (next milestone, same table formats) |
| averages.out | `averages.dat` |
| fort.9 / fort.8 (wave function) | `checkpoint.bin` (or fort.9 with `CheckpointFormat`) |
| fort.21 / fort.22 (BK coefficients) | `bk_coefficients.bin` (or fort.21 with `CheckpointFormat`) |

## Checkpoint files (binary)
Native files are in host byte order (little-endian on x86-64). Both start with a header: a
16-byte magic string (`SPODVR-CHECKPNT` / `SPODVR-BKCOEFFS`), a `uint32` format version (1), the
SystemID (`uint64` length + characters) and the list of restart parameters (`uint64` count, then
name and `double` value per entry). A restart compares this list with its own input.

`checkpoint.bin` (written to `checkpoint.bin.tmp`, then renamed):
| Field | Type |
|---|---|
| header | see above |
| time | `double` |
| n_rot, n_xy, n_zr, n_r, n_zsp | 5 × `int32` |
| analysis surface: index, previous index (0-based), fixed, anlyTol | `int32`, `int32`, `uint8`, `double` |
| specular wave function (ir fastest) | n_r · n_zsp × `complex<double>` |
| scattering wave function, layout (irot, ixy, izr) | n_rot · n_xy · n_zr × `complex<double>` |
| end marker | `uint64` |

`bk_coefficients.bin` (appended and flushed every analysis step):
| Field | Type |
|---|---|
| header | see above |
| nstata, nviba | 2 × `int32` |
| labels j, mJ, n, m, order of every state | nstata × 5 × `int32` |
| channel energies (hartree), (is, v) with is fastest | nstata · nviba × `double` |
| per analysis step: time, coefficients (is, v) | `double`, nstata · nviba × `complex<double>` |

The Fortran-layout files fort.9 and fort.21 (gfortran unformatted sequential records) are
described in `src/io/fortran_files.hpp`.
