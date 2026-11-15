# Output files

All files are written to `--output-dir`. Numbers that are compared with the Fortran code are
written in full precision (`%.16e`). Lengths are in bohr and times in atomic units everywhere;
energies are in eV in the result tables (`.tsv`, `analysis.json`) and in hartree in the `.dat`
files, unless a column says otherwise (each file's `#` header names the units).

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

Closed channels (energetically inaccessible at that energy) are omitted from `s_matrix.tsv`;
in the summed tables they contribute 0. The diffraction order depends on the surface point group
(|n| + |m| for C2v and C4v).

`P_reaction` is not measured but derived: 1 − P_scattered. A run that ends before the scattered
packet has passed the analysis lines therefore reports P_reaction close to 1 at the energies whose
flux has not arrived yet. The means `E_vib_mean`, `E_rot_mean` and `E_dif_mean` are averages
over the scattered molecules (weighted by their probabilities and divided by P_scattered): where
P_scattered is tiny they are dominated by noise and mean nothing.

`analysis.json` holds the run metadata: program version, title, SystemID, units, initial state,
basis sizes (`n_vibrational_levels` = nviba, `j_max` = jana, `max_diffraction_order` = mdifa,
`n_rotational_states` = nana, `n_diffraction_channels` = ndifa, `n_states` = nstata), the
analysis surface, the energy range and the list of result files. A `CRPDistrtSurface` run adds
`scm`, its SCM identity, so that the runs of an ensemble can be told apart:

| Field | Content |
|---|---|
| `crp_data_set`, `parameter_file`, `position_file` | the keys `CRPlibDataSet`, `SCMParFilepath`, `PosFilepath` as given |
| `slab_tag`, `position_file_version` | the first field of line 1 of the position file; its version (3 or 2) |
| `slab_number`, `cutoff`, `rescaling` | `PosFileSlabNum`, `SCMCutoff` (bohr), `EnableSurfExpan` |
| `slab_lattice_constant`, `slab_atoms` | a_T (bohr, line 3 of the position file); the atoms of layers 1–2 that ΔV uses |
| `rms_displacement`, `top_layer_shift`, `highest_atom` | of those atoms (bohr): the rms of distorted − ideal, its mean [x, y, z] over layer 1, the largest distorted z |
| `parameters` | P1–P18 of layer 1, then of layer 2, as read from the parameter file |

The set-up output of a `CRPDistrtSurface` run has, after the CRP data set, the block "SCM DISTORTED
SURFACE": the files, the snapshot and its atoms per layer, the lattice constants a_T and a0 and the
rescaling s, the displacements (rms, largest, the mean shift of layer 1, the highest atom), and the
checks of [fortran_deviations.md](fortran_deviations.md#scm-distorted-surface-pes-srcpotentialscm-systemid-crpdistrtsurface):
the ideal atoms against the CRP lattice (SCM D9), V − V_CRP at the gas-phase points (V_asy and,
with switching, the switch's V_gas) and on the plane Z = z_asymptotic (SCM D3; a changed gas-phase
point stops the run), a table of the jump of V across the cell edges (SCM D10), the lowest V − V_CRP
on the Z planes below the CRP's Zmin at Rmin ≤ r ≤ Rmax (SCM D15), each with its warning, and the
lowest V − V_CRP at r outside Rmin .. Rmax on every Z plane, with whether it lies below Rmin or
inside the r absorber (SCM D15, no warning). An analysis-only run checks only the gas-phase points.

## Data files (space-separated)
The `.dat` files hold plain columns of numbers separated by spaces, after `#` comment lines that
name the columns and units; `numpy.loadtxt(f)` reads them.

| File | One row per | Columns |
|---|---|---|
| `averages.dat` | averages time (every dtAve) | t (a.u.), norm, ⟨Z⟩ (bohr), ⟨r⟩ (bohr), rest norm (the Fortran's averages.out) |
| `asymptotic_potential.dat` | r grid point | index ir, r (bohr), V_asy (hartree), V_asy (eV): the gas-phase potential at Z = zasymp, θ = π/2 |
| `gas_phase_states.dat` | rovibrational state (v, j) | v, j, E (hartree), ψ(r) at every r grid point (a `#` line lists the grid) |

## Correspondence with the Fortran output
| Fortran | C++ |
|---|---|
| stdout per-energy blocks (`smatrix_combine`) | screen summary tables; `scattering_summary.tsv`, `probabilities_rovib.tsv`, `probabilities_diffraction.tsv` |
| fort.25 (state-to-state S-matrix, `writeProbs`) | `s_matrix.tsv` |
| fil10 / unit 3 (fine analysis, spectrum) | `fine_analysis/spectrum/` |
| stdout blocks and fort.25 of a fine analysis (interpolated) | `fine_analysis/scattering_summary.tsv`, `s_matrix.tsv`, … |
| averages.out | `averages.dat` |
| fort.9 / fort.8 (wave function) | `checkpoint.bin` (or fort.9 with `CheckpointFormat`) |
| fort.21 / fort.22 (BK coefficients) | `bk_coefficients.bin` (or fort.21 with `CheckpointFormat`) |

## Fine analysis (`onlyAnalysis`, `alsoAnalysis`)
Written to `fine_analysis/` in the output directory, in the same table formats as above:

| Path | Content |
|---|---|
| `fine_analysis/spectrum/` | the four tables + analysis.json at every point E_k = k π / T of the spectral grid inside the 95% confidence interval (Fortran fil10 and unit 3). `E_z` is k·spres; `s_matrix.tsv` only with `writeSMatrix` |
| `fine_analysis/*.tsv`, `analysis.json` | the spectrum interpolated linearly to the NrOfEs analysis energies (Fortran: stdout and fort.25 of the analysis run). Energies outside the grid or the confidence interval are left out |
| `fine_analysis/grid.json` | records in the file / up to tEnd / used (ndused), transform length, time step, effective time, spectral resolution (eV), first and last grid index, number of points (nen95), skipped analysis energies with the reason, and the count of large phase steps |

The interpolated S-matrix elements and phases differ from the Fortran fort.25 because two Fortran
bugs are fixed; the probabilities are the same. See `docs/fortran_fine_analysis_bugs.md`.

## Checkpoint files (binary)
Native files are in host byte order (little-endian on x86-64). Both start with a header: a
16-byte magic string (`SPODVR-CHECKPNT` / `SPODVR-BKCOEFFS`), a `uint32` format version (2; version 1
lacked the number of phi points), the
SystemID (`uint64` length + characters) and the list of restart parameters (`uint64` count, then
name and `double` value per entry). A restart compares this list with its own input. A
`CRPDistrtSurface` run appends its SCM identity: `PosFileSlabNum`, `SCMCutoff`, `EnableSurfExpan`
(0 or 1), `PosFilepath (a_T)`, `PosFilepath (atoms used)`, `PosFilepath (rms displacement)`,
`PosFilepath (layer 1 mean shift x)` (and `y`, `z`), `PosFilepath (highest atom)` and
`SCMParFilepath (layer L, Pk)` for L = 1, 2 and k = 1 … 18; the files of other PESs keep their list.

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
