# Output files

All files are written to `--output-dir`. Numbers that are compared with the Fortran code are
written in full precision (`%.16e`). Lengths are in bohr and times in atomic units everywhere;
energies are in eV in the result tables (`.tsv`, `analysis.json`) and in hartree in the `.dat`
files and `run_info.json`, unless a column or field name says otherwise (each file's `#` header names
the units).

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

## Run metadata

`run_info.json` is written by every run that reaches its output directory (all run modes), next to
the gas-phase files; an analysis-only run (`onlyAnalysis`) writes it into `fine_analysis/`, next to its
results, so that it never replaces the file of the propagation it analyses (`RestartDir` may be the
output directory). It holds what a plot or a reader needs to interpret the other files, in one
place. Atomic units (bohr, hartree, atomic units of time, electron masses) unless a field name ends
in a unit (`_ev`, `_amu`, `_deg`, `_rad`); `units` repeats this rule.

| Field | Content |
|---|---|
| `format_version` | 1 |
| `program` | `name`, `version`, `git`: the commit the `spodvr` executable was built from (`git describe --always --dirty`: `-dirty` when tracked files had uncommitted changes, new untracked files are not noticed; "unknown" outside a git checkout) |
| `title`, `system_id`, `input_file` | the input's `Title` and `SystemID`; the input file as given on the command line |
| `run_mode` | `code` (0 new, 1 restart, 2 analysis only, 6 potential only) and `name` |
| `new_calculation`, `restart_dir` | `NewCalc`; `RestartDir` as given (null if absent) |
| `times` | `t_start`, `t_end`, `dt_step`, `n_steps`, `dt_analysis`, `steps_per_analysis`, `dt_averages`, `analyses_per_average` |
| `molecule` | `mass1_amu`, `mass2_amu`, `homonuclear`, `total_mass`, `reduced_mass` (electron masses) |
| `initial_state` | `v`, `j`, `mJ`, `energy` and `energy_ev`: E(vin, jin) of the gas phase |
| `wave_packet` | `e_min_ev`, `e_max_ev` (`EMinWP`, `EMaxWP`: normal energies), `translational_energy_ev`, `z0`, `alpha` (the packet is (2α/π)^¼ exp(−α(Z − Z0)² − i kz Z)), `sigma_z` = 1/(2√α) (the width of \|ψ\|²), `k_total`, `kz`, `kx` and `ky` (components along the reciprocal lattice axes), `k_parallel` (Cartesian magnitude), `theta_incidence_rad`, `phi_incidence_rad`, `e_parallel_ev`, `phi_parallel_deg` |
| `analysis_energies` | `collision_ev` and `total_ev` (= collision + E(vin, jin)): the energies of the on-the-fly S matrix |
| `grids` | `z_start`, `dz`, `nz`, `nz_specular`, `r_start`, `dr`, `nr`, `nx`, `ny`, `dx`, `dy`, `j_max`, `mj_max` (the effective value: `MJMax` < 1 means `JMax`), `n_theta`, `n_phi`, `nz_potential` (Z points below z_asymptotic), `xy_volume_element`, `unit_cell_area`, `unit_cell_norm` |
| `rotational_basis` | `states`: [J, mJ] of every rotational state, in the order of the wave function's first index; `initial_index` |
| `lattice` | `a`, `b`, `skew_angle_rad`, `point_group`, `delta_kx`, `delta_ky` (reciprocal lattice spacings) |
| `z_asymptotic`, `z_asymptotic_unswitched` | where the gas phase begins (with switching: `ZmaxSwitch`), and the PES's own value |
| `analysis` | `z_ana` (on the grid) and its 0-based index `iz_ana`, `z_ana_specular` (`ZAnaSp`, null if not given), `dynamic_specular_analysis`, `fixed_z_analysis_specular` |
| `absorbers` | `power`, and `list`: per absorber its `name` (`ZOpt`, `ROpt`, `ZSpOpt`, and `ZOptLow`, `ZSpOptLow` when on), its `grid` (scattering, specular, both), `coordinate`, `start`, `end`, effective `strength` (hartree, after the energy table) and the grid `points` inside it (`first`, `last`, 0-based) |
| `potential` | `use_switching`, `z_min_switch`, `z_max_switch`, `impose_v_min`, `v_min`, `v_min_ev`, `impose_v_max`, `v_max`, `v_max_ev`, `use_symmetry`, `symmetry_tolerance_ev` |
| `unit_factors` | `hartree_to_ev`, `amu_to_me`, `bohr_to_angstrom`, `hartree_to_wavenumber`: the conversion factors the program uses |

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
