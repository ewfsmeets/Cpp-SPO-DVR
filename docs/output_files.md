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
| `plot_data` | only when a plot-data key is set: `write_densities`, `write_potential`, `dt_frames` (null without frames), `sites` (per site of `PlotSites` its `name`, `x`, `y` and its grid point `ix`, `iy`), `norm_budget` (true; only with `NormBudget`) and `frame_extras` (the names; only with `FrameExtras`) |
| `pes_cut` | only in a `write2DPES` run: `file` (`pes_cut_2d.dat`), `x`, `y` (bohr), `theta_deg`, `phi_deg`. Such a run stops before the gas phase: `initial_state` `energy` and the `analysis_energies` are null; `run_mode` stays the input's (new, restart, …) although the run does none of it |

## Data files (space-separated)

The `.dat` files hold plain columns of numbers separated by spaces, after `#` comment lines that
name the columns and units; `numpy.loadtxt(f)` reads them.

| File | One row per | Columns |
|---|---|---|
| `averages.dat` | averages time (every dtAve) | t (a.u.), norm, ⟨Z⟩ (bohr), ⟨r⟩ (bohr), rest norm (the Fortran's averages.out) |
| `asymptotic_potential.dat` | r grid point | index ir, r (bohr), V_asy (hartree), V_asy (eV): the gas-phase potential at Z = zasymp, θ = π/2 |
| `gas_phase_states.dat` | rovibrational state (v, j) | v, j, E (hartree), ψ(r) at every r grid point (a `#` line lists the grid) |
| `pes_cut_2d.dat` | (r, Z) point of a `write2DPES` run | r (bohr), Z (bohr), V (hartree), V (eV): the PES itself (no switching, no VMin/VMax) at the fixed x, y, θ, φ that the `#` lines give; Z fastest over the scattering grid's nZ points, a blank line after each r (the Fortran's order). The run writes it with `asymptotic_potential.dat` and `run_info.json`, then stops |

## Series directories (binary, plot data)

A series directory holds arrays recorded during the run, for the plot tool (`tools/plots/README.md`).
It is written only when its key is set (input_schema.md, "Plot data keys"), and it holds:

- `manifest.json`: `format_version` (1), `kind` (the directory's name), `byte_order` (`little`),
  `memory_order` (C order: the last axis fastest), `time_dependent`, `times` (file, dtype, units) and,
  under `datasets`, per dataset its `file`, `dtype` (`<f8` double, `<i4` 32-bit integer), the `shape` of
  one record (`[]` for a scalar), `axes`, `units` and `description`; then the fields of the series (grids,
  record interval);
- `times.bin`: the time of each record (`<f8`, atomic units);
- one `<dataset>.bin` per dataset (the dataset's name is the file name without `.bin`): its records
  one after another, without a header.

The records are appended during the run: for each record every dataset file is written and flushed
first, the time last. A reader takes n = min(number of times, size of each data file / bytes per
record) records, so a run that is still going, or that stopped, never gives a time without its data.
In numpy, for the dataset `prob_scat_zr` of the directory `d`:

```python
cells = np.fromfile(d / "prob_scat_zr.bin", "<f8").reshape(-1, nz, nr)[:n]
```

A restart writes the records of its own segment into its own output directory, from its tStart. When
that time is on the record grid of the earlier run, its first record repeats that run's last one (the
same wave function: bitwise equal on the same device and build, equal to rounding when the restart runs on
another device); the plot tool joins the segments and drops the repeat after checking that it equals the earlier
record to 1e-12 of its largest value.

### `densities/` (`writeDensities`)

One record at tStart and at every dtAve, the times of `averages.dat`:

| File | Type and shape | Meaning |
|---|---|---|
| `prob_scat_zr.bin` | `<f8` [nZ, nR] | probability per (Z, r) cell of the scattering wave, summed over x, y and the rotational states: Σ\|ψ\|²·dZ·dr·dx·dy·sin(skew) |
| `prob_spec_zr.bin` | `<f8` [nZSp, nR] | probability per (Z, r) cell of the specular wave χ (the initial channel, incoming and specularly reflected, every v): \|χ\|²·dZ·dr |
| `surface_index.bin` | `<i4` [] | 0-based Z index of the specular analysis surface (it moves until it is fixed, `DynamZSpAna`) |
| `norm.bin` | `<f8` [] | the norm of the `averages.dat` row at this time |

The manifest adds `record_interval_au` (dtAve) and `grids` (`z_start`, `dz`, `nz`, `nz_specular`,
`r_start`, `dr`, `nr`; both Z grids start at `z_start` with the same spacing). The records are the
numbers the averages are computed from, so the sum of both datasets over all cells is the norm (to
rounding), and writing them adds no computation to the time loop. A probability per cell divided by
dZ·dr is a density in bohr⁻².

### `frames/` (`dtFrames`)

One record at tStart and at every dtFrames, of the whole wave packet ψ_full = ψ_scat + s2s·χ, the scattering
wave with the specular wave χ added in the initial rotational channel (s2s(x, y) = exp(i k∥·R)·unit-cell norm,
the factor of the two-grid projection), reduced over the axes each dataset leaves out:

| File | Type and shape | Meaning |
|---|---|---|
| `prob_zrxy.bin` | `<f8` [nZ, nR, nX, nY] | probability per (Z, r, x, y) cell of ψ_full, summed over the rotational states: Σ\|ψ_full\|²·dZ·dr·dx·dy·sin(skew) |
| `prob_tail_zr.bin` | `<f8` [nZSp − nZ, nR] | probability per (Z, r) cell of χ beyond the scattering grid (its first row is Z index nZ): \|χ\|²·dZ·dr |
| `pop_jmj_z.bin` | `<f8` [nZ, nrot] | probability per (Z, rotational state) of ψ_full, summed over r, x, y; the states in the order of `rotational_states` |
| `surface_index.bin` | `<i4` [] | 0-based Z index of the specular analysis surface |

With `FrameExtras` the record holds, after those, the datasets of the extras that are on (in this order):

| File | Extra | Type and shape | Meaning |
|---|---|---|---|
| `prob_angle_z.bin` | `angular` | `<f8` [nZ, nθ, nφ] | probability per DVR node (θ_t, φ_k) of ψ_full, summed over r, x, y; divided by w_t·2π/nφ (`theta_weights`) a density per steradian (φ from the first lattice vector; for a homonuclear molecule the θ grid covers [0, π/2]: the other hemisphere is its mirror (π − θ, φ + π), and the weights of the nodes below π/2 count their mirror nodes too, the node at π/2, whose mirror lies in its own row, not) |
| `pop_nm_z.bin` | `diffraction` | `<f8` [nZ, nX, nY] | probability per parallel momentum K = k∥ + n b₁ + m b₂ of ψ_full (the x-y FFT of exp(−i k∥·R) ψ_full, the BK analysis's convention), indices in FFT order (n for n ≤ nX/2, else n − nX, as the analysis labels the orders), summed over r and the rotational states; only open orders at Z ≥ zasymp are diffraction channels |
| `pop_vj_z.bin` | `vibration` | `<f8` [nZ, nR, nJ] | probability per gas-phase state φ_vJ(r), all nR of each J of the basis (`vibration_j_values`), summed over x, y and mJ: a projection on the asymptotic states, also where the molecule is not free |
| `p2_numerator_z.bin`, `rot_numerator_z.bin` | `rotation_moments` | `<f8` [nZ] | ⟨P₂(cos θ)⟩ times the probability at each Z (analytic J → J, J ± 2 elements), and ⟨J(J+1)/(2μr²)⟩ (hartree) times it: divided by the probability at that Z, the averages |
| `prob_kr.bin`, `prob_kz_spec.bin` | `momentum` | `<f8` [nR], [nZSp] | probability per k_r point of ψ_full on the scattering grid (FFT order, spacing 2π/(nR dr)), and per k_Z point of the specular wave χ (FFT order, spacing 2π/(nZSp dZ)) |

At each Z the first three sum to the frame's probability at that Z (`pop_jmj_z` summed over the states);
`prob_kr` sums to the probability on the scattering grid (`prob_zrxy` summed), `prob_kz_spec` to the norm of the
specular wave (Parseval, and completeness for `vibration`). The manifest adds
`record_interval_au` (dtFrames), `grids` (`z_start`, `dz`, `nz`, `nz_specular`, `r_start`,
`dr`, `nr`, `nx`, `ny`), `rotational_states` ([J, mJ] per state), `initial_index` (the state that holds χ),
`sites` (as in `potential/`) and, with `FrameExtras`, `extras`: `names`, with `angular` its `theta_rad`, `theta_weights`, `phi_rad` and
`homonuclear`, with `vibration` its `vibration_j_values`. The sum of `prob_zrxy` and `prob_tail_zr` is the norm (ψ_scat holds no
initial-channel component: the cross terms with χ vanish in the sum over x, y, not point by point).

### `norm_budget/` (`NormBudget`)

One record per analysis window (t − dtAna, t], at its end time t, from tStart + dtAna on (a restart's first window
starts at its tStart, so no record repeats the run before): what each absorber removed in the window.

| File | Type and shape | Meaning |
|---|---|---|
| `loss.bin` | `<f8` [n_absorbers] | each absorber's loss in the window, summed over its cells |
| `loss_error.bin` | `<f8` [n_absorbers] | an error indicator of `loss`, summed over the same cells |
| `norm_start.bin` | `<f8` [] | the probability on both grids at the window's start |
| `norm_end.bin` | `<f8` [] | the probability on both grids at the window's end |

For each (Z, r) cell of a grid, ρ is its probability (the `densities/` cells, from the same reduction) at the
window's start and end. The absorbers acting on the cell shrink its density by f² per time step (f = exp(−V_opt dt),
the propagation's own tables), F = exp(−γ) over the window's n steps, γ = 2n V_opt dt. With a constant inflow J
into the cell during the window (dρ/dt = −κρ + J), the cell loses exactly

```text
loss = ρ_start (1 − F) + (ρ_end − F ρ_start) (γ/(1 − F) − 1),
```

which reduces to ρ_start (1 − F) without motion and to ρ γ for a steady cell. The indicator is γ²/6 of the cell's
inflow ρ_end − ρ_start + loss: the leading distance between the loss and the mean of the two one-sided estimates
ρ_start (1 − F) and ρ_end (1/F − 1), so it falls with dtAna squared. Where a Z absorber and the r absorber overlap, the
cell's loss and indicator are shared in the ratio of their V_opt. The manifest adds `record_interval_au` (dtAna),
`window` (`steps` n, `dt_step`) and `absorbers`, one object per entry of the loss datasets: `name` (`ZOpt`, `ZOptLow`,
`ROpt`, `ZSpOpt`, `ZSpOptLow`) and `grid` (`scattering` or `specular`; `ROpt` acts on both and has one entry per
grid).

Summed over the absorbers and windows, the losses account for the norm that left the grids, norm_start of the first
window minus norm_end of the last, to within the summed indicator and the time stepping's own error (the absorbers act
once per time step, between the kinetic half steps, not continuously: on the hard-wall test grid, where by t = 6000
a.u. the r absorber has eaten the tail of the bound state on the specular grid, 4.1e-7 of the norm at dtStep 2 and
1.05e-7 at dtStep 1). The plot tool shows the residual against the indicator, and its verdict allows 1e-4 of the
initial norm beyond it. A larger residual is probability
that left another way. Apart from a diverging run, the only other non-unitary step
of the propagation is the angular transform at Z < zasymp: its DVR → FBR keeps only the rotational basis's part of
exp(−iV dt)ψ (`AngularTransform.AnisotropicPhaseLosesTheNormOutsideTheBasis`), the Fortran's method as well; the
budget's residual is attributed to it. On a 4 × 4 Dai–Light grid with JMax 2 the norm fell by 4.2e-3 by t = 6000
a.u. with every absorber switched off (4.6e-3 at normal incidence, 4.5e-3 of residual with JMax 4), while the packet
was at the surface; halving dtStep reduced it to 3.2e-3. The hard wall, without angular coupling, loses nothing.
That probability counts as reaction in P_reaction = 1 − P_scattered.

### `potential/` (`writePotential`)

Written once at set-up (also by a potential-only run), not a time series: `manifest.json` has `time_dependent`
false and no `times.bin`, and each dataset file holds one record. The potential is the one the propagation uses
(switching to the gas phase and the VMin/VMax limits included), at the grid points of the interaction region
Z < zasymp (`nz_potential` points), in hartree, read through the same surface-symmetry maps as the propagation.

| File | Type and shape | Meaning |
|---|---|---|
| `v_min.bin` | `<f8` [nZpot, nR, nX, nY] | the minimum of V over the DVR orientations (the θ and φ nodes) at each grid point: the lower envelope over the run's nodes (each point its own orientation; no rotational or zero-point energy) |
| `v_min_node.bin` | `<i4` [nZpot, nR, nX, nY, 2] | the orientation of that minimum, (iθ, iφ): indices into `theta_nodes_rad` and `phi_nodes_rad` (the first in iθ, then iφ order; a NaN anywhere gives NaN and its node) |
| `v_sites.bin` | `<f8` [nSites, nZpot, nR, nθ, nφ] | V at every orientation at the (x, y) grid point of each site of `PlotSites` (only with sites) |

The manifest adds `grids` (`z_start`, `dz`, `nz_potential`, `r_start`, `dr`, `nr`, `nx`, `ny`),
`theta_nodes_rad`, `theta_weights`, `phi_nodes_rad`, `homonuclear` (θ ≤ π/2 only), `sites` (per site `name`,
`x`, `y`, `ix`, `iy` and `ixy` = iy + nY·ix), `limits` (switching and the VMin/VMax settings, the limits in
hartree) and `hartree_to_ev`. The files are written one (Z, r) point at a time: no memory beyond one row (no
memory-estimate item). With the full grid of a run, `v_min.bin` takes nZpot·nR·nX·nY·8 bytes.

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
| 2dpes.dat (`write2DPES`: r, Z, V in eV) | `pes_cut_2d.dat` (r, Z, V in hartree and eV) |

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
