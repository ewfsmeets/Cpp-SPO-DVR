# The method as the code implements it

What the program computes and in which order, with the C++ class or function that does each part and
its Fortran routine. The theory is that of M. F. Somers, *Six-Dimensional Quantum Dynamics of
Dissociative Chemisorption of Hydrogen on Metal Surfaces*, PhD thesis, Leiden University (2004) (README,
"How to cite"); equation numbers below refer to its chapter 2. How the code is organised:
[developer_guide.md](developer_guide.md). Where the C++ deliberately differs from the Fortran:
[fortran_deviations.md](fortran_deviations.md).

## Coordinates and units

- The molecule's centre of mass (x, y, Z) and its bond (r, θ, φ). x and y are **skewed** coordinates
  along the surface lattice vectors a1 = alat (1, 0) and a2 = blat (cos γ, sin γ), γ the skew angle
  (90° for square cells, 60° for hexagonal ones); Z is the height above the surface plane; θ is
  measured from the surface normal, φ from a1 towards a2 (`potential/potential.hpp`,
  `symmetry/surface_lattice.hpp`).
- The PESs evaluate V with the atoms at ±r/2 about (x, y, Z), so for unequal masses (HD) the point
  the kinetic energy refers to and the point the PES sees differ (fortran_deviations.md, SCM D13).
- Atomic units inside the program; the input is in bohr, eV, amu, degrees and a.u. of time.
  Conversions go through `units::legacy` (the factors the Fortran uses almost everywhere: 27.2114 eV
  per hartree, 1822.89 mₑ per amu, 0.529177 Å per bohr) or `units::precise` (the Fortran's flux
  analysis, which is not ported, and the SCM correction), never through a factor written into the
  code (`core/unit_conversion.hpp`).
- Every PES returns V in eV.

## Hamiltonian and grids

The nuclear Hamiltonian of a diatomic molecule above a static surface (thesis (2.9)–(2.12)):

  H = −∇²_R / 2M + H_vib(r) + H_rot(r, θ, φ) + V(x, y, Z, r, θ, φ),

with M the total mass, µ the reduced mass, H_rot = j² / (2µ r²), and H_vib = −(1/2µ) ∂²/∂r²: the
stored function is r ψ, so T_r has this simple form.

- Z, r, x, y: uniform grids with FFT kinetic operators (thesis (2.55)). The x and y momenta are the
  components along the reciprocal lattice vectors, k_x = k_x,0 + n · 2π/(alat sin γ) and likewise
  k_y, shifted by the incident parallel momentum (`input/process.cpp`); with them
  T_xy = (k_x² + k_y² − 2 cos γ k_x k_y) / 2M (`propagation/propagation_tables.cpp`).
- θ, φ: a discrete variable representation (DVR), Gauss–Legendre points in θ and an FFT grid in φ;
  the finite basis (FBR) is the rotational basis |J mJ⟩ (`grids/rotational_basis.hpp`; for
  homonuclear molecules J in steps of 2). The transform between them is a Legendre sum per mJ block
  and an FFT over φ (`AngularTransform`; thesis (2.61)–(2.63); Fortran: `transform.f90`).
- The wave function on the **scattering grid** is stored as (rotational state, x y point, r Z point),
  first index fastest. Rotation, absorbers and the gas phase are applied in the FBR; the potential in
  the DVR.

## Two grids

The incoming packet starts far from the surface. To keep the large scattering grid short, the
program carries the incoming channel on a second, longer (Z, r) grid, the **specular grid**
(Fortran: `TwoGrids.f90`, ProjOp = 2, the only option the C++ supports).

- The specular wave χ(r, Z) holds the component of the initial rotational state (j_in, mJ_in) in the
  initial diffraction channel (0, 0), for every vibrational state. It starts as the initial packet
  (`initial_wave_packet.cpp`; thesis (2.33)–(2.34): a Gaussian in Z times the rovibrational state)
  and is propagated with its own kinetic operator, the gas-phase V(r) beyond the interaction region
  and its own absorbers (`SpecularOperator`).
- Inside the interaction region (Z below z_asymptotic) χ feels the PES only through the scattering
  grid: it is added to the initial channel of the scattering grid before the potential step and
  projected back out after it, which replaces χ there (`PotentialOperator::apply_block`; Fortran:
  `AddSpecToScat`, `SubSpecFromScat`). The initial channel thus lives on the specular grid
  everywhere; the scattering grid holds every other channel.
- The split is exact only up to the grids: the two Z grids are periodic over different lengths, and
  the specular grid can have a low-Z absorber of its own (`ZSpOptLow*`).

## One analysis step

The time loop advances one analysis step of n time steps per call (`WavePacketPropagation`,
`propagation/propagation.cpp`), n = int(dtAna / dtStep), truncated as in the Fortran
(`input/process.cpp`). The split-operator propagator (thesis (2.35)) with adjacent half steps of T
merged is

  [T(dt/2) V(dt) T(dt/2)]ⁿ = T(dt/2) V (T V)ⁿ⁻¹ T(dt/2),

the list `analysis_step_schedule` (`propagation/split_operator_schedule.hpp`; Fortran:
`SplitOperator.f90 Propagate`), which the CPU and the GPU backend both run. Each entry is one T block,
followed by one V block except for the last.

T is T_Z + T_r + T_xy only. The rotational energy j²/2µr² is diagonal in (r, J), so it is applied in
the V block as rot(dt/2) · exp(−iV dt) · rot(dt/2), together with the absorbers and the gas phase.

| Stage | What it does | C++ | Fortran |
|---|---|---|---|
| specular grid, T | 2D FFT in (Z, r), the kinetic factors, back | `SpecularOperator::apply_kinetic` | `KineticOp.f90 PerformEXPKinOpSpec` |
| specular grid, its V part | the gas phase beyond the interaction region (rot(dt/2) exp(−i V(r) dt) rot(dt/2)), then the Z and r absorbers; skipped on the final half step, and run before the scattering grid's T of the same block | `SpecularOperator::apply_gas_phase`, `apply_absorbers` | `GasPhase.f90 PerformEXPGasPhaseOnSpecGrid`, `Absorption.f90` |
| T_Z + T_r | one 2D FFT per work unit of the scattering wave | `ScatteringKineticOperator::apply_zr` | `KineticOp.f90` |
| T_xy | per (Z, r) point: the momentum shift of oblique incidence, 2D FFT in (x, y), exp(−i T_xy dt), back, unshift | `ScatteringKineticOperator::apply_xy_group` | `SplitOperator.f90 PerformEXPKineticOp` |
| V, Z ≥ z_asymptotic | absorbers and the free molecule, one scalar factor per (rotational state, r) | `PotentialOperator::apply_block` | `PerformEXPOpticalZ/R`, `PerformEXPGasPhaseOnScatGrid` |
| V, Z < z_asymptotic | the absorbers; the specular wave added to the initial channel; rot(dt/2) (the code merges these factors, so the added term carries its own rot(dt/2)); FBR → DVR, exp(i A), DVR → FBR; rot(dt/2); the initial channel projected back onto the specular grid | `PotentialOperator::apply_block`, `apply_interaction` | `Potential.f90 PerformEXPPotentialOp`, `PerformEXPInteractionOp` |

- The CPU fuses T_xy and V per group of (Z, r) points (`SplitOperatorPropagator::apply_xy_and_potential`),
  which gives the same bits as the two stages one after the other.
- The stored potential is A = −V dt, on the DVR grid and only for Z below z_asymptotic (the
  `nz_potential` points; `potential/potential_grid.cpp`); exp(i A) is precomputed when that table fits
  in 4 GiB, else evaluated in every step, with the same results.
- With a switching function the potential is blended into the gas-phase V(r) between ZminSwitch and
  ZmaxSwitch, and z_asymptotic becomes ZmaxSwitch (`potential_grid.cpp`, `input/process.cpp`;
  Fortran: `switchPotentialToGasPhase`).

## FFT conventions

The Fortran calls the transform from coordinates to momenta "Backward"
(`KineticOp.f90 TransToMomenta`: `BackwardFFTx`, `BackwardFFTy`). In the C++ it is `FftSign::Minus`
with the scale d/√(2π) per Z or r dimension and √(d/n) per x or y dimension, and the way back is
`FftSign::Plus` with the scale 1/(n · scale) (`make_propagation_tables`, `ScatteringKineticOperator`).
`fft/fft.hpp` states the Fortran's sign convention.

## Absorbers

The Z ends of both grids and the large-r end carry quadratic negative imaginary potentials
(thesis (2.46)), applied as exp(−V_opt dt) factors in the V block (on the specular grid in its V part,
after its T). Their strengths follow from a target energy and the absorber length through the table
of Vibók and Balint-Kurti (`optical_potential_strength`, `input/wave_packet_params.cpp`; Fortran:
`MyFunctions.f90`). Since 5.24 there are optional absorbers at the low-Z end (`ZOptLow*`,
`ZSpOptLow*`). How to choose them, and the error they leave: [inputs/README.md](../inputs/README.md),
sections 7 and 8.

## Asymptotic analysis and the S matrix

The scattering amplitude formalism of Balint-Kurti (thesis (2.37)–(2.43)):

- **Coefficients.** After every analysis step the wave packet is projected onto the asymptotic
  channels (v, j, mJ, n, m) at an analysis line: the scattering wave at Z = `Zana`, and the specular
  wave, which holds the initial rotational and diffraction state for every v, at its own line Z_sp on
  the specular grid once that line is fixed; until then those coefficients are zero
  (`BalintKurtiAnalyzer`, `analysis/balint_kurti.cpp`; thesis (2.41); Fortran:
  `Analyzor.f90 BalintKurtiCoeffs`).
- **The specular line.** With `DynamZSpAna` it starts at the specular absorber and follows the tail of
  the incoming packet inwards: it is the outermost point with less than a tolerance of norm between it
  and the absorber (the tolerance is max(`TolSp`, the initial norm in the absorber)). It is fixed when
  it would have to move outwards again (the reflected packet arrives) or when it reaches `Zana`
  (`AnalysisSurface::update`; Fortran: `FindNewAnalysisPoint`).
- **Time to energy.** A(E) = Σ_t C(t) e^{iEt} Δt_ana, a rectangle sum over the analysis steps with
  Δt_ana = n · dtStep (`SMatrixAccumulator::add`; thesis (2.28), (2.40)).
- **S matrix.** For an open channel with perpendicular momentum k_z at total energy E
  (`SMatrixAccumulator::result_at`; thesis (2.39); Fortran: `smatrix_calc`):

    S = A(E) · (α / 2π)^{1/4} / M · exp(i (−k_z Z_m + (k_z,0 − k_0) Z_0 − π)) · √(k_z k_0 / w),

  with k_0 the incident perpendicular momentum at E, k_z,0 and α the packet's central momentum and
  width parameter, Z_0 its start, w = exp(−(k_z,0 − k_0)² / 2α) the packet's weight at k_0, and Z_m
  the line where the coefficient was taken (`measurement_line`): Z_sp for the initial rotational and
  diffraction state once its line is fixed (and the run's checkpoint is at hand), `Zana` otherwise and
  for all other channels. The Fortran uses `Zana` for every channel
  ([analysis D2](fortran_deviations.md#analysis-d2)).
- **Probabilities.** P = |S|² (thesis (2.42)). The reaction probability is 1 − Σ P over the analysed
  open channels (thesis (2.43)), so it counts everything not scattered into an analysed channel:
  dissociation (taken out by the large-r absorber), flux into states outside the analysis basis,
  whatever the low-Z absorbers take out, and norm that has not crossed an analysis line by `tEnd`.
  The accuracy of the probabilities: [inputs/README.md](../inputs/README.md), section 8.

## Fine analysis

A run writes its BK coefficients; the fine analysis (`onlyAnalysis`, or `alsoAnalysis` after a run)
transforms the whole series at once: a zero-padded FFT onto the spectral grid E_k = kπ/T (E measured
from the initial channel energy), kept within the 95 % energy range of the packet, then linearly
interpolated to the `NrOfEs` analysis energies (`analysis/fine_analysis.hpp`; Fortran:
`AsympInfo.f90`, `Analyzor.f90`). The on-the-fly analysis instead sums at the exact analysis energies
(above). The Fortran's interpolation has errors that the C++ fixes:
[fortran_fine_analysis_bugs.md](fortran_fine_analysis_bugs.md).

## Gas phase

The rovibrational states of the free molecule come from a Fourier grid Hamiltonian in r with V(r) at
z_asymptotic, diagonalised with LAPACK (`gasphase/`; Fortran: `FourierGridHamiltonian.f90`,
`GasPhase.f90`, where `jacobi.f90` diagonalises). They define the initial state and the vibrational
channels of the analysis.
