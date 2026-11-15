#pragma once
// The SCM set-up diagnostics, none of which the Fortran has: the slab against the CRP lattice (SCM
// D9), ΔV at the gas-phase potential (SCM D3), the jump of V at the cell edges (SCM D10) and ΔV where
// the CRP clamps Z or r (SCM D15), plus the slab's statistics. check_* returns a result as data for
// the set-up report, require_* throws an InputError naming the input key when it is unacceptable.

#include <optional>
#include <vector>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/crp_potential.hpp"
#include "potential/potential.hpp"
#include "potential/scm/scm_files.hpp"

namespace spodvr {

class ScmCrpPotential;

// Ideal atoms farther than this from the CRP lattice give a warning (bohr); the SRP48 example slab
// lies within 1.2e-6 bohr (SCM D9).
inline constexpr double scm_lattice_warning_distance = 1e-5;
// Ideal atoms farther than this from the CRP lattice are an error (bohr): ΔV would subtract the
// pair energies of another lattice than the one the CRP potential describes.
inline constexpr double scm_lattice_error_distance = 1e-3;
// |ΔV| on the plane Z = z_asymptotic above this gives a warning (eV): there the potential should be
// the gas phase's (SCM D3; the SRP48 example on the grid of crp_cu111_srp48_small from Z = 0.5 bohr:
// 4.841e-05 eV at Z = 11.9 bohr, r < 4.2 bohr, potential-only run 2026-11-14).
inline constexpr double scm_plane_warning_energy = 1e-3;
// A cell-edge point counts as accessible when V lies below min V_asy + EMaxWP + this (eV): the
// zero-point energy of H₂ (0.27 eV) plus a margin (SCM D10).
inline constexpr double scm_accessible_margin = 0.5;

// How far the slab's ideal atom positions lie from the sites of the CRP lattice (SCM D9). The
// deviation of an atom is its distance (bohr) from the nearest site of its layer.
struct ScmLatticeCheck {
  int top_layer_atoms = 0;              // layer-1 atoms checked
  int second_layer_atoms = 0;           // layer-2 atoms checked (none for a C6v data set)
  double top_layer_deviation = 0.0;     // bohr, the largest of layer 1
  double second_layer_deviation = 0.0;  // bohr, the largest of layer 2
  int worst_atom = -1;                  // index in ScmSlab::atoms of the largest deviation; −1: none checked

  // The largest deviation of both layers (NaN if one is NaN).
  [[nodiscard]] double max_deviation() const;
  // True if the largest deviation is above scm_lattice_warning_distance (or NaN). Out of line, so
  // that the NaN test is compiled with strict floating point in every build.
  [[nodiscard]] bool warning() const;
};

// SCM D9: the ideal layer 1 on the sites of the CRP lattice (a₀, γ of `crp`, unrescaled) at z = 0;
// with a C3v data set layer 2 at z = −d₁ under the data set's second-layer site (not checked for
// C6v, a one-layer PES).
[[nodiscard]] ScmLatticeCheck check_slab_lattice(const ScmSlab& slab, const SystemParameters& crp,
                                                 const CrpDataSetSummary& data_set);

// Throws InputError("PosFilepath", ...), naming the worst atom, if an ideal atom lies farther than
// scm_lattice_error_distance from its site (or the deviation is NaN).
void require_slab_lattice(const ScmLatticeCheck& check, const ScmSlab& slab, const SystemParameters& crp);

// The snapshot in numbers: atoms per layer and their displacements d = distorted − ideal (bohr).
struct ScmSlabStatistics {
  int top_layer_atoms = 0;
  int second_layer_atoms = 0;
  double rms_displacement = 0.0;  // √(mean |d|²) over the atoms used
  double max_displacement = 0.0;  // the largest |d|
  double highest_atom = 0.0;      // the largest distorted z: the atom nearest to the gas phase
  Vec3 top_layer_shift;           // the mean d of layer 1: the rigid part of its distortion
};

// Summed in file order, so the result does not depend on the build or the threads.
[[nodiscard]] ScmSlabStatistics slab_statistics(const ScmSlab& slab);

// A point of the 6D grid: lattice coordinates (u, v), Z and r (bohr), θ and φ (rad).
struct ScmPoint {
  double u = 0.0, v = 0.0, z = 0.0, r = 0.0, theta = 0.0, phi = 0.0;
};

// The extreme value of a scan and the first point (in scan order) where it occurs; a NaN wins and
// stays, so that it is reported instead of hidden.
struct ScmExtreme {
  double value = 0.0;    // eV
  ScmPoint at;
  long long points = 0;  // points scanned
};

// SCM D3: how ΔV = V − V_CRP(u/s, v/s) changes the gas-phase potential, at its points (0, 0, Z, r,
// π/2, 0) for every r of the grid, and on the plane Z = z_asymptotic. The points lie at z_asymptotic
// (V_asy) and, with UseSwitching, also at the PES's own zasymp, where the switch takes its V_gas.
struct ScmAsymptoticCheck {
  double z_asymptotic = 0.0;          // bohr (ZmaxSwitch with UseSwitching)
  std::vector<double> point_heights;  // bohr: z_asymptotic, then the switch's zasymp if it differs
  double gas_phase_minimum = 0.0;     // eV: the lowest V at the V_asy points (Z = z_asymptotic)
  int changed_points = 0;             // gas-phase points where ΔV ≠ 0 (or NaN): an error
  ScmExtreme largest_at_points;       // the largest |ΔV| at the gas-phase points
  // The largest |ΔV| on the plane over the grid's (x, y, θ, φ) and its r < ROptStart (absent in an
  // analysis-only run: no 6D potential).
  std::optional<ScmExtreme> largest_on_plane;
  double plane_r_limit = 0.0;       // ROptStart (bohr)

  // True if |ΔV| on the plane is above scm_plane_warning_energy (or NaN).
  [[nodiscard]] bool plane_warning() const;
};

// SCM D3 (data). With `with_plane` false only the V_asy points are checked. The plane scan runs over
// OpenMP; its result does not depend on the threads. `crp_scale` is s (1 without rescaling).
[[nodiscard]] ScmAsymptoticCheck check_asymptotic_potential(const Potential& scm, const Potential& crp,
                                                            double crp_scale, const Config& config,
                                                            const SpatialGrids& grids, const AngularGrid& angular,
                                                            bool with_plane);

// SCM D3: throws InputError("PosFilepath", ...) if ΔV changes the gas-phase potential at one of its
// points (or is NaN there): V_asy, the gas-phase states and the switch must not depend on the slab.
void require_unchanged_asymptotic_potential(const ScmAsymptoticCheck& check, const ScmSlab& slab,
                                            const ScmSlabStatistics& statistics);

// The jump of V across the cell edges at one (Z, r) and θ = 90°.
struct ScmEdgeJump {
  double z = 0.0;                    // bohr
  double r = 0.0;                    // bohr
  double max_jump = 0.0;             // eV, over every edge pair (NaN if a V is NaN)
  int accessible_pairs = 0;          // pairs with V below the accessible limit on both sides
  double max_accessible_jump = 0.0;  // eV, over those pairs (0 if there is none)
};

// SCM D10: |V(t·a, 0) − V(t·a, a)| and |V(0, t·a) − V(a, t·a)|, the two values the periodic grid joins,
// at the edge points t = k/n and the orientations φ of edge_jump_sampling, θ = 90°, for each (Z, r)
// of the sampling; a is the lattice constant of the run.
struct ScmEdgeJumps {
  double accessible_limit = 0.0;  // eV: min V_asy + EMaxWP + scm_accessible_margin
  int points_per_edge = 0;
  int phi_points = 0;             // φ = k·π/phi_points (H₂ is homonuclear: φ + π is the same geometry)
  std::vector<ScmEdgeJump> rows;  // per Z of edge_heights, each r of edge_bond_lengths (scm_diagnostics.cpp)
};

// SCM D10 (data), a few thousand evaluations of V (not on the grid; iz = ir = 0).
[[nodiscard]] ScmEdgeJumps check_edge_jumps(const Potential& scm, double accessible_limit);

// The CRP library's clamp margin (bohr): it evaluates Z < Zmin + δ at Zmin + δ, r < Rmin + δ at
// Rmin + δ and r > Rmax − δ at Rmax − δ (crplib qd_wrapper.cpp QdWrapper::clamped_r, clamped_z).
inline constexpr double crp_clamp_margin = 1e-12;

// True if the CRP library clamps this Z (or r) of the data set, with its own comparisons.
[[nodiscard]] bool crp_clamps_z(double z, const CrpDataSetSummary& data_set);
[[nodiscard]] bool crp_clamps_r(double r, const CrpDataSetSummary& data_set);

// SCM D15, the Z region: ΔV = V − V_CRP(u/s, v/s) on the Z planes of the 6D potential where the CRP
// clamps Z (below Zmin) and ΔV does not, at the r it does not clamp (Rmin .. Rmax). The r it clamps
// belong to ScmOutsideRRangeCheck, so the two scans cover every clamped point once.
struct ScmHoleCheck {
  double z_min = 0.0;   // the CRP data set's Zmin (bohr)
  double r_min = 0.0;   // the CRP data set's Rmin and Rmax (bohr)
  double r_max = 0.0;
  double z_start = 0.0; // the first Z point of the grid (bohr)
  int planes = 0;       // Z planes of the 6D potential below Zmin
  ScmExtreme lowest;    // the lowest ΔV over the grid points of these planes with Rmin ≤ r ≤ Rmax (planes > 0)

  // The Z grid extends below Zmin (SCM D15).
  [[nodiscard]] bool warning() const { return planes > 0; }
};

// SCM D15 (data). The scan runs over OpenMP; its result does not depend on the threads.
[[nodiscard]] ScmHoleCheck check_holes_below_z_min(const Potential& scm, const Potential& crp, double crp_scale,
                                                   const CrpDataSetSummary& data_set, const SpatialGrids& grids,
                                                   const AngularGrid& angular);

// SCM D15, the r region: ΔV = V − V_CRP(u/s, v/s) on every Z plane of the 6D potential at the r the CRP
// clamps (below Rmin or beyond Rmax), where ΔV uses the true r. Reported without a warning: grids
// usually reach beyond Rmax, because the r absorber lies at large r (SRP48: Rmax 4.35 bohr, the
// example's absorber 4.2–9.8 bohr); the report says whether the lowest point lies inside it.
struct ScmOutsideRRangeCheck {
  double r_min = 0.0;          // the CRP data set's Rmin and Rmax (bohr)
  double r_max = 0.0;
  int theta_points = 0;        // the θ of the angular grid scanned
  int phi_points = 0;          // the φ of the angular grid scanned
  ScmExtreme lowest;           // the lowest ΔV over the points scanned (none: no r outside the range)
  // The r absorber as applied (specular_operator.cpp, specular_kernels.cu; the 6D grid's factors are 1
  // outside it): the grid points of SpatialGrids::r_optical, none with nR = 1.
  int r_absorber_points = 0;
  double r_absorber_first = 0.0;  // bohr: its first and last grid point (r_absorber_points > 0)
  double r_absorber_last = 0.0;
  bool in_r_absorber = false;     // the lowest point is one of them
};

// SCM D15 (data), over every (x, y) and θ of the grid but only its first φ, to keep the cost small.
// The scan runs over OpenMP; its result does not depend on the threads.
[[nodiscard]] ScmOutsideRRangeCheck check_holes_outside_r_range(const Potential& scm, const Potential& crp,
                                                                double crp_scale, const CrpDataSetSummary& data_set,
                                                                const SpatialGrids& grids, const AngularGrid& angular);

// Everything the SCM set-up report shows besides the input and the slab.
struct ScmDiagnostics {
  double crp_lattice_constant = 0.0;        // a₀ (bohr)
  std::optional<double> rescaling_factor;   // s = a_T/a₀ with EnableSurfExpan
  ScmSlabStatistics slab;
  ScmLatticeCheck lattice;                  // SCM D9
  ScmAsymptoticCheck asymptotic;            // SCM D3
  std::optional<ScmEdgeJumps> edge_jumps;   // SCM D10; absent in an analysis-only run
  std::optional<ScmHoleCheck> holes;        // SCM D15; absent in an analysis-only run
  std::optional<ScmOutsideRRangeCheck> outside_r_range;  // SCM D15; absent in an analysis-only run
};

// All of the above for the run's potential and grids.
[[nodiscard]] ScmDiagnostics diagnose_scm_potential(const ScmCrpPotential& potential, const Config& config,
                                                    const SpatialGrids& grids, const AngularGrid& angular);

}  // namespace spodvr
