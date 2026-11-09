#pragma once
// The molecule-surface PES of a CRP data set: V of a diatomic molecule at a surface as a function of
// both atom positions, with the gradient (Corrugation Reducing Procedure: V = Va + Vb + I, the
// interpolation function I interpolated over the molecular cuts, switched to a 2D gas-phase
// potential far from the surface). Fortran: CRP_C3v/C4v/C6v_MoleculeSurface.f90 (POT6D, RRINTII,
// RRINT), CRP_GasPhase.f90.

#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "crplib/atom_surface.hpp"
#include "crplib/dataset.hpp"
#include "crplib/spline.hpp"

namespace crplib {

// V (Hartree) and its gradient with respect to both atom positions (Hartree/bohr).
struct MoleculeSurfaceValue {
  double v;
  Vec3 gradient_a;
  Vec3 gradient_b;
};

class MoleculeSurface {
 public:
  // Converts the data to internal units, reduces the cuts by the atom-surface PES (when the data
  // set holds total energies) and builds the splines. The options select bitwise-neutral speed-ups
  // (only the tests switch them off). Fortran: Initialize_CRP_*_MoleculeSurface_PES,
  // Read_MoleculeSurface_Data. Throws DataError.
  explicit MoleculeSurface(const DataSet& data, EvaluationOptions options = {});

  // V and gradient for atoms a and b (bohr, Cartesian; the surface atoms of the top layer at
  // z = 0). nullopt where the Fortran gives up silently (Fortran: switch = .TRUE., V = 0): the
  // molecule's (r, Z) outside a cut's grid below ZsoMax, or r outside the free-molecule potential
  // above it. Throws EvaluationError where the Fortran stops. Fortran: POT6D
  [[nodiscard]] std::optional<MoleculeSurfaceValue> evaluate(const Vec3& a, const Vec3& b) const;

  [[nodiscard]] Symmetry symmetry() const { return symmetry_; }
  [[nodiscard]] double lattice_constant() const { return lattice_constant_; }  // bohr
  // The (r, Z) range covered by all molecular cuts (the intersection). Fortran: Rmin, Rmax, Zmin
  [[nodiscard]] double r_min() const { return r_min_; }
  [[nodiscard]] double r_max() const { return r_max_; }
  [[nodiscard]] double z_min() const { return z_min_; }
  // The last Z of the gas-phase Z dependence: beyond it the PES no longer changes with Z.
  // Fortran: Zasy
  [[nodiscard]] double z_asymptotic() const { return z_dependence_.x().back(); }
  // The free-molecule potential at its last r (Hartree). Fortran: Rdep%F(Rdep%NX)
  [[nodiscard]] double h2_asymptote() const { return h2_.y().back(); }
  [[nodiscard]] const AtomSurface& atom_surface() const { return *atom_; }
  [[nodiscard]] std::vector<std::string> warnings() const;

  // Switching between the full and the gas-phase potential: F = 1 below ZsoMin, 0 above ZsoMax.
  // Fortran: Switch_Function
  struct Switch {
    double f, df_dz;
  };
  [[nodiscard]] Switch switch_function(double z) const;
  // The 2D gas-phase potential VH2(r) + Zdep(Z) and its derivatives; throws EvaluationError for
  // Z < ZsoMin (Fortran: STOP 'ERROR IN V2DINF'). Fortran: Calculate_CRP_GasPhase_Potential
  struct GasPhaseValue {
    double v, d_z, d_r;
  };
  [[nodiscard]] GasPhaseValue gas_phase(double z, double r) const;
  [[nodiscard]] double zso_min() const { return zso_min_; }
  [[nodiscard]] double zso_max() const { return zso_max_; }
  // Whether the molecule at (r, Zcc) would be evaluated (not the silent V = 0 branch).
  [[nodiscard]] bool covers(double r, double z) const;
  // The number of different (r, Z) grids among the molecular cuts used (S4: one lookup each).
  [[nodiscard]] int cut_grids() const { return static_cast<int>(grid_cut_.size()); }

 private:
  Symmetry symmetry_;
  double lattice_constant_;
  std::unique_ptr<const AtomSurface> atom_;
  double zso_min_, zso_max_;
  detail::CubicSpline1D z_dependence_;  // Fortran: GasPhase%Zdep
  detail::CubicSpline1D h2_;            // Fortran: GasPhase%Rdep
  std::vector<detail::BicubicSpline> cuts_;
  // S4: the used cuts grouped by bitwise identical (r, Z) grids. grid_cut_[g]: the first cut of
  // group g; group_of_cut_[k]: the group of cut k.
  std::vector<int> grid_cut_;
  std::vector<int> group_of_cut_;
  bool share_grid_lookup_;

  [[nodiscard]] bool evaluate_cuts(double r, double z, std::span<detail::BicubicValue> values) const;
  double r_min_, r_max_, z_min_;
  std::vector<std::string> warnings_;
};

}  // namespace crplib
