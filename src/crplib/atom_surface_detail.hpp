#pragma once
// The atom-surface implementations of the CRP library, visible to their tests. Fortran:
// CRP_C6v_AtomSurface.f90 and CRP_C3v_AtomSurface.f90 (hexagonal: ATII, FOUR3II, MURII, ORIGEN2),
// CRP_C4v_AtomSurface.f90 (square: int_PES6CuHATOMCu100, ORIGEN, REGION, TRANSF, MURO, VREP1D).

#include <array>
#include <span>
#include <utility>
#include <vector>

#include "crplib/atom_surface.hpp"
#include "crplib/periodic_spline_3d.hpp"
#include "crplib/spline.hpp"

namespace crplib::detail {

// The repulsive sum and its gradient. Fortran: MURII / MURO outputs v, dvdx, dvdy, dvdz
struct RepulsiveValue {
  double v, dvdx, dvdy, dvdz;
};

// The data set's slab geometry in internal units. Fortran: TYPE(Slab_Geometry) after
// Convert_Length_Units
struct SlabGeometry {
  double lattice_constant;
  std::vector<double> interlayer_spacings;
};
SlabGeometry slab_geometry(const DataSet& data);

// The corrugation function (Atom/intrep.dat) in internal units, clamped at both ends with the given
// slopes. Fortran: Read_AtomSurface_Data, CorrugationFunction1D
CubicSpline1D corrugation_function(const DataSet& data);

// Shifts (x, y) into the unit cell with a corner at the origin; Cartesian in and out. delta is the
// lattice constant. Fortran: ORIGEN2 (hexagonal, skewed cell)
std::pair<double, double> fold_hexagonal(double x, double y, double delta);
// Fortran: ORIGEN (square cell)
std::pair<double, double> fold_square(double x, double y, double delta);
// Maps a point of the square unit cell into the triangle 0 <= y <= x <= delta/2 by one of the 8
// symmetry operations; returns the region number 1..8 (needed to map the gradient back).
// Fortran: REGION
int fold_into_triangle(double& x, double& y, double delta);
// The gradient of region `region` from the gradient in the triangle. Fortran: TRANSF
void unfold_gradient(int region, double& dvdx, double& dvdy);

// Distance from the point (x, y) (Cartesian, bohr) to the closed polygon `corners` (convex,
// counter-clockwise).
double distance_to_polygon(double x, double y, std::span<const std::pair<double, double>> corners);

// One atom of the MURII / MURO cluster at its position in internal units, computed as the Fortran
// loops compute it.
struct ClusterAtom {
  double x, y, z;
  int layer;  // 0: top layer
};

// Speed-up S2 (docs/crp_library.md "Speed"): the cluster atoms that can reach a point of the fold
// cell. An atom farther than the last point of the corrugation function adds exactly +0 to the
// repulsive sum and +-0 to its gradient; the sums start at +0 and can never become -0, so leaving
// such an atom out changes nothing, bitwise. Both tests below are conservative: a computed
// distance at or within rounding of the cut-off is never skipped.
class ReachableCluster {
 public:
  ReachableCluster() = default;
  // The atoms (in cluster order) whose in-plane distance to the fold cell `cell` (convex,
  // counter-clockwise) is at most reach() + 2 cell_margin().
  ReachableCluster(std::span<const ClusterAtom> cluster, std::span<const std::pair<double, double>> cell, double cutoff,
                   double lattice_constant);

  [[nodiscard]] const std::vector<ClusterAtom>& atoms() const { return atoms_; }
  // The cut-off with a relative margin of 1e-10 (infinite for a cut-off <= 0): an atom with
  // |z - z_atom| beyond it is out of reach (the per-call layer test).
  [[nodiscard]] double reach() const { return reach_; }
  // 1e-8 lattice constants: the folded point is used with the reachable atoms only if it lies
  // within cell_margin() / 2 of the cell (else, e.g. after the rounding of a fold far from the
  // origin, the full cluster is used).
  [[nodiscard]] double cell_margin() const { return cell_margin_; }

 private:
  std::vector<ClusterAtom> atoms_;
  double reach_ = 0.0;
  double cell_margin_ = 0.0;
};

// C3v and C6v. Fortran: ATII with FOUR3II (7 or 10 basis functions), MURII and ORIGEN2
class HexagonalAtomSurface final : public AtomSurface {
 public:
  explicit HexagonalAtomSurface(const DataSet& data, EvaluationOptions options = {});
  [[nodiscard]] AtomSurfaceValue evaluate(const Vec3& position) const override;
  [[nodiscard]] std::vector<std::string> warnings() const override;
  [[nodiscard]] ClusterSize cluster_size() const override;

  // The sum of the corrugation function over the surface atoms near (x, y, z): C6v 30 atoms of the
  // top layer, C3v an 8 x 8 x 4 cluster of four ABC-stacked layers. Fortran: MURII
  [[nodiscard]] RepulsiveValue repulsive_sum(double x, double y, double z) const;
  // The same sum over the atoms within reach (S2), for (x, y) folded into the unit cell; bitwise
  // equal to repulsive_sum.
  [[nodiscard]] RepulsiveValue reachable_repulsive_sum(double x, double y, double z) const;
  [[nodiscard]] const std::vector<ClusterAtom>& cluster() const { return cluster_; }
  [[nodiscard]] const ReachableCluster& reachable_cluster() const { return reachable_; }
  // The Fourier basis functions and their x, y derivatives at (x, y). Fortran: FOUR3II
  void fourier_basis(double x, double y, std::span<double> f, std::span<double> dfdx, std::span<double> dfdy) const;
  // Radius around the folded unit cell within which every surface atom is in the cluster.
  [[nodiscard]] double coverage_radius() const;
  [[nodiscard]] double lattice_constant() const { return geometry_.lattice_constant; }

 private:
  [[nodiscard]] RepulsiveValue sum_over(std::span<const ClusterAtom> atoms, double x, double y, double z,
                                        double z_reach) const;
  void add_atom(double dx, double dy, double dz, int layer, RepulsiveValue& sum) const;
  [[nodiscard]] bool in_fold_cell(double x, double y) const;

  Symmetry symmetry_;
  SlabGeometry geometry_;
  CubicSpline1D corrugation_;
  std::vector<ClusterAtom> cluster_;  // MURII, in the Fortran loop order
  ReachableCluster reachable_;        // S2
  bool prune_cluster_;
  std::vector<CubicSpline1D> cuts_;  // reduced: data minus repulsive sum
  std::span<const double> table_;    // fouri3D (n x n)
  int n_basis_;
};

// C4v. Fortran: int_PES6CuHATOMCu100 with ORIGEN, REGION, MURO, the periodic 3D spline, TRANSF
class SquareAtomSurface final : public AtomSurface {
 public:
  explicit SquareAtomSurface(const DataSet& data, EvaluationOptions options = {});
  [[nodiscard]] AtomSurfaceValue evaluate(const Vec3& position) const override;
  [[nodiscard]] std::vector<std::string> warnings() const override;
  [[nodiscard]] ClusterSize cluster_size() const override;

  // The sum of the corrugation function over 32 atoms of the first two layers. Fortran: MURO
  [[nodiscard]] RepulsiveValue repulsive_sum(double x, double y, double z) const;
  // The same sum over the atoms within reach (S2), for (x, y) folded into the triangle of REGION;
  // bitwise equal to repulsive_sum.
  [[nodiscard]] RepulsiveValue reachable_repulsive_sum(double x, double y, double z) const;
  [[nodiscard]] const std::vector<ClusterAtom>& cluster() const { return cluster_; }
  [[nodiscard]] const ReachableCluster& reachable_cluster() const { return reachable_; }
  [[nodiscard]] double coverage_radius() const;
  [[nodiscard]] const PeriodicSpline3D& spline() const { return spline_; }

 private:
  [[nodiscard]] RepulsiveValue sum_over(std::span<const ClusterAtom> atoms, double x, double y, double z,
                                        double z_reach) const;
  void add_atom(double x, double y, double z, const ClusterAtom& atom, RepulsiveValue& sum) const;
  [[nodiscard]] bool in_fold_triangle(double x, double y) const;

  SlabGeometry geometry_;
  CubicSpline1D corrugation_;
  std::vector<ClusterAtom> cluster_;  // MURO, in the Fortran order
  ReachableCluster reachable_;        // S2
  bool prune_cluster_;
  PeriodicSpline3D spline_;
};

}  // namespace crplib::detail
