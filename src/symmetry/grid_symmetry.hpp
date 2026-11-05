#pragma once

#include <string>
#include <vector>

#include "symmetry/point_group.hpp"
#include "symmetry/surface_lattice.hpp"

namespace spodvr {

// The surface point group acting on the (x, y, phi) grid, and the reduction of the potential it
// allows. Built once at setup.
//
// An operation g is usable on the grid if it maps grid points onto grid points:
//   (x, y): its lattice matrix A sends (ix, iy) -> (a11 ix + a12 (nx/ny) iy, a21 (ny/nx) ix + a22 iy)
//           mod (nx, ny), which needs a12 nx and a21 ny to be multiples of ny and nx;
//   phi:    phi_i = 2 pi i / nphi -> s phi_i + 2 pi k/12, which needs k nphi to be a multiple of 12.
// The usable operations form a subgroup (Fortran: triangle A instead of triangle 1 when nphi is
// not a multiple of 6). The (x, y) points fall into orbits; the first point of each orbit in
// memory order (ixy = iy + ny ix) is its representative. For every point p the operation h with
// h p = rep(p) is stored; since V(h q) = V(q),
//   V(p, theta, phi_i) = V(rep(p), theta, phi_{phi_index(h)[i]}),
// so the potential is needed at the representatives only. theta is never changed.
// Fortran: C2v/C4v/C6vSurfaceSymmetry.f90 setupReducedXYGrids*, setup*SymmetryIndices
// (generalised: the Fortran stores one index per (phi, theta, x, y) point instead)
class GridSymmetry {
 public:
  // An operation of the group that cannot be used on this grid, and why.
  struct DroppedOperation {
    SurfaceOperation operation;
    std::string reason;
  };

  // Throws InputError (parameter "UseSymmetry") if an operation of the group is not a symmetry of
  // the lattice (wrong point group for the cell) or if nx or ny is 1.
  GridSymmetry(const SurfaceLattice& lattice, const PointGroup& group, int nx, int ny, int nphi);

  // No symmetry: every (x, y) point is its own representative.
  GridSymmetry(int nx, int ny, int nphi);

  [[nodiscard]] const std::string& point_group() const { return group_name_; }
  [[nodiscard]] int group_order() const { return group_order_; }
  [[nodiscard]] const std::vector<SurfaceOperation>& operations() const { return operations_; }  // usable
  [[nodiscard]] const std::vector<DroppedOperation>& dropped_operations() const { return dropped_; }
  // nphi must be a multiple of this for all operations of the group to be usable.
  [[nodiscard]] int phi_points_multiple() const { return phi_multiple_; }

  [[nodiscard]] int nx() const { return nx_; }
  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] int nphi() const { return nphi_; }
  [[nodiscard]] int n_xy() const { return nx_ * ny_; }
  [[nodiscard]] int n_representatives() const { return static_cast<int>(representatives_.size()); }
  // Full-grid ixy of representative irep.
  [[nodiscard]] int representative_xy(int irep) const { return representatives_[irep]; }
  // Representative (0 .. n_representatives - 1) of the (x, y) point ixy, and the index (into
  // operations()) of the operation h with h ixy = representative.
  [[nodiscard]] int representative_of(int ixy) const { return representative_of_[ixy]; }
  [[nodiscard]] int operation_of(int ixy) const { return operation_of_[ixy]; }
  // phi index of the image of phi_i under operations()[iop].
  [[nodiscard]] const std::vector<int>& phi_index(int iop) const { return phi_index_[iop]; }

  // (x, y) grid index of the image of (ix, iy) under operations()[iop].
  [[nodiscard]] int apply_to_xy(int iop, int ix, int iy) const;

  // True if the reduction is trivial (only the identity is usable).
  [[nodiscard]] bool is_trivial() const { return operations_.size() == 1; }

 private:
  void build_orbits();

  std::string group_name_ = "none";
  int group_order_ = 1;
  int phi_multiple_ = 1;
  int nx_;
  int ny_;
  int nphi_;
  std::vector<SurfaceOperation> operations_;
  std::vector<IntegerMatrix2> matrices_;  // lattice matrix of each usable operation
  std::vector<DroppedOperation> dropped_;
  std::vector<int> representatives_;
  std::vector<int> representative_of_;
  std::vector<int> operation_of_;
  std::vector<std::vector<int>> phi_index_;
};

}  // namespace spodvr
