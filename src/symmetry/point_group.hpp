#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "symmetry/surface_operation.hpp"

namespace spodvr {

// A surface point group (about the origin of the PES, x = y = 0) as the list of its operations:
// identity first, then the rotations by increasing angle, then the mirrors by increasing angle.
// Built from generators by closure:
//   C2v  rotation 180 deg, mirror at 0 deg            (rectangular cell)
//   C4v  rotation  90 deg, mirror at 0 deg            (square cell)
//   C6v  rotation  60 deg, mirror at 0 deg            (hexagonal cell)
//   C3v  rotation 120 deg, mirror at 90 deg           (hexagonal cell, top site at the origin: the
//                                                       mirrors through the top site that keep the
//                                                       fcc and hcp hollows in place; p3m1)
//   C3   rotation 120 deg
//   none the identity only
// Names are compared like the Fortran (.SEQ.: case-insensitive, blanks ignored at both ends).
// Throws InputError (parameter "SystemID") for an unknown name.
// Fortran: CxvSurfaceSymmetry.f90 (only C6v is used there; C2v, C4v are switched off and C3v, C3
// have no implementation)
class PointGroup {
 public:
  explicit PointGroup(std::string_view name);

  [[nodiscard]] const std::string& name() const { return name_; }
  [[nodiscard]] const std::vector<SurfaceOperation>& operations() const { return operations_; }
  [[nodiscard]] const std::vector<SurfaceOperation>& generators() const { return generators_; }
  [[nodiscard]] int order() const { return static_cast<int>(operations_.size()); }

 private:
  std::string name_;
  std::vector<SurfaceOperation> generators_;
  std::vector<SurfaceOperation> operations_;
};

// The closure of a set of operations under composition, in the canonical order (identity,
// rotations by increasing angle, mirrors by increasing angle).
[[nodiscard]] std::vector<SurfaceOperation> group_closure(const std::vector<SurfaceOperation>& generators);

// True if the operations form a group (contain the identity, closed under composition).
[[nodiscard]] bool is_group(const std::vector<SurfaceOperation>& operations);

// The smallest number of phi points m for which every operation maps the phi grid
// 2 pi i / nphi onto itself whenever nphi is a multiple of m: 12 / gcd(12, all angles).
// Fortran: getCxvNumPhisSmallestMultiple (2, 4, 6 for C2v, C4v, C6v)
[[nodiscard]] int phi_points_multiple(const std::vector<SurfaceOperation>& operations);

}  // namespace spodvr
