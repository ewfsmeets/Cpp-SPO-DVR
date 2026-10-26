#pragma once

#include <vector>

#include "core/index_range.hpp"
#include "grids/angular_grid.hpp"

namespace spodvr {

// The rotational (FBR) basis |J mJ> of the scattering wave function. Ordering as in Fortran:
// mJ = -MjMax .. MjMax (outer), J = JMin .. JMax in steps of dJ (inner), only |mJ| <= J.
// All indices are 0-based. Fortran: transform.f90 InitGenComplexTransform (JBasisSet,
// MjBasisSet, TheTransformIndices), RotDifBasis.f90 RotationBasis (nrotw, ichn)
class RotationalBasis {
 public:
  RotationalBasis(const AngularGrid& grid, int j_in, int mj_in);

  [[nodiscard]] int size() const { return static_cast<int>(j_.size()); }  // Fortran: nrotw
  [[nodiscard]] int j(int irot) const { return j_[irot]; }
  [[nodiscard]] int mj(int irot) const { return mj_[irot]; }
  [[nodiscard]] int mj_max() const { return mj_max_; }     // Fortran: TheTransform%MjMax
  [[nodiscard]] int initial_index() const { return initial_; }  // (JIn, MJIn). Fortran: ichn - 1
  // Indices of the basis functions with the given mJ (consecutive). Fortran: TheTransformIndices
  [[nodiscard]] IndexRange mj_block(int mj) const { return blocks_[mj + mj_max_]; }

 private:
  int mj_max_;
  int initial_;
  std::vector<int> j_;
  std::vector<int> mj_;
  std::vector<IndexRange> blocks_;  // per mJ = -MjMax .. MjMax
};

}  // namespace spodvr
