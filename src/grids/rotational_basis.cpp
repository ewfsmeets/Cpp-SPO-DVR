// Rotational basis |J mJ> and the index of the initial state.
// Fortran: transform.f90 InitGenComplexTransform, RotDifBasis.f90 RotationBasis
#include "grids/rotational_basis.hpp"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>

namespace spodvr {

RotationalBasis::RotationalBasis(const AngularGrid& grid, int j_in, int mj_in)
    : mj_max_(std::min(grid.mj_max(), grid.j_max())), initial_(-1) {
  for (int m = -mj_max_; m <= mj_max_; ++m) {
    const int first = size();
    for (int j = grid.j_min(); j <= grid.j_max(); j += grid.delta_j()) {
      if (std::abs(m) > j) continue;
      if (j == j_in && m == mj_in) initial_ = size();
      j_.push_back(j);
      mj_.push_back(m);
    }
    // Every mJ has a block: J = JMax >= MjMax is always in the basis.
    blocks_.push_back({first, size() - 1});
  }
  if (initial_ < 0) throw std::invalid_argument("RotationalBasis: the initial state (JIn, MJIn) is not in the basis");
}

}  // namespace spodvr
