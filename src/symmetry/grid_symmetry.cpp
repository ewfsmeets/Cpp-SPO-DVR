// The point group on the (x, y, phi) grid: usable operations, orbits, representatives.
// Fortran: C2v/C4v/C6vSurfaceSymmetry.f90 (generalised)
#include "symmetry/grid_symmetry.hpp"

#include <numeric>
#include <utility>

#include "core/error.hpp"

namespace spodvr {

namespace {

int modulo(int a, int n) { return ((a % n) + n) % n; }

// (k nphi / 12) if k nphi is a multiple of 12, else -1.
int phi_shift(const SurfaceOperation& g, int nphi) {
  const int product = g.angle_twelfths * nphi;
  return product % 12 == 0 ? product / 12 : -1;
}

}  // namespace

GridSymmetry::GridSymmetry(int nx, int ny, int nphi) : nx_(nx), ny_(ny), nphi_(nphi) {
  operations_ = {identity_operation()};
  matrices_ = {IntegerMatrix2{1, 0, 0, 1}};
  build_orbits();
}

GridSymmetry::GridSymmetry(const SurfaceLattice& lattice, const PointGroup& group, int nx, int ny, int nphi)
    : group_name_(group.name()),
      group_order_(group.order()),
      phi_multiple_(spodvr::phi_points_multiple(group.operations())),
      nx_(nx),
      ny_(ny),
      nphi_(nphi) {
  if (nx == 1 || ny == 1)
    throw InputError("UseSymmetry", "surface symmetry needs both the x and the y degree of freedom (nX, nY > 1).");

  for (const SurfaceOperation& g : group.operations()) {
    const std::optional<IntegerMatrix2> a = lattice.lattice_matrix(g);
    if (!a)
      throw InputError("UseSymmetry", "point group " + group.name() + " does not fit the surface cell (" + describe(g) +
                                          " is not a lattice symmetry).");
    if (phi_shift(g, nphi) < 0) {
      dropped_.push_back({g, "nphi = " + std::to_string(nphi) + " is not a multiple of " +
                                 std::to_string(12 / std::gcd(12, g.angle_twelfths))});
      continue;
    }
    if ((a->a12 * nx) % ny != 0 || (a->a21 * ny) % nx != 0) {
      dropped_.push_back(
          {g, "it does not map the " + std::to_string(nx) + " x " + std::to_string(ny) + " (x, y) grid onto itself"});
      continue;
    }
    operations_.push_back(g);
    matrices_.push_back(*a);
  }
  build_orbits();
}

int GridSymmetry::apply_to_xy(int iop, int ix, int iy) const {
  const IntegerMatrix2& a = matrices_[iop];
  const int ix_image = modulo(a.a11 * ix + (a.a12 * nx_ / ny_) * iy, nx_);
  const int iy_image = modulo((a.a21 * ny_ / nx_) * ix + a.a22 * iy, ny_);
  return iy_image + ny_ * ix_image;
}

// Orbits in memory order: the first unassigned point starts a new orbit and is its
// representative; every image g rep of it gets the operation g^-1 (the first g in the
// canonical order that reaches it).
void GridSymmetry::build_orbits() {
  const int n_ops = static_cast<int>(operations_.size());
  representative_of_.assign(n_xy(), -1);
  operation_of_.assign(n_xy(), -1);
  std::vector<int> inverse_index(n_ops, -1);
  for (int i = 0; i < n_ops; ++i)
    for (int j = 0; j < n_ops; ++j)
      if (operations_[j] == inverse(operations_[i])) inverse_index[i] = j;

  for (int ix = 0; ix < nx_; ++ix)
    for (int iy = 0; iy < ny_; ++iy) {
      const int ixy = iy + ny_ * ix;
      if (representative_of_[ixy] >= 0) continue;
      const int irep = static_cast<int>(representatives_.size());
      representatives_.push_back(ixy);
      for (int iop = 0; iop < n_ops; ++iop) {
        const int image = apply_to_xy(iop, ix, iy);
        if (representative_of_[image] >= 0) continue;
        representative_of_[image] = irep;
        operation_of_[image] = inverse_index[iop];
      }
    }

  // phi_i -> s phi_i + 2 pi k/12 on the grid 2 pi i / nphi.
  phi_index_.clear();
  for (const SurfaceOperation& g : operations_) {
    const int shift = phi_shift(g, nphi_);
    std::vector<int> index(nphi_);
    for (int i = 0; i < nphi_; ++i) index[i] = modulo((g.mirror ? -i : i) + shift, nphi_);
    phi_index_.push_back(std::move(index));
  }
}

}  // namespace spodvr
