// The atom-surface PES: shared parts (unit conversion of the slab and the corrugation function,
// the unit-cell folds, polygon distances) and the choice of the implementation.
// Fortran: CRP_Wrapper.f90 SetUp_CRP_AtomSurface_PES, CRP_C*v_AtomSurface.f90 ORIGEN, ORIGEN2,
// REGION, TRANSF.
#include "crplib/atom_surface.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "crplib/atom_surface_detail.hpp"
#include "crplib/constants.hpp"
#include "crplib/errors.hpp"

namespace crplib {

std::unique_ptr<const AtomSurface> make_atom_surface(const DataSet& data, EvaluationOptions options) {
  switch (data.symmetry) {
    case Symmetry::c3v:
    case Symmetry::c6v: return std::make_unique<detail::HexagonalAtomSurface>(data, options);
    case Symmetry::c4v: return std::make_unique<detail::SquareAtomSurface>(data, options);
  }
  throw DataError("unknown symmetry");
}

}  // namespace crplib

namespace crplib::detail {

// Fortran: SetUp_CRP_AtomSurface_PES (LatticeConstant, InterlayerSpacings converted to bohr)
SlabGeometry slab_geometry(const DataSet& data) {
  SlabGeometry g;
  g.lattice_constant = convert_length(data.lattice_constant, data.units.length, LengthUnit::bohr);
  for (double spacing : data.interlayer_spacings)
    g.interlayer_spacings.push_back(convert_length(spacing, data.units.length, LengthUnit::bohr));
  return g;
}

// Fortran: Read_AtomSurface_Data, reading Atom/intrep.dat and Compute_Spline1D(..., 2, dz1, 2, dz2)
CubicSpline1D corrugation_function(const DataSet& data) {
  const double slope_first = convert_force(data.corrugation_slope_first, data.units, internal_units);
  const double slope_last = convert_force(data.corrugation_slope_last, data.units, internal_units);
  std::vector<double> x, f;
  for (std::size_t i = 0; i < data.corrugation.x.size(); ++i) {
    x.push_back(convert_length(data.corrugation.x[i], data.units.length, LengthUnit::bohr));
    f.push_back(convert_energy(data.corrugation.f[i], data.units.energy, EnergyUnit::hartree));
  }
  return CubicSpline1D(x, f, {EndCondition::first_derivative, slope_first},
                       {EndCondition::first_derivative, slope_last});
}

namespace {

// The Fortran's INTEGER = REAL assignment (truncation towards zero) of the cell index, with the
// extra -1 for negative values: floor, except one cell further for negative integers (quirk kept).
// Fortran: ORIGEN / ORIGEN2, ind(i) = ori(i) or ori(i) - 1
int cell_index(double ori) {
  constexpr double limit = 2.0e9;  // INT of a larger value is undefined behaviour in C++
  if (!(std::abs(ori) < limit)) {
    throw EvaluationError("CRP atom-surface PES: coordinate " + std::to_string(ori) +
                          " lattice constants from the origin");
  }
  return ori >= 0.0 ? static_cast<int>(ori) : static_cast<int>(ori - 1.0);
}

}  // namespace

// Fortran: ORIGEN2
std::pair<double, double> fold_hexagonal(double x, double y, double delta) {
  double bigx = x - y * usq3;
  double bigy = 2.00 * y * usq3;
  const int ind1 = cell_index(bigx / delta);
  const int ind2 = cell_index(bigy / delta);
  bigx = bigx - ind1 * delta;
  bigy = bigy - ind2 * delta;
  return {bigx + bigy * 0.500, sq3 * bigy * 0.500};
}

// Fortran: ORIGEN
std::pair<double, double> fold_square(double x, double y, double delta) {
  const int ind1 = cell_index(x / delta);
  const int ind2 = cell_index(y / delta);
  return {x - ind1 * delta, y - ind2 * delta};
}

// Fortran: REGION. The comparisons (.LE., .GE., .GT.) decide the boundaries and are kept as written.
int fold_into_triangle(double& x, double& y, double delta) {
  if (y <= delta / 2.0) {  // lower half of the unit cell
    if (x <= delta / 2.0) {
      if (x >= y) return 1;
      const double aux = x;
      x = y;
      y = aux;
      return 3;
    }
    if (y >= (delta - x)) {
      const double aux = x;
      x = y;
      y = delta - aux;
      return 4;
    }
    x = delta - x;
    return 2;
  }
  if (x <= delta / 2.0) {  // upper half
    if (y > (delta - x)) {
      y = delta - y;
      return 7;
    }
    const double aux = x;
    x = delta - y;
    y = aux;
    return 5;
  }
  if (x > y) {
    const double aux = x;
    x = delta - y;
    y = delta - aux;
    return 6;
  }
  x = delta - x;
  y = delta - y;
  return 8;
}

// Fortran: TRANSF
void unfold_gradient(int region, double& dvdx, double& dvdy) {
  double aux = 0.0;
  switch (region) {
    case 1: return;
    case 2: dvdx = -dvdx; return;
    case 3:
      aux = dvdx;
      dvdx = dvdy;
      dvdy = aux;
      return;
    case 4:
      aux = dvdx;
      dvdx = -dvdy;
      dvdy = aux;
      return;
    case 5:
      aux = dvdx;
      dvdx = dvdy;
      dvdy = -aux;
      return;
    case 6:
      aux = dvdx;
      dvdx = -dvdy;
      dvdy = -aux;
      return;
    case 7: dvdy = -dvdy; return;
    case 8:
      dvdx = -dvdx;
      dvdy = -dvdy;
      return;
    default: throw EvaluationError("CRP C4v atom-surface PES: unknown region " + std::to_string(region));
  }
}

double distance_to_polygon(double x, double y, std::span<const std::pair<double, double>> corners) {
  const std::size_t n = corners.size();
  bool inside = true;
  double distance = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; ++i) {
    const auto [ax, ay] = corners[i];
    const auto [bx, by] = corners[(i + 1) % n];
    const double ex = bx - ax, ey = by - ay;
    if (ex * (y - ay) - ey * (x - ax) < 0.0) inside = false;  // right of a counter-clockwise edge
    const double t = std::clamp(((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey), 0.0, 1.0);
    distance = std::min(distance, std::hypot(x - (ax + t * ex), y - (ay + t * ey)));
  }
  return inside ? 0.0 : distance;
}

// S2. In-plane, an atom is kept if a point within 2 cell_margin of the cell can be within reach of
// it (the folded points are checked against cell_margin / 2 per call). Every computed distance r
// satisfies r >= (1 - 3 eps) times the in-plane or the z distance, so an atom left out has a
// computed r > cutoff and contributes zero in the full sum too.
ReachableCluster::ReachableCluster(std::span<const ClusterAtom> cluster,
                                   std::span<const std::pair<double, double>> cell, double cutoff,
                                   double lattice_constant)
    : reach_(cutoff > 0.0 ? cutoff * (1.0 + 1.0e-10) : std::numeric_limits<double>::infinity()),
      cell_margin_(1.0e-8 * lattice_constant) {
  for (const ClusterAtom& atom : cluster)
    if (distance_to_polygon(atom.x, atom.y, cell) <= reach_ + 2.0 * cell_margin_) atoms_.push_back(atom);
}

}  // namespace crplib::detail
