// Atom-surface PES of the C4v CRP data sets (fcc(100)). Fortran: CRP_C4v_AtomSurface.f90
// (Read_AtomSurface_Data, int_PES6CuHATOMCu100, MURO, VREP1D). The six site cuts fill a 5 x 5 grid
// over the unit cell by symmetry; a 3D spline periodic in x, y and clamped in z interpolates the
// reduced potential in the irreducible triangle, to which every point is folded (REGION).
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

#include "crplib/atom_surface_detail.hpp"
#include "crplib/errors.hpp"

namespace crplib::detail {

namespace {

// The 32 atoms of MURO in lattice units: 16 of the top layer and 16 of the second layer (z offset
// -1 = one interlayer spacing down). Fortran: VEC(i, 1..3)
struct ClusterAtom {
  double x, y, layer;
};
constexpr std::array<ClusterAtom, 32> muro_atoms = {{
    {0., -2., 0.},    {-1., -1., 0.},    {0., -1., 0.},     {1., -1., 0.},    {2., -1., 0.},     {-2., 0., 0.},
    {-1., 0., 0.},    {0., 0., 0.},      {1., 0., 0.},      {2., 0., 0.},     {-1., 1., 0.},     {0., 1., 0.},
    {1., 1., 0.},     {2., 1., 0.},      {0., 2., 0.},      {1., 2., 0.},     {-0.5, -1.5, -1.}, {0.5, -1.5, -1.},
    {1.5, -1.5, -1.}, {-1.5, -0.5, -1.}, {-0.5, -0.5, -1.}, {0.5, -0.5, -1.}, {1.5, -0.5, -1.},  {-1.5, 0.5, -1.},
    {-0.5, 0.5, -1.}, {0.5, 0.5, -1.},   {1.5, 0.5, -1.},   {2.5, 0.5, -1.},  {-0.5, 1.5, -1.},  {0.5, 1.5, -1.},
    {1.5, 1.5, -1.},  {0.5, 2.5, -1.},
}};

// Grid points (i, j) of the 5 x 5 spline grid that each site cut fills (1-based as in the Fortran;
// the sites in CRPData.3D order: top, hollow, bridge, t2h, b2h, t2b).
struct GridPoint {
  int i, j;
};
const std::array<std::vector<GridPoint>, 6> site_grid_points = {{
    {{1, 1}, {1, 5}, {5, 1}, {5, 5}},                                  // top
    {{3, 3}},                                                          // hollow
    {{3, 1}, {1, 3}, {5, 3}, {3, 5}},                                  // bridge
    {{2, 2}, {4, 4}, {4, 2}, {2, 4}},                                  // t2h
    {{3, 2}, {2, 3}, {4, 3}, {3, 4}},                                  // b2h
    {{2, 1}, {4, 1}, {1, 2}, {5, 2}, {1, 4}, {5, 4}, {2, 5}, {4, 5}},  // t2b
}};

}  // namespace

// Fortran: MURO and VREP1D
RepulsiveValue SquareAtomSurface::repulsive_sum(double x, double y, double z) const {
  const double a = geometry_.lattice_constant;
  const double spacing = geometry_.interlayer_spacings.front();
  RepulsiveValue sum{0.0, 0.0, 0.0, 0.0};
  for (const ClusterAtom& atom : muro_atoms) {
    const double xij = atom.x * a;
    const double yij = atom.y * a;
    // Only the first sub-layer is included. A negative first spacing is a relative value (in units
    // of a / sqrt(2)); it was unit-converted as a length all the same (Fortran quirk, kept).
    const double zij = spacing < 0.0 ? -atom.layer * spacing * a / std::sqrt(2.0) : atom.layer * spacing;
    const double r = std::sqrt((x - xij) * (x - xij) + (y - yij) * (y - yij) + (z - zij) * (z - zij));
    // VREP1D: zero beyond the corrugation function; no lower limit (SPLINT extrapolates).
    double vr = 0.0, dvdr = 0.0;
    if (!(r > corrugation_.x().back())) {
      const auto vs = corrugation_.value_and_slope_splint(r);
      vr = vs.value;
      dvdr = vs.slope;
    }
    sum.v = sum.v + vr;
    if (r > 0.0) {
      sum.dvdx = sum.dvdx + dvdr * (x - xij) / r;
      sum.dvdy = sum.dvdy + dvdr * (y - yij) / r;
      sum.dvdz = sum.dvdz + dvdr * (z - zij) / r;
    }
  }
  return sum;
}

// Fortran: Read_AtomSurface_Data. The reduced cuts (V minus MURO at the site) and their reduced end
// slopes fill the spline grid; x and y run over one unit cell (5 points, period a).
SquareAtomSurface::SquareAtomSurface(const DataSet& data)
    : geometry_(slab_geometry(data)), corrugation_(corrugation_function(data)) {
  const double a = geometry_.lattice_constant;
  const int nz = static_cast<int>(data.atom_cuts.front().potential.x.size());
  std::vector<double> z(nz);
  std::vector<double> reduced(static_cast<std::size_t>(6) * nz), first(6), last(6);
  for (int site = 0; site < 6; ++site) {
    const AtomCut& cut = data.atom_cuts[site];
    const double slope_first = convert_force(cut.slope_first, data.units, internal_units);
    const double slope_last = convert_force(cut.slope_last, data.units, internal_units);
    for (int i = 0; i < nz; ++i) {
      z[i] = convert_length(cut.potential.x[i], data.units.length, LengthUnit::bohr);
      const double pot = convert_energy(cut.potential.f[i], data.units.energy, EnergyUnit::hartree);
      const RepulsiveValue m = repulsive_sum(cut.x * a, cut.y * a, z[i]);
      reduced[site * nz + i] = pot - m.v;
      if (i == 0) first[site] = slope_first - m.dvdz;
      if (i == nz - 1) last[site] = slope_last - m.dvdz;
    }
  }
  // The 1D splines of the cuts (Compute_Spline1D in the Fortran) are not used by the evaluation.

  std::vector<double> grid(5);
  for (int i = 0; i < 5; ++i) grid[i] = static_cast<double>(i) * a * 0.25;  // Fortran: FLOAT(i-1) * a * 0.25
  std::vector<double> f(static_cast<std::size_t>(25) * nz), dfdz_first(25), dfdz_last(25);
  for (int site = 0; site < 6; ++site) {
    for (const GridPoint& p : site_grid_points[site]) {
      const int ij = (p.i - 1) + 5 * (p.j - 1);
      for (int k = 0; k < nz; ++k) f[ij + 25 * k] = reduced[site * nz + k];
      dfdz_first[ij] = first[site];
      dfdz_last[ij] = last[site];
    }
  }
  spline_ = PeriodicSpline3D(grid, grid, z, f, dfdz_first, dfdz_last);
}

// Fortran: int_PES6CuHATOMCu100
AtomSurfaceValue SquareAtomSurface::evaluate(const Vec3& p) const {
  if (p.z >= spline_.z().back()) return {0.0, {0.0, 0.0, 0.0}};
  const double a = geometry_.lattice_constant;
  auto [x0c, y0c] = fold_square(p.x, p.y, a);
  const int region = fold_into_triangle(x0c, y0c, a);
  const RepulsiveValue rep = repulsive_sum(x0c, y0c, p.z);
  const PeriodicSpline3D::Value s = spline_.evaluate(x0c, y0c, p.z);
  double dvdx = s.d_dx + rep.dvdx;
  double dvdy = s.d_dy + rep.dvdy;
  unfold_gradient(region, dvdx, dvdy);
  return {s.value + rep.v, {dvdx, dvdy, s.d_dz + rep.dvdz}};
}

// The distance from the folded triangle to the nearest atom of the first two layers that is not in
// the MURO cluster.
double SquareAtomSurface::coverage_radius() const {
  const double a = geometry_.lattice_constant;
  const std::array<std::pair<double, double>, 3> triangle = {{{0.0, 0.0}, {0.5 * a, 0.0}, {0.5 * a, 0.5 * a}}};
  double radius = std::numeric_limits<double>::infinity();
  for (int layer = 0; layer < 2; ++layer) {
    const double shift = layer == 0 ? 0.0 : 0.5;
    for (int i = -20; i <= 20; ++i) {
      for (int j = -20; j <= 20; ++j) {
        const double x = i + shift, y = j + shift;
        const bool in_cluster = std::any_of(muro_atoms.begin(), muro_atoms.end(), [&](const ClusterAtom& c) {
          return c.x == x && c.y == y && c.layer == -static_cast<double>(layer);
        });
        if (!in_cluster) radius = std::min(radius, distance_to_polygon(x * a, y * a, triangle));
      }
    }
  }
  return radius;
}

std::vector<std::string> SquareAtomSurface::warnings() const {
  std::vector<std::string> w;
  const double cutoff = corrugation_.x().back();
  if (cutoff > coverage_radius()) {
    w.push_back("the corrugation function reaches " + std::to_string(cutoff) + " bohr, beyond the " +
                std::to_string(coverage_radius()) +
                " bohr covered by the MURO atom cluster: the repulsive sum is not periodic (V jumps at cell edges)");
  }
  return w;
}

}  // namespace crplib::detail
