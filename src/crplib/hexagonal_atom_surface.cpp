// Atom-surface PES of the C3v and C6v CRP data sets. Fortran: CRP_C6v_AtomSurface.f90 and
// CRP_C3v_AtomSurface.f90 (Read_AtomSurface_Data, ATII, FOUR3II, MURII, ORIGEN2). The two differ in
// the Fourier basis (C3v adds the sine combinations) and in MURII (C6v: 30 atoms of the top layer;
// C3v: an 8 x 8 x 4 cluster of four ABC-stacked layers, E.W.F. Smeets 2018).
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

#include "crplib/atom_surface_detail.hpp"
#include "crplib/constants.hpp"
#include "crplib/errors.hpp"
#include "crplib/fourier_tables.hpp"

namespace crplib::detail {

namespace {

// Crystal coordinates (i, j) of the 30 top-layer atoms of the C6v MURII. Fortran: DATA ivec
constexpr std::array<int, 30> c6v_cluster_i = {0,  1,  2, 3, -1, 0, 1,  2,  3, -2, -1, 0,  1,  2, 3,
                                               -2, -1, 0, 1, 2,  3, -2, -1, 0, 1,  2,  -2, -1, 0, 1};
constexpr std::array<int, 30> c6v_cluster_j = {-2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0,
                                               1,  1,  1,  1,  1,  1,  2,  2,  2,  2, 2, 3, 3, 3, 3};

// The C3v MURII cluster: i, j = -n .. n+1 in every layer k = 0 .. m. Fortran: PARAMETER(n=3, m=3)
constexpr int c3v_cluster_n = 3;
constexpr int c3v_cluster_layers = 4;

// The atom positions of MURII, in its loop order, computed as the Fortran loops compute them.
std::vector<ClusterAtom> murii_cluster(Symmetry symmetry, const SlabGeometry& geometry) {
  const double a = geometry.lattice_constant;
  std::vector<ClusterAtom> cluster;
  if (symmetry == Symmetry::c6v) {  // CRP_C6v_AtomSurface.f90 MURII
    for (std::size_t i = 0; i < c6v_cluster_i.size(); ++i) {
      const double bigxij = c6v_cluster_i[i] * a;
      const double bigyij = c6v_cluster_j[i] * a;
      const double yij = bigyij * sq3 * 0.500;
      const double xij = bigxij + bigyij * 0.500;
      cluster.push_back({xij, yij, 0.0, 0});
    }
    return cluster;
  }
  // CRP_C3v_AtomSurface.f90 MURII (loop order i, j, k)
  for (int i = -c3v_cluster_n; i <= c3v_cluster_n + 1; ++i) {
    for (int j = -c3v_cluster_n; j <= c3v_cluster_n + 1; ++j) {
      for (int k = 0; k < c3v_cluster_layers; ++k) {
        // fcc(111) ABC stacking: layer k shifted by (k mod 3)/3 along both lattice vectors
        const double uuc = (static_cast<double>(i) + (k % 3) * (1.0 / 3.0)) * a;
        const double vuc = (static_cast<double>(j) + (k % 3) * (1.0 / 3.0)) * a;
        const double xijk = uuc + vuc * cos_pi_3;
        const double yijk = vuc * sin_pi_3;
        double zijk = 0.0;
        for (int l = 0; l < k; ++l) zijk = zijk - geometry.interlayer_spacings[l];
        cluster.push_back({xijk, yijk, zijk, k});
      }
    }
  }
  return cluster;
}

// The STOPs of MURII, out of line (keeps add_atom small enough to be inlined into the cluster loops).
[[noreturn]] void throw_atom_too_close(double r) {
  throw EvaluationError("CRP atom-surface PES: atom " + std::to_string(r) +
                        " bohr from a surface atom, closer than the corrugation function (Fortran MURII: R too small)");
}
[[noreturn]] void throw_atom_on_subsurface_atom() {
  throw EvaluationError("CRP atom-surface PES: atom on a sub-surface atom (Fortran MURII: r smaller then 0)");
}

// The unit cell into which ORIGEN2 folds (closed rhombus, counter-clockwise), Cartesian.
std::array<std::pair<double, double>, 4> hexagonal_cell(double a) {
  const double h = a * sq3 * 0.5;
  return {{{0.0, 0.0}, {a, 0.0}, {1.5 * a, h}, {0.5 * a, h}}};
}

}  // namespace

HexagonalAtomSurface::HexagonalAtomSurface(const DataSet& data, EvaluationOptions options)
    : symmetry_(data.symmetry),
      geometry_(slab_geometry(data)),
      corrugation_(corrugation_function(data)),
      cluster_(murii_cluster(data.symmetry, geometry_)),
      reachable_(cluster_, hexagonal_cell(geometry_.lattice_constant), corrugation_.x().back(),
                 geometry_.lattice_constant),
      prune_cluster_(options.prune_cluster),
      table_(data.symmetry == Symmetry::c3v ? std::span<const double>(fourier_table_atom_c3v())
                                            : std::span<const double>(fourier_table_atom_c6v())),
      n_basis_(data.symmetry == Symmetry::c3v ? 10 : 7) {
  // Fortran: Read_AtomSurface_Data. Each cut is stored reduced: V minus the repulsive sum at the
  // site; the end slopes are reduced in the same way.
  const double a = geometry_.lattice_constant;
  for (const AtomCut& cut : data.atom_cuts) {
    const double slope_first = convert_force(cut.slope_first, data.units, internal_units);
    const double slope_last = convert_force(cut.slope_last, data.units, internal_units);
    const std::size_t n = cut.potential.x.size();
    std::vector<double> z(n), reduced(n);
    double der1 = 0.0, der2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      z[i] = convert_length(cut.potential.x[i], data.units.length, LengthUnit::bohr);
      const double pot = convert_energy(cut.potential.f[i], data.units.energy, EnergyUnit::hartree);
      RepulsiveValue m{};
      try {
        m = repulsive_sum(cut.x * a, cut.y * a, z[i]);
      } catch (const EvaluationError& e) {
        throw DataError("Atom/" + cut.file_name + ": " + e.what());
      }
      reduced[i] = pot - m.v;
      if (i == 0) der1 = slope_first - m.dvdz;
      if (i == n - 1) der2 = slope_last - m.dvdz;
    }
    cuts_.emplace_back(z, reduced, SplineEnd{EndCondition::first_derivative, der1},
                       SplineEnd{EndCondition::first_derivative, der2});
  }
}

// One atom's contribution to MURII (the corrugation function at its distance). layer > 0 marks an
// atom below the top layer (an atom at distance 0 there is an error).
void HexagonalAtomSurface::add_atom(double dx, double dy, double dz, int layer, RepulsiveValue& sum, Want want) const {
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  double vr = 0.0, dvdr = 0.0;
  if (!(r > corrugation_.x().back())) {
    if (r < corrugation_.x().front()) throw_atom_too_close(r);
    if (want == Want::value) {
      vr = corrugation_.value_splint(r);
    } else {
      const auto vs = corrugation_.value_and_slope_splint(r);
      vr = vs.value;
      dvdr = vs.slope;
    }
  }
  sum.v = sum.v + vr;
  if (!(r > 0.0) && layer != 0) throw_atom_on_subsurface_atom();
  if (want == Want::value) return;
  if (r > 0.0) {
    const double dvdrsr = dvdr / r;
    sum.dvdx = sum.dvdx + dvdrsr * dx;
    sum.dvdy = sum.dvdy + dvdrsr * dy;
    sum.dvdz = sum.dvdz + dvdrsr * dz;
  } else {  // an atom of the top layer at the surface atom: only the z derivative
    sum.dvdz = sum.dvdz + dvdr;
  }
}

// Fortran: MURII (C6v: the 30 atoms of the top layer; C3v: the 8 x 8 x 4 cluster)
RepulsiveValue HexagonalAtomSurface::repulsive_sum(double x, double y, double z, Want want) const {
  return sum_over(cluster_, x, y, z, std::numeric_limits<double>::infinity(), want);
}

// S2: the atoms out of reach left out (in the plane once, in z per call); the rest in cluster order.
RepulsiveValue HexagonalAtomSurface::reachable_repulsive_sum(double x, double y, double z, Want want) const {
  if (!in_fold_cell(x, y)) return repulsive_sum(x, y, z, want);
  return sum_over(reachable_.atoms(), x, y, z, reachable_.reach(), want);
}

// The sum over `atoms` in their order, skipping those with |z - z_atom| > z_reach (S2).
RepulsiveValue HexagonalAtomSurface::sum_over(std::span<const ClusterAtom> atoms, double x, double y, double z,
                                              double z_reach, Want want) const {
  RepulsiveValue sum{0.0, 0.0, 0.0, 0.0};
  for (const ClusterAtom& atom : atoms) {
    if (std::abs(z - atom.z) > z_reach) continue;
    add_atom(x - atom.x, y - atom.y, z - atom.z, atom.layer, sum, want);
  }
  return sum;
}

// (x, y) within half the cell margin of the unit cell, in skewed coordinates (bohr).
bool HexagonalAtomSurface::in_fold_cell(double x, double y) const {
  const double a = geometry_.lattice_constant;
  const double limit = 0.5 * reachable_.cell_margin();
  const double u = x - y * usq3;
  const double v = 2.0 * y * usq3;
  return u >= -limit && u <= a + limit && v >= -limit && v <= a + limit;
}

ClusterSize HexagonalAtomSurface::cluster_size() const {
  return {static_cast<int>(cluster_.size()), static_cast<int>(reachable_.atoms().size())};
}

// Fortran: FOUR3II. Basis functions of the hexagonal cell (sums over stars of reciprocal lattice
// vectors) in crystal coordinates; the y derivative is transformed to Cartesian at the end.
void HexagonalAtomSurface::fourier_basis(double xa, double ya, std::span<double> fxy, std::span<double> dfdx,
                                         std::span<double> dfdy, Want want) const {
  const bool c3v = symmetry_ == Symmetry::c3v;
  const double fac = 2.00 * pi / geometry_.lattice_constant;
  const double bigx = xa - ya * usq3;
  const double bigy = 2.00 * ya * usq3;
  const double fx = fac * bigx;  // Fortran: fbigx
  const double fy = fac * bigy;  // Fortran: fbigy

  // The sines enter the C3v sine stars (2, 4, 5) and the derivatives; Want::value (S3) skips the
  // others (then dfdx, dfdy are meaningless).
  const bool gradient = want == Want::value_and_gradient;
  const bool star_sines = c3v || gradient;
  auto sine = [](bool needed, double angle) { return needed ? std::sin(angle) : 0.0; };
  const double c2a = std::cos(fx), c2b = std::cos(fy), c2c = std::cos(fx + fy);
  const double s2a = sine(star_sines, fx), s2b = sine(star_sines, fy), s2c = sine(star_sines, fx + fy);
  const double c3a = std::cos(2.00 * fx + fy), c3b = std::cos(fx + 2.00 * fy), c3c = std::cos(fx - fy);
  const double s3a = sine(gradient, 2.00 * fx + fy), s3b = sine(gradient, fx + 2.00 * fy),
               s3c = sine(gradient, fx - fy);
  const double c4a = std::cos(2.00 * fx), c4b = std::cos(2.00 * fy), c4c = std::cos(2.00 * (fx + fy));
  const double s4a = sine(star_sines, 2.00 * fx), s4b = sine(star_sines, 2.00 * fy);
  const double s4c = sine(star_sines, 2.00 * (fx + fy));
  const double c5a = std::cos(3.00 * fx + fy), c5b = std::cos(fx + 3.00 * fy), c5c = std::cos(2.00 * fx + 3.00 * fy);
  const double c5d = std::cos(3.00 * fx + 2.00 * fy), c5e = std::cos(-fx + 2.00 * fy), c5f = std::cos(2.00 * fx - fy);
  const double s5a = sine(star_sines, 3.00 * fx + fy), s5b = sine(star_sines, fx + 3.00 * fy);
  const double s5c = sine(star_sines, 2.00 * fx + 3.00 * fy), s5d = sine(star_sines, 3.00 * fx + 2.00 * fy);
  const double s5e = sine(star_sines, -fx + 2.00 * fy), s5f = sine(star_sines, 2.00 * fx - fy);
  const double c6a = std::cos(3.00 * fx), c6b = std::cos(3.00 * fy), c6c = std::cos(3.00 * (fx + fy));
  const double s6a = sine(gradient, 3.00 * fx), s6b = sine(gradient, 3.00 * fy), s6c = sine(gradient, 3.00 * (fx + fy));
  const double c7a = std::cos(2.00 * (2.00 * fx + fy)), c7b = std::cos(2.00 * (fx + 2.00 * fy));
  const double c7c = std::cos(2.00 * (fx - fy));
  const double s7a = sine(gradient, 2.00 * (2.00 * fx + fy)), s7b = sine(gradient, 2.00 * (fx + 2.00 * fy));
  const double s7c = sine(gradient, 2.00 * (fx - fy));

  int k = 0;
  auto put = [&](double f, double dx, double dy) {
    fxy[k] = f;
    dfdx[k] = dx;
    dfdy[k] = dy;
    ++k;
  };
  put(1.00, 0.00, 0.00);
  put(c2a + c2b + c2c, -fac * (s2a + s2c), -fac * (s2b + s2c));
  if (c3v) put(-s2a - s2b + s2c, -fac * (c2a - c2c), -fac * (c2b - c2c));
  put(c3a + c3b + c3c, -fac * (2.00 * s3a + s3b + s3c), -fac * (s3a + 2.00 * s3b - s3c));
  put(c4a + c4b + c4c, -2.00 * fac * (s4a + s4c), -2.00 * fac * (s4b + s4c));
  if (c3v) put(-s4a - s4b + s4c, -2.00 * fac * (c4a - c4c), -2.00 * fac * (c4b - c4c));
  put(c5a + c5b + c5c + c5d + c5e + c5f, -fac * (3.00 * s5a + s5b + 2.00 * s5c + 3.00 * s5d - s5e + 2.00 * s5f),
      -fac * (s5a + 3.00 * s5b + 3.00 * s5c + 2.00 * s5d + 2.00 * s5e - s5f));
  if (c3v) {
    put(-s5a - s5b + s5c + s5d - s5e - s5f, fac * (-3.00 * c5a - c5b + 2.00 * c5c + 3.00 * c5d + c5e - 2.00 * c5f),
        fac * (-c5a - 3.00 * c5b + 3.00 * c5c + 2.00 * c5d - 2.00 * c5e + c5f));
  }
  put(c6a + c6b + c6c, -3.00 * fac * (s6a + s6c), -3.00 * fac * (s6b + s6c));
  put(c7a + c7b + c7c, -2.00 * fac * (2.00 * s7a + s7b + s7c), -2.00 * fac * (s7a + 2.00 * s7b - s7c));

  // Derivatives with respect to Cartesian coordinates
  if (!gradient) return;
  for (int i = 0; i < n_basis_; ++i) dfdy[i] = (-dfdx[i] + 2.00 * dfdy[i]) * usq3;
}

AtomSurfaceValue HexagonalAtomSurface::evaluate(const Vec3& p) const { return compute(p, Want::value_and_gradient); }

double HexagonalAtomSurface::value(const Vec3& p) const { return compute(p, Want::value).v; }

// Fortran: ATII. Want::value: V only (the gradient is meaningless; S3).
AtomSurfaceValue HexagonalAtomSurface::compute(const Vec3& p, Want want) const {
  // Contribution from the 1D repulsive terms, at the position folded into the unit cell
  const auto [x0c, y0c] = fold_hexagonal(p.x, p.y, geometry_.lattice_constant);
  const RepulsiveValue rep =
      prune_cluster_ ? reachable_repulsive_sum(x0c, y0c, p.z, want) : repulsive_sum(x0c, y0c, p.z, want);

  // Interpolation over z for each site (all sites share the z range of the first)
  const std::vector<double>& z = cuts_.front().x();
  if (p.z < z.front() || p.z > z.back()) {
    throw EvaluationError("CRP atom-surface PES: z = " + std::to_string(p.z) + " bohr outside the data range [" +
                          std::to_string(z.front()) + ", " + std::to_string(z.back()) + "] (Fortran ATII)");
  }
  std::array<double, 10> vs{}, dvsdz{};
  for (int i = 0; i < n_basis_; ++i) {
    vs[i] = cuts_[i].value_dpl(p.z);
    if (want == Want::value_and_gradient) dvsdz[i] = cuts_[i].slope_dplp(p.z);
  }

  std::array<double, 10> f{}, dfdx{}, dfdy{};
  fourier_basis(p.x, p.y, f, dfdx, dfdy, want);

  // a(j), b(j): Fourier coefficients of the potential and of its z derivative
  double v = 0.00, dvdx = 0.00, dvdy = 0.00, dvdz = 0.00;
  for (int j = 0; j < n_basis_; ++j) {
    double a = 0.00, b = 0.00;
    for (int ksit = 0; ksit < n_basis_; ++ksit) {
      a = a + table_[ksit + n_basis_ * j] * vs[ksit];
      b = b + table_[ksit + n_basis_ * j] * dvsdz[ksit];
    }
    v = v + a * f[j];
    dvdx = dvdx + a * dfdx[j];
    dvdy = dvdy + a * dfdy[j];
    dvdz = dvdz + b * f[j];
  }
  return {v + rep.v, {dvdx + rep.dvdx, dvdy + rep.dvdy, dvdz + rep.dvdz}};
}

// The distance from the folded unit cell (closed rhombus) to the nearest surface atom that is not in
// the cluster: within it the cluster sum equals the sum over the infinite surface.
double HexagonalAtomSurface::coverage_radius() const {
  const double a = geometry_.lattice_constant;
  const double h = a * sq3 * 0.5;
  const std::array<std::pair<double, double>, 4> cell = hexagonal_cell(a);
  const bool c3v = symmetry_ == Symmetry::c3v;
  double radius = std::numeric_limits<double>::infinity();
  const int layers = c3v ? c3v_cluster_layers : 1;
  for (int k = 0; k < layers; ++k) {
    const double shift = (k % 3) / 3.0;
    for (int i = -20; i <= 20; ++i) {
      for (int j = -20; j <= 20; ++j) {
        bool in_cluster = false;
        if (c3v) {
          in_cluster = i >= -c3v_cluster_n && i <= c3v_cluster_n + 1 && j >= -c3v_cluster_n && j <= c3v_cluster_n + 1;
        } else {
          for (std::size_t c = 0; c < c6v_cluster_i.size(); ++c)
            if (c6v_cluster_i[c] == i && c6v_cluster_j[c] == j) in_cluster = true;
        }
        if (in_cluster) continue;
        const double u = i + shift, v = j + shift;
        radius = std::min(radius, distance_to_polygon((u + 0.5 * v) * a, v * h, cell));
      }
    }
  }
  return radius;
}

std::vector<std::string> HexagonalAtomSurface::warnings() const {
  std::vector<std::string> w;
  const double cutoff = corrugation_.x().back();
  if (cutoff > coverage_radius()) {
    w.push_back("the corrugation function reaches " + std::to_string(cutoff) + " bohr, beyond the " +
                std::to_string(coverage_radius()) +
                " bohr covered by the MURII atom cluster: the repulsive sum is not periodic (V jumps at cell edges)");
  }
  return w;
}

}  // namespace crplib::detail
