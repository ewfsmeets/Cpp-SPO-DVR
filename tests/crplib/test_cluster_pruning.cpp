// Tests: speed-up S2, the atoms of the MURII / MURO cluster beyond the reach of the corrugation
// function left out (crplib::detail::ReachableCluster, reachable_repulsive_sum,
// EvaluationOptions::prune_cluster). Every comparison is bitwise (the sign of zero counts): the
// pruned and the full repulsive sums at dense point sets in and on the fold cell, for cut-offs from
// small to beyond the cluster and every cluster variant; at points constructed to lie at the
// cut-off distance of an atom (to the last ulp); the PES far from the origin; the QD wrapper on a
// dynamics grid; and (SPODVR_LONG_TESTS) ten million points per variant.
// The corrugation function used here does not vanish at its cut-off, so an atom at exactly the
// cut-off distance changes the sum: a wrong boundary decision shows up.
#include <gtest/gtest.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/atom_surface.hpp"
#include "crplib/atom_surface_detail.hpp"
#include "crplib/constants.hpp"
#include "crplib/errors.hpp"
#include "crplib/qd_wrapper.hpp"

using namespace crplib;
using namespace crplib::detail;

namespace {

bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

bool same_bits(const RepulsiveValue& a, const RepulsiveValue& b) {
  return same_bits(a.v, b.v) && same_bits(a.dvdx, b.dvdx) && same_bits(a.dvdy, b.dvdy) && same_bits(a.dvdz, b.dvdz);
}

enum class Variant { c6v, c3v, c4v, c4v_negative_spacing };

std::string name(Variant v) {
  constexpr const char* names[] = {"C6v (1 layer)", "C3v (4 layers)", "C4v (2 layers)", "C4v (relative spacing)"};
  return names[static_cast<int>(v)];
}

DataSet variant_data_set(Variant v) {
  switch (v) {
    case Variant::c6v: return crplib_test::pairwise_data_set(Symmetry::c6v);
    case Variant::c3v: return crplib_test::pairwise_data_set(Symmetry::c3v);
    case Variant::c4v: return crplib_test::pairwise_data_set(Symmetry::c4v);
    case Variant::c4v_negative_spacing: return crplib_test::case_data_set("c4v_negative_spacing");
  }
  return {};
}

// The corrugation function 0.03 exp(-1.2 rho) on [0, cutoff] (25 knots, the last exactly cutoff).
DataSet with_cutoff(DataSet d, double cutoff) {
  const int n = 25;
  d.corrugation.x.assign(n, 0.0);
  d.corrugation.f.assign(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const double rho = i == n - 1 ? cutoff : cutoff * i / (n - 1);
    d.corrugation.x[i] = rho;
    d.corrugation.f[i] = 0.03 * std::exp(-1.2 * rho);
  }
  d.corrugation_slope_first = -0.036;
  d.corrugation_slope_last = -0.036 * std::exp(-1.2 * cutoff);
  return d;
}

// The full and the pruned repulsive sum of one data set, at points of its fold cell.
struct Sums {
  std::function<RepulsiveValue(double, double, double)> full, reachable;
  std::vector<ClusterAtom> cluster;
  int within_reach = 0;
  std::vector<std::pair<double, double>> cell_corners;  // the fold cell (polygon)
  std::vector<std::pair<double, double>> points;        // in and on the fold cell
};

Sums make_sums(const DataSet& d) {
  Sums s;
  const double a = d.lattice_constant;
  if (d.symmetry == Symmetry::c4v) {
    auto surface = std::make_shared<SquareAtomSurface>(d);
    s.full = [surface](double x, double y, double z) { return surface->repulsive_sum(x, y, z); };
    s.reachable = [surface](double x, double y, double z) { return surface->reachable_repulsive_sum(x, y, z); };
    s.cluster = surface->cluster();
    s.within_reach = static_cast<int>(surface->reachable_cluster().atoms().size());
    s.cell_corners = {{0.0, 0.0}, {0.5 * a, 0.0}, {0.5 * a, 0.5 * a}};
    for (int i = 0; i <= 20; ++i)
      for (int j = 0; j <= i; ++j) s.points.emplace_back(0.025 * i * a, 0.025 * j * a);
    // REGION outputs of lattice points, cell edges and diagonals (incl. the negative-integer quirk)
    for (int i = -3; i <= 3; ++i)
      for (double t : {0.0, 0.25, 0.5, 0.75, 1e-13, 1.0 - 1e-13}) {
        for (auto [x, y] : {std::pair{i * a, t * a}, {t * a, i * a}, {i * a + t * a, i * a + t * a}}) {
          auto [x0, y0] = fold_square(x, y, a);
          (void)fold_into_triangle(x0, y0, a);
          s.points.emplace_back(x0, y0);
        }
      }
    return s;
  }
  auto surface = std::make_shared<HexagonalAtomSurface>(d);
  s.full = [surface](double x, double y, double z) { return surface->repulsive_sum(x, y, z); };
  s.reachable = [surface](double x, double y, double z) { return surface->reachable_repulsive_sum(x, y, z); };
  s.cluster = surface->cluster();
  s.within_reach = static_cast<int>(surface->reachable_cluster().atoms().size());
  const double h = a * sq3 * 0.5;
  s.cell_corners = {{0.0, 0.0}, {a, 0.0}, {1.5 * a, h}, {0.5 * a, h}};
  for (int i = 0; i <= 20; ++i)
    for (int j = 0; j <= 20; ++j) {
      const double bigx = 0.05 * i * a, bigy = 0.05 * j * a;
      s.points.emplace_back(bigx + bigy * 0.500, sq3 * bigy * 0.500);  // as ORIGEN2 writes its output
    }
  // ORIGEN2 outputs of lattice points (incl. the negative-integer quirk) and of points next to them
  for (int i = -3; i <= 3; ++i)
    for (int j = -3; j <= 3; ++j)
      for (double t : {0.0, 1e-13, -1e-13, 0.5}) s.points.push_back(fold_hexagonal((i + 0.5 * j + t) * a, j * h, a));
  return s;
}

std::vector<double> heights() {
  std::vector<double> z;
  for (int k = 0; k <= 80; ++k) z.push_back(-12.0 + 0.25 * k);  // below the deepest layer to above the reach
  return z;
}

// The sum, or nullopt where MURII stops (a point exactly on a sub-surface atom).
std::optional<RepulsiveValue> sum_or_stop(const std::function<RepulsiveValue(double, double, double)>& sum, double x,
                                          double y, double z) {
  try {
    return sum(x, y, z);
  } catch (const EvaluationError&) {
    return std::nullopt;
  }
}

// Both sums equal bitwise, or both stop. Returns 1 (a comparison).
int compare_at(const Sums& s, double x, double y, double z, const std::string& where) {
  const std::optional<RepulsiveValue> full = sum_or_stop(s.full, x, y, z);
  const std::optional<RepulsiveValue> reachable = sum_or_stop(s.reachable, x, y, z);
  const std::string at =
      where + " at x = " + std::to_string(x) + ", y = " + std::to_string(y) + ", z = " + std::to_string(z);
  EXPECT_EQ(full.has_value(), reachable.has_value()) << at;
  if (full && reachable) {
    EXPECT_TRUE(same_bits(*full, *reachable)) << at << ": v " << full->v << " vs " << reachable->v;
  }
  return 1;
}

}  // namespace

// The atoms within reach are a subsequence of the cluster, in its order.
TEST(ClusterPruning, ReachableAtomsKeepTheClusterOrder) {
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v, Variant::c4v_negative_spacing}) {
    const DataSet d = variant_data_set(v);
    const Sums s = make_sums(d);
    const ReachableCluster reachable(s.cluster, s.cell_corners, d.corrugation.x.back(), d.lattice_constant);
    std::size_t next = 0;
    for (const ClusterAtom& atom : reachable.atoms()) {
      while (next < s.cluster.size() && !(s.cluster[next].x == atom.x && s.cluster[next].y == atom.y &&
                                          s.cluster[next].z == atom.z && s.cluster[next].layer == atom.layer))
        ++next;
      ASSERT_LT(next, s.cluster.size()) << name(v) << ": not a subsequence in cluster order";
      ++next;
    }
    EXPECT_EQ(static_cast<int>(reachable.atoms().size()), s.within_reach) << name(v);
  }
  // The synthetic data sets (cut-off 6 bohr): MURII 30 / 256 atoms, MURO 32.
  EXPECT_EQ(make_atom_surface(variant_data_set(Variant::c6v))->cluster_size().atoms, 30);
  EXPECT_EQ(make_atom_surface(variant_data_set(Variant::c3v))->cluster_size().atoms, 256);
  EXPECT_EQ(make_atom_surface(variant_data_set(Variant::c3v))->cluster_size().within_reach, 54);
  EXPECT_EQ(make_atom_surface(variant_data_set(Variant::c4v))->cluster_size().atoms, 32);
}

// A cut-off of zero or below disables the pruning (no reach to compute).
TEST(ClusterPruning, NonPositiveCutoffKeepsEverything) {
  const std::vector<ClusterAtom> cluster = {{0.0, 0.0, 0.0, 0}, {50.0, 0.0, 0.0, 0}, {0.0, 0.0, -90.0, 1}};
  const std::vector<std::pair<double, double>> cell = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}};
  EXPECT_EQ(ReachableCluster(cluster, cell, 0.0, 1.0).atoms().size(), 3U);
  EXPECT_EQ(ReachableCluster(cluster, cell, -2.0, 1.0).atoms().size(), 3U);
  EXPECT_EQ(ReachableCluster(cluster, cell, 2.0, 1.0).atoms().size(), 2U);  // in the plane only; z per call
}

// Dense points in and on the fold cell, z from below the deepest layer to above the reach, for
// cut-offs from a third of a lattice constant to beyond the cluster (where nothing is pruned).
TEST(ClusterPruning, BitwiseOnDenseCellPoints) {
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v, Variant::c4v_negative_spacing}) {
    const DataSet base = variant_data_set(v);
    const double a = base.lattice_constant;
    const double coverage = base.symmetry == Symmetry::c4v ? SquareAtomSurface(base).coverage_radius()
                                                           : HexagonalAtomSurface(base).coverage_radius();
    for (double cutoff : {0.35 * a, 1.0 * a, 6.0, coverage, coverage + 0.5 * a, 10.0 * a}) {
      const Sums s = make_sums(with_cutoff(base, cutoff));
      if (cutoff < coverage) {
        EXPECT_LT(s.within_reach, static_cast<int>(s.cluster.size())) << name(v);
      }
      if (cutoff == 10.0 * a) {
        EXPECT_EQ(s.within_reach, static_cast<int>(s.cluster.size())) << name(v);
      }
      int n = 0;
      for (const auto& [x, y] : s.points)
        for (double z : heights()) n += compare_at(s, x, y, z, name(v) + ", cut-off " + std::to_string(cutoff));
      EXPECT_GT(n, 10000);
    }
  }
}

// Points at the cut-off distance of an atom: above and below each atom at a cell corner or edge
// point, at the z where the computed distance is the cut-off, and 1 to 6 ulps either side (the
// per-call z test and the in-plane selection both decide there).
TEST(ClusterPruning, BitwiseAtTheCutoffDistance) {
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v, Variant::c4v_negative_spacing}) {
    const DataSet base = variant_data_set(v);
    for (double cutoff : {1.0 * base.lattice_constant, 1.7 * base.lattice_constant}) {
      const Sums s = make_sums(with_cutoff(base, cutoff));
      std::vector<std::pair<double, double>> probes = s.cell_corners;
      for (std::size_t c = 0; c < s.cell_corners.size(); ++c) {  // edge midpoints
        const auto [ax, ay] = s.cell_corners[c];
        const auto [bx, by] = s.cell_corners[(c + 1) % s.cell_corners.size()];
        probes.emplace_back(0.5 * (ax + bx), 0.5 * (ay + by));
      }
      int n = 0;
      for (const ClusterAtom& atom : s.cluster)
        for (const auto& [x, y] : probes) {
          const double rho2 = (x - atom.x) * (x - atom.x) + (y - atom.y) * (y - atom.y);
          if (!(rho2 <= cutoff * cutoff)) continue;
          const double dz = std::sqrt(cutoff * cutoff - rho2);
          for (double z0 : {atom.z + dz, atom.z - dz}) {
            double z = z0;
            for (int k = 0; k < 6; ++k) z = std::nextafter(z, -1e300);
            for (int k = -6; k <= 6; ++k, z = std::nextafter(z, 1e300)) n += compare_at(s, x, y, z, name(v));
          }
        }
      EXPECT_GT(n, 100) << name(v);
    }
  }
}

// The cut-off equal to the computed in-plane distance from a cell corner to a top-layer atom: at the
// corner (z = 0) that atom is exactly at the cut-off (it counts in the full sum); then 1 to 3 ulps
// away from the corner in x and y. Tests the in-plane selection with its margins.
TEST(ClusterPruning, BitwiseWithTheCutoffAtAnAtomDistance) {
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v}) {
    const DataSet base = variant_data_set(v);
    const double a = base.lattice_constant;
    const Sums plain = make_sums(base);
    int cases = 0;
    for (const auto& [px, py] : plain.cell_corners)
      for (const ClusterAtom& atom : plain.cluster) {
        if (atom.layer != 0) continue;
        const double dx = px - atom.x, dy = py - atom.y;
        const double rho = std::sqrt(dx * dx + dy * dy + 0.0 * 0.0);  // as the sums compute r at z = 0
        if (rho < 0.4 * a || rho > 2.5 * a || cases >= 24) continue;
        ++cases;
        const Sums s = make_sums(with_cutoff(base, rho));
        double x = px, y = py;
        for (int k = 0; k < 3; ++k) x = std::nextafter(x, -1e300);
        for (int i = -3; i <= 3; ++i, x = std::nextafter(x, 1e300)) {
          compare_at(s, x, py, 0.0, name(v) + ", cut-off at an atom distance");
          y = py;
          for (int k = 0; k < 3; ++k) y = std::nextafter(y, -1e300);
          for (int j = -3; j <= 3; ++j, y = std::nextafter(y, 1e300))
            compare_at(s, x, y, 0.0, name(v) + ", cut-off at an atom distance");
        }
      }
    EXPECT_GE(cases, 10) << name(v);
  }
}

// The PES with and without pruning, also far from the origin, where the fold's rounding can push a
// point out of the cell (then the full cluster is used).
TEST(ClusterPruning, AtomSurfaceBitwiseNearAndFarFromTheOrigin) {
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v, Variant::c4v_negative_spacing}) {
    const DataSet d = variant_data_set(v);
    const auto pruned = make_atom_surface(d);
    const auto full = make_atom_surface(d, EvaluationOptions{.prune_cluster = false});
    const double a = d.lattice_constant;
    for (double scale : {0.0, 1.0, 1e3, 1e6, 1e8})
      for (int i = 0; i < 40; ++i) {
        const double x = scale * a * (0.37 + 0.11 * i) * (i % 2 == 0 ? 1.0 : -1.0) + 0.013 * i;
        const double y = scale * a * (0.71 - 0.07 * i) + 0.29 * i;
        for (double z : {-1.0, 0.3, 1.7, 4.0, 8.5}) {
          const AtomSurfaceValue p = pruned->evaluate({x, y, z}), f = full->evaluate({x, y, z});
          EXPECT_TRUE(same_bits(p.v, f.v) && same_bits(p.gradient.x, f.gradient.x) &&
                      same_bits(p.gradient.y, f.gradient.y) && same_bits(p.gradient.z, f.gradient.z))
              << name(v) << " at " << x << ", " << y << ", " << z;
        }
      }
  }
}

// The PES as SPO-DVR evaluates it, on a dynamics grid (6 x 6 cell points, Z = -2..8, r = 0.8..3.3,
// 18 orientations), with and without pruning.
TEST(ClusterPruning, QdWrapperBitwiseOnADynamicsGrid) {
  for (Symmetry symmetry : {Symmetry::c6v, Symmetry::c3v, Symmetry::c4v}) {
    const DataSet d = crplib_test::pairwise_data_set(symmetry);
    const QdWrapper pruned(d), full(d, EvaluationOptions{.prune_cluster = false});
    const double a = d.lattice_constant;
    int differences = 0, n = 0;
    for (int iu = 0; iu < 6; ++iu)
      for (int iv = 0; iv < 6; ++iv)
        for (int iz = 0; iz <= 40; ++iz)
          for (int ir = 0; ir <= 10; ++ir)
            for (double theta : {0.4, 1.0, std::numbers::pi / 2})
              for (int ip = 0; ip < 6; ++ip) {
                const double u = iu * a / 6, v = iv * a / 6, z = -2.0 + 0.25 * iz, r = 0.8 + 0.25 * ir;
                const double phi = ip * std::numbers::pi / 3;
                if (!same_bits(pruned.potential_ev(u, v, z, r, theta, phi), full.potential_ev(u, v, z, r, theta, phi)))
                  ++differences;
                ++n;
              }
    EXPECT_EQ(differences, 0) << symmetry_name(symmetry) << ": " << differences << " of " << n << " points differ";
  }
}

// SPODVR_LONG_TESTS: ten million points per cluster variant (a fixed pseudo-random sequence over a
// region of 40 x 40 lattice constants and z = -15..10), folded as the PES folds them.
TEST(ClusterPruningLong, TenMillionPoints) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  for (Variant v : {Variant::c6v, Variant::c3v, Variant::c4v, Variant::c4v_negative_spacing}) {
    const DataSet d = variant_data_set(v);
    const Sums s = make_sums(d);
    const double a = d.lattice_constant;
    std::uint64_t state = 0x9E3779B97F4A7C15ULL;
    auto next = [&state] {  // 53 random bits in [0, 1)
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      return static_cast<double>(state >> 11) * 0x1.0p-53;
    };
    int differences = 0;
    for (int i = 0; i < 10'000'000; ++i) {
      const double x = (next() - 0.5) * 40.0 * a, y = (next() - 0.5) * 40.0 * a, z = -15.0 + 25.0 * next();
      double x0 = 0.0, y0 = 0.0;
      if (d.symmetry == Symmetry::c4v) {
        std::tie(x0, y0) = fold_square(x, y, a);
        (void)fold_into_triangle(x0, y0, a);
      } else {
        std::tie(x0, y0) = fold_hexagonal(x, y, a);
      }
      if (!same_bits(s.full(x0, y0, z), s.reachable(x0, y0, z))) ++differences;
    }
    EXPECT_EQ(differences, 0) << name(v);
  }
}
