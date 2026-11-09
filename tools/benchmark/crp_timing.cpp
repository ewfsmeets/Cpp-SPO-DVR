// Timing of the CRP library (docs/performance.md "CRP library"): the cost of one evaluation of the
// PES as SPO-DVR calls it (QD wrapper, V in eV), on one thread, at points spread like a dynamics
// grid over the cell, Z, r and the orientations; with the speed-ups (crplib::EvaluationOptions)
// off and on. Every V must be bitwise the same.
// Usage: spodvr_crp_timing <data set directory> [repetitions]   (best of the repetitions)
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "crplib/dataset.hpp"
#include "crplib/qd_wrapper.hpp"

namespace {

struct Point {
  double u, v, z, r, theta, phi;
};

// 36 centres in the cell, Z from -1 to 8 bohr, r from 0.6 to 4.1 bohr, 3 x 6 orientations.
std::vector<Point> dynamics_points(double lattice_constant) {
  std::vector<Point> points;
  for (int iu = 0; iu < 6; ++iu)
    for (int iv = 0; iv < 6; ++iv)
      for (int iz = 0; iz < 19; ++iz)
        for (int ir = 0; ir < 8; ++ir)
          for (int it = 0; it < 3; ++it)
            for (int ip = 0; ip < 6; ++ip)
              points.push_back({iu * lattice_constant / 6.0, iv * lattice_constant / 6.0, -1.0 + 0.5 * iz,
                                0.6 + 0.5 * ir, (it + 0.5) * std::numbers::pi / 6.0, ip * std::numbers::pi / 3.0});
  return points;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2 || argc > 3) {
    std::fprintf(stderr, "usage: spodvr_crp_timing <data set directory> [repetitions]\n");
    return 1;
  }
  try {
    const crplib::DataSet data = crplib::read_data_set(argv[1]);
    const int repetitions = argc == 3 ? std::stoi(argv[2]) : 3;
    std::vector<double> reference;
    using Options = crplib::EvaluationOptions;
    for (const auto& [name, options] : {std::pair{"plain", Options{.prune_cluster = false, .share_grid_lookup = false}},
                                        std::pair{"S2", Options{.prune_cluster = true, .share_grid_lookup = false}},
                                        std::pair{"S4", Options{.prune_cluster = false, .share_grid_lookup = true}},
                                        std::pair{"default", Options{}}}) {
      const crplib::QdWrapper wrapper(data, options);
      const std::vector<Point> points = dynamics_points(wrapper.parameters().lattice_constant);
      std::vector<double> v(points.size());
      double best = 1e300;
      for (int k = 0; k < repetitions; ++k) {
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < points.size(); ++i) {
          const Point& p = points[i];
          v[i] = wrapper.potential_ev(p.u, p.v, p.z, p.r, p.theta, p.phi);
        }
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
        best = std::min(best, elapsed.count());
      }
      if (reference.empty()) reference = v;
      const bool same = std::equal(v.begin(), v.end(), reference.begin(), [](double a, double b) {
        return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
      });
      std::printf("%s %-8s %zu evaluations, %.3f us per evaluation (best of %d)%s\n",
                  crplib::symmetry_name(data.symmetry), name, points.size(),
                  1e6 * best / static_cast<double>(points.size()), repetitions, same ? "" : "  V DIFFERS FROM plain");
      if (!same) return 2;
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
