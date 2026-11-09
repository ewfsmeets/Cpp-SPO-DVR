// Timing of the CRP library (docs/performance.md "CRP library"): the cost of one evaluation of the
// PES as SPO-DVR calls it (QD wrapper, V in eV), on one thread, at points spread like a dynamics
// grid over the cell, Z, r and the orientations.
// Usage: spodvr_crp_timing <data set directory> [repetitions]   (best of the repetitions)
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <numbers>
#include <string>
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
    const crplib::QdWrapper wrapper(data);
    const std::vector<Point> points = dynamics_points(wrapper.parameters().lattice_constant);
    const int repetitions = argc == 3 ? std::stoi(argv[2]) : 3;
    double best = 1e300, checksum = 0.0;
    for (int k = 0; k < repetitions; ++k) {
      const auto start = std::chrono::steady_clock::now();
      double sum = 0.0;
      for (const Point& p : points) sum += wrapper.potential_ev(p.u, p.v, p.z, p.r, p.theta, p.phi);
      const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
      best = std::min(best, elapsed.count());
      checksum = sum;
    }
    std::printf("%s: %zu evaluations, %.3f us per evaluation (best of %d), checksum %.17g\n",
                crplib::symmetry_name(data.symmetry), points.size(), 1e6 * best / static_cast<double>(points.size()),
                repetitions, checksum);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}
