#pragma once
// The runs of the Fortran SCM reference driver in tests/scm/reference/
// (tools/fortran_reference/run_scm_reference.sh): each run's output and the files it read, and
// the evaluation points it was given (points.in of its case).

#include <array>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "scm_driver_output.hpp"

namespace scm_test {

// tests/scm/reference: the driver outputs and the synthetic inputs.
inline const std::filesystem::path reference_directory = SPODVR_SCM_REFERENCE_DIR;
// potentials/H2Cu111SRP_C3v.dir: the example params.dat and examplePosfile.dat.
inline const std::filesystem::path example_directory = SPODVR_SCM_EXAMPLE_DIR;

// One driver run: its output and the files it read.
struct ReferenceRun {
  std::filesystem::path output;      // relative to tests/scm/reference
  std::filesystem::path parameters;  // params.dat as run
  std::filesystem::path slab;        // the position file as run
};

// The four runs: the example slab with rescaling on and off, and the synthetic slabs of versions
// 3 and 2 (a configuration other than the first, layer-2 parameters unlike layer 1's).
inline const std::vector<ReferenceRun>& reference_runs() {
  static const std::vector<ReferenceRun> runs = {
      {"example/rescaling_on.out", example_directory / "params.dat", example_directory / "examplePosfile.dat"},
      {"example/rescaling_off.out", example_directory / "params.dat", example_directory / "examplePosfile.dat"},
      {"synthetic/v3_slab3_rescaling_on.out", reference_directory / "synthetic/params.in",
       reference_directory / "synthetic/slab_v3.in"},
      {"synthetic/v2_slab2_rescaling_off.out", reference_directory / "synthetic/params.in",
       reference_directory / "synthetic/slab_v2.in"},
  };
  return runs;
}

// The point list of a run's case: points.in next to its output.
inline std::filesystem::path points_file(const ReferenceRun& run) {
  return reference_directory / run.output.parent_path() / "points.in";
}

// One evaluation point as the driver read it: u, v (skewed lattice coordinates), Z, r (bohr), θ, φ.
struct EvaluationPoint {
  double u = 0.0, v = 0.0, z = 0.0, r = 0.0, theta = 0.0, phi = 0.0;
};

// Reads a point list: the number of points, then one line of six 16-digit hex words per point.
// Throws std::runtime_error naming the file for a missing file, a short list or a bad word.
inline std::vector<EvaluationPoint> read_points(const std::filesystem::path& file) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("cannot open " + file.string());
  int count = -1;
  in >> count;
  if (count <= 0) throw std::runtime_error(file.string() + ": no point count on line 1");
  std::vector<EvaluationPoint> points;
  for (int k = 0; k < count; ++k) {
    std::array<double, 6> values{};
    for (double& value : values) {
      std::string word;
      if (!(in >> word) || !is_hex_word(word))
        throw std::runtime_error(file.string() + ": point " + std::to_string(k + 1) + " is not six hex words");
      value = from_hex_word(word);
    }
    points.push_back({values[0], values[1], values[2], values[3], values[4], values[5]});
  }
  return points;
}

}  // namespace scm_test
