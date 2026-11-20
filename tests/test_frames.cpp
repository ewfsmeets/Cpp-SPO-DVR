// Tests: frames/ (dtFrames; frames/frame_reducers.*, frames/frame_source.*, report/frame_files.*): the reductions
// of one plane against a hand computation, with a negative test; the same bytes with 1 and 8 threads; the
// projected size; the memory-estimate item. The physics checks (sums, orthogonality, numpy against the
// checkpoints, restart joins) run on the plot fixtures (tools/plots/tests/test_runs_frames.py).
#include <gtest/gtest.h>

#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "counting_backend.hpp"
#include "frames/frame_source.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/memory_estimate.hpp"
#include "report/frame_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;
using Complex = std::complex<double>;

namespace {

std::string file_bytes(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream content;
  content << file.rdbuf();
  return content.str();
}

// The small hard-wall input of the call tests with frames every dtAve (100 a.u.).
nlohmann::json frames_input() {
  nlohmann::json json = small_hardwall_input();
  json["dtFrames"] = 100.0;
  return json;
}

// Runs the program on `json` with `threads` OpenMP threads into `dir`.
std::string run_into(const nlohmann::json& json, const fs::path& dir, int threads) {
  fs::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir;
  cl.threads = threads;
  std::ostringstream out;
  std::ostringstream err;
  const int previous = omp_get_max_threads();
  const int status = run_program(cl, out, err, false);
  omp_set_num_threads(previous);
  EXPECT_EQ(status, 0) << err.str();
  return out.str();
}

// A small layout: 3 rotational states (the initial one is 1), 2 x 2 (x, y) points, 2 r points, 3 + 2 Z points.
struct SyntheticFrame {
  std::vector<Complex> s2s = {{0.5, 0.1}, {-0.2, 0.4}, {0.3, -0.3}, {0.1, 0.2}};
  FrameLayout layout{3, 5, 2, 4, 3, 1, 0.25, 0.5, s2s};
  SpecularWave specular{2, 5};
  std::vector<Complex> plane = std::vector<Complex>(layout.plane().size());

  SyntheticFrame() {
    for (std::size_t k = 0; k < plane.size(); ++k) plane[k] = {0.1 * double(k % 7) - 0.3, 0.05 * double(k % 5)};
    for (int iz = 0; iz < 5; ++iz)
      for (int ir = 0; ir < 2; ++ir) specular(ir, iz) = {0.2 + 0.1 * iz, -0.1 * ir + 0.05 * iz};
    specular(1, 2) = {0.0, 0.0};  // one point without chi
  }
};

}  // namespace

// One plane by hand: prob_zrxy = |psi_scat + s2s chi (initial state only)|^2 dZ dr dx dy sin, pop_jmj_z its sum
// over r and (x, y) per state, the tail |chi|^2 dZ dr beyond the scattering grid (|chi|^2 as abs(chi)^2, as the
// densities' specular cells, bitwise). Negative test: chi left out of
// the initial state (a reducer that forgets psi_full) changes prob_zrxy at exactly the (r, x, y) points of the
// plane where chi is not zero.
TEST(FrameReducers, OnePlaneByHand) {
  SyntheticFrame f;
  const int iz = 2;
  FrameRecord record = make_frame_record(f.layout);
  std::vector<double> scratch(std::size_t(f.layout.nr) * f.layout.n_rot);
  reduce_scattering_plane(f.layout, iz, f.plane, f.specular, scratch, record);
  reduce_specular_tail(f.layout, f.specular, record);
  std::vector<double> pop(3, 0.0);
  for (int ir = 0; ir < 2; ++ir)
    for (int ixy = 0; ixy < 4; ++ixy) {
      double cell = 0.0;
      for (int irot = 0; irot < 3; ++irot) {
        Complex psi = f.plane[f.layout.plane()(irot, ixy, ir)];
        if (irot == 1) psi += f.s2s[ixy] * f.specular(ir, iz);
        cell += std::norm(psi);
        pop[irot] += std::norm(psi) * 0.25;
      }
      EXPECT_NEAR(record.prob_zrxy[f.layout.zrxy()(ixy, ir, iz)], cell * 0.25, 1e-15) << ir << " " << ixy;
    }
  for (int irot = 0; irot < 3; ++irot) EXPECT_NEAR(record.pop_jmj_z[f.layout.jmj_z()(irot, iz)], pop[irot], 1e-15);
  for (int izt = 3; izt < 5; ++izt)
    for (int ir = 0; ir < 2; ++ir)
      EXPECT_EQ(record.prob_tail_zr[f.layout.tail_zr()(ir, izt - 3)],
                std::abs(f.specular(ir, izt)) * std::abs(f.specular(ir, izt)) * 0.5);  // as densities/ squares

  FrameLayout slip = f.layout;
  slip.initial_channel = -1;  // chi never added
  FrameRecord without_chi = make_frame_record(slip);
  reduce_scattering_plane(slip, iz, f.plane, f.specular, scratch, without_chi);
  int changed = 0;
  for (int ir = 0; ir < 2; ++ir)
    for (int ixy = 0; ixy < 4; ++ixy) {
      const std::size_t k = f.layout.zrxy()(ixy, ir, iz);
      changed += (std::abs(without_chi.prob_zrxy[k] - record.prob_zrxy[k]) <= 1e-15) ? 0 : 1;
    }
  EXPECT_EQ(changed, 4);  // r index 0 only: chi(1, 2) = 0
}

// The reductions do not depend on the thread count: a run with 1 and with 8 threads writes the same bytes.
TEST(FrameFiles, OneAndEightThreadsWriteTheSameBytes) {
  TemporaryDirectory dir("frames_threads");
  run_into(frames_input(), dir.path() / "one", 1);
  const std::string out = run_into(frames_input(), dir.path() / "eight", 8);
  EXPECT_NE(out.find("Frames (dtFrames): 11 written"), std::string::npos) << out;
  for (const std::string name :
       {"manifest.json", "times.bin", "prob_zrxy.bin", "prob_tail_zr.bin", "pop_jmj_z.bin", "surface_index.bin"}) {
    const std::string one = file_bytes(dir.path() / "one" / "frames" / name);
    EXPECT_FALSE(one.empty()) << name;
    EXPECT_EQ(one, file_bytes(dir.path() / "eight" / "frames" / name)) << name;
  }
}

// The projected size is the written size: 11 records (tStart and every 100 a.u. to 1000) of the datasets and the
// time; the PLOT DATA line names the key.
TEST(FrameFiles, ProjectedSizeIsTheWrittenSize) {
  TemporaryDirectory dir("frames_size");
  const std::string out = run_into(frames_input(), dir.path(), 8);
  EXPECT_NE(out.find("dtFrames"), std::string::npos);
  const Config config = build_config(parse_input(frames_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const auto size = frames_data_size(config, grids, basis);
  ASSERT_TRUE(size.has_value());
  std::uintmax_t written = 0;
  for (const auto& entry : fs::directory_iterator(dir.path() / "frames"))
    if (entry.path().extension() == ".bin") written += fs::file_size(entry.path());
  EXPECT_EQ(size->bytes, written);
}

// With dtFrames the memory estimate lists the frame record and its per-r sums: the sizes FrameSource allocates.
TEST(FrameFiles, MemoryEstimateItem) {
  const Config config = build_config(parse_input(frames_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const MemoryEstimate e = estimate_memory(config, grids, angular, make_grid_symmetry(config, angular), 1);
  ASSERT_EQ(e.items.back().name, "Frame record (dtFrames)");
  const FrameLayout layout = frame_layout(config, grids, basis, {});
  const FrameRecord record = make_frame_record(layout);
  EXPECT_EQ(e.items.back().values, record.prob_zrxy.size() + record.prob_tail_zr.size() + record.pop_jmj_z.size() +
                                       std::size_t(layout.nr) * layout.n_rot);
  const Config without = build_config(parse_input(small_hardwall_input()));
  const MemoryEstimate base = estimate_memory(without, grids, angular, make_grid_symmetry(without, angular), 1);
  EXPECT_EQ(base.items.size() + 1, e.items.size());
}
