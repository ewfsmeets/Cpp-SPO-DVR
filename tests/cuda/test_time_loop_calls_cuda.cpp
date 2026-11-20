// Tests: the backend calls of the time loop with the GPU backend (the twin of test_time_loop_calls.cpp;
// the counts come from the program's time loop, so they must equal the CPU's).
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "counting_backend.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda_test_setup.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;

namespace {

std::string file_bytes(const std::filesystem::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream content;
  content << file.rdbuf();
  return content.str();
}

}  // namespace

TEST(CudaTimeLoopCalls, WithoutPlotKeysAndWithWriteDensities) {
  SPODVR_REQUIRE_GPU();
  for (const bool densities : {false, true}) {
    TemporaryDirectory dir("cuda_time_loop_calls");
    BackendCallCounts counts;
    nlohmann::json json = small_hardwall_input();
    if (densities) json["writeDensities"] = true;
    ASSERT_EQ(run_counted(json, dir.path(), cuda::make_cuda_backend, counts), 0);
    EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps) << densities;
    EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps) << densities;
    EXPECT_EQ(counts.block_densities, small_hardwall_averages_rows) << densities;
    EXPECT_EQ(counts.wave_packet_on_host, 1) << densities;
    EXPECT_EQ(counts.read_wave_packet, 0) << densities;
  }
}

// dtFrames: every Z plane of the scattering grid once per frame (11 frames), as on the CPU.
TEST(CudaTimeLoopCalls, FramesReadEachPlaneOncePerFrame) {
  SPODVR_REQUIRE_GPU();
  TemporaryDirectory dir("cuda_time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["dtFrames"] = 100.0;
  ASSERT_EQ(run_counted(json, dir.path(), cuda::make_cuda_backend, counts), 0);
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps + json["nZ"].get<int>() * 11);
  EXPECT_EQ(counts.block_densities, small_hardwall_averages_rows);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
}

// NormBudget: the block densities once per analysis step and once at tStart, as on the CPU; the averages use those
// densities, so averages.dat is bytewise the GPU run's without the key.
TEST(CudaTimeLoopCalls, NormBudgetReducesOncePerAnalysisStep) {
  SPODVR_REQUIRE_GPU();
  TemporaryDirectory dir("cuda_time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  BackendCallCounts plain_counts;
  ASSERT_EQ(run_counted(json, dir.path() / "plain", cuda::make_cuda_backend, plain_counts), 0);
  json["NormBudget"] = true;
  ASSERT_EQ(run_counted(json, dir.path(), cuda::make_cuda_backend, counts), 0);
  EXPECT_EQ(file_bytes(dir.path() / "averages.dat"), file_bytes(dir.path() / "plain" / "averages.dat"));
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.block_densities, small_hardwall_analysis_steps + 1);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
}
