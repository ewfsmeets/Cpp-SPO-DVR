// Tests: the backend calls of the time loop with the GPU backend (the twin of test_time_loop_calls.cpp;
// the counts come from the program's time loop, so they must equal the CPU's).
#include <gtest/gtest.h>

#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "comparison.hpp"
#include "counting_backend.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda_test_setup.hpp"
#include "propagation/cpu_backend.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using spodvr::testing::nan_aware_max;
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

// EnergyResolvedSites/Energies: one site column per site and analysis step, as on the CPU; the whole run's sums (of
// the scattering and the specular wave) agree with the CPU's to 1e-12 of their largest magnitude (CPU vs GPU drift
// about 1e-13 per 1000 steps, docs/cuda.md).
TEST(CudaTimeLoopCalls, EnergyResolvedSumsReadOneColumnPerSiteAndStepAndMatchTheCpu) {
  SPODVR_REQUIRE_GPU();
  nlohmann::json json = small_hardwall_input();
  json["PlotSites"] = {{{"name", "a"}, {"x", 0.0}, {"y", 0.0}}, {{"name", "b"}, {"x", 0.5}, {"y", 0.0}}};
  json["EnergyResolvedSites"] = {"a", "b"};
  json["EnergyResolvedEnergies"] = {0.3, 0.4, 0.5};
  TemporaryDirectory gpu_dir("cuda_time_loop_calls");
  BackendCallCounts counts;
  ASSERT_EQ(run_counted(json, gpu_dir.path(), cuda::make_cuda_backend, counts), 0);
  EXPECT_EQ(counts.site_column, 2 * small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
  TemporaryDirectory cpu_dir("cuda_time_loop_calls");
  BackendCallCounts cpu_counts;
  ASSERT_EQ(run_counted(json, cpu_dir.path(), make_cpu_backend, cpu_counts), 0);
  // Both files are parts of one wave packet: their differences are compared with the largest magnitude of either
  // (on the hard wall nothing couples into the scattering wave, whose sums are rounding noise, about 1e-27).
  std::vector<std::vector<std::complex<double>>> cpu_values, gpu_values;
  for (const std::string name : {"sums.bin", "specular_sums.bin"}) {
    const std::string a = file_bytes(cpu_dir.path() / "energy_resolved" / name);
    const std::string b = file_bytes(gpu_dir.path() / "energy_resolved" / name);
    ASSERT_EQ(a.size(), b.size()) << name;
    cpu_values.emplace_back(a.size() / 16);
    gpu_values.emplace_back(b.size() / 16);
    std::memcpy(cpu_values.back().data(), a.data(), a.size());
    std::memcpy(gpu_values.back().data(), b.data(), b.size());
  }
  double largest = 0.0, worst = 0.0;
  for (std::size_t f = 0; f < cpu_values.size(); ++f)
    for (std::size_t k = 0; k < cpu_values[f].size(); ++k) {
      largest = nan_aware_max(largest, std::abs(cpu_values[f][k]));
      worst = nan_aware_max(worst, std::abs(cpu_values[f][k] - gpu_values[f][k]));
    }
  EXPECT_GT(largest, 1.0);
  EXPECT_LE(worst, 1e-12 * largest);
}
