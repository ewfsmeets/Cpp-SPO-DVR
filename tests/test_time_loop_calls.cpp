// Tests: the backend calls of the time loop, counted through run_program (tests/counting_backend.hpp):
// without any plot key exactly the calls of the program before the plot data existed, and each plot key
// adds only the calls its design lists (docs/CONTRIBUTING.md, "Additions to the time loop").
#include <gtest/gtest.h>

#include "counting_backend.hpp"
#include "propagation/cpu_backend.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;

namespace {

// The calls without plot keys: one propagation and one BK plane per analysis step, the block densities
// once per averages row, and the wave packet on the host once after the loop (the native checkpoint).
void expect_base_counts(const BackendCallCounts& c) {
  EXPECT_EQ(c.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(c.z_plane, small_hardwall_analysis_steps);
  EXPECT_EQ(c.block_densities, small_hardwall_averages_rows);
  EXPECT_EQ(c.wave_packet_on_host, 1);
  EXPECT_EQ(c.read_wave_packet, 0);
}

}  // namespace

TEST(TimeLoopCalls, WithoutPlotKeys) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  ASSERT_EQ(run_counted(small_hardwall_input(), dir.path(), make_cpu_backend, counts), 0);
  expect_base_counts(counts);
}

// writeDensities writes the densities the averages already computed: no extra call.
TEST(TimeLoopCalls, WriteDensitiesAddsNoCall) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["writeDensities"] = true;
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  expect_base_counts(counts);
  EXPECT_TRUE(std::filesystem::exists(dir.path() / "densities" / "times.bin"));
}

// dtFrames: each frame (tStart and every 100 a.u.: 11) reads every Z plane of the scattering grid once, and
// nothing else; the wave packet is never copied to the host in the loop.
TEST(TimeLoopCalls, FramesReadEachPlaneOncePerFrame) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["dtFrames"] = 100.0;
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  const int frames = 11;
  const int nz = json["nZ"].get<int>();
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps + nz * frames);
  EXPECT_EQ(counts.block_densities, small_hardwall_averages_rows);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
  EXPECT_EQ(counts.read_wave_packet, 0);
}

// NormBudget: the block densities once per analysis step and once at tStart, whatever dtAve is (the averages rows
// use the same densities: one reduction, not two); nothing else.
TEST(TimeLoopCalls, NormBudgetReducesOncePerAnalysisStep) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["NormBudget"] = true;
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.block_densities, small_hardwall_analysis_steps + 1);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
  EXPECT_EQ(counts.read_wave_packet, 0);
  EXPECT_TRUE(std::filesystem::exists(dir.path() / "norm_budget" / "times.bin"));
}

// FrameExtras work inside the frames' pass over the planes: no backend call beyond the frames'.
TEST(TimeLoopCalls, FrameExtrasAddNoCall) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["dtFrames"] = 100.0;
  json["FrameExtras"] = {"angular", "diffraction", "vibration", "rotation_moments", "momentum"};
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
  EXPECT_EQ(counts.z_plane, small_hardwall_analysis_steps + json["nZ"].get<int>() * 11);
  EXPECT_EQ(counts.block_densities, small_hardwall_averages_rows);
  EXPECT_EQ(counts.wave_packet_on_host, 1);
}

// The Fortran checkpoint files: the wave packet on the host twice after the loop, never in it.
TEST(TimeLoopCalls, BothCheckpointFormatsReadTheWavePacketTwiceAfterTheLoop) {
  TemporaryDirectory dir("time_loop_calls");
  BackendCallCounts counts;
  nlohmann::json json = small_hardwall_input();
  json["CheckpointFormat"] = "both";
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  EXPECT_EQ(counts.wave_packet_on_host, 2);
  EXPECT_EQ(counts.propagate_analysis_step, small_hardwall_analysis_steps);
}
