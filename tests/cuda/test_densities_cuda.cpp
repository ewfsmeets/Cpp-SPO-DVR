// Tests: densities/ (writeDensities) on the GPU. The records come from the device reduction
// block_densities() and the specular wave (copied exactly), so a run on the GPU writes the CPU's records to
// the CPU-vs-GPU level of the propagation (rule "CPU = GPU"), which grows with the number of steps.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "comparison.hpp"
#include "counting_backend.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda_test_setup.hpp"
#include "io/binary_file.hpp"
#include "propagation/cpu_backend.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

constexpr std::size_t records = 31;            // tStart and every dtAve (100 a.u.) up to tEnd 3000
constexpr double steps_per_record = 50.0;       // dtAve / dtStep

// The tolerance of record k relative to its largest cell: the CPU and GPU propagations drift apart by about
// 1e-13 of the norm per 1000 steps (cuFFT vs oneMKL rounding, docs/cuda.md), and a density doubles the
// relative difference of the amplitude; at least 1e-13 (one analysis step). Measured: 2.2e-13 at 1500 steps.
double record_tolerance(std::size_t k) {
  return std::max(1e-13, 2e-13 * steps_per_record * static_cast<double>(k) / 1000.0);
}

// A coupled case (the hard wall's scattering wave is only rounding noise): cu111_dailight_small on a 4 x 4
// (x, y) grid with JMax = MJMax = 2, until tEnd 3000 a.u., when the packet has reached the surface.
nlohmann::json coupled_input() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["nX"] = 4;
  json["nY"] = 4;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  json["tEnd"] = 3000.0;
  json["nStep"] = 1500;
  json["writeDensities"] = true;
  return json;
}

// Every value of a little-endian binary file of T (a dataset of densities/).
template <typename T>
std::vector<T> read_values(const fs::path& path) {
  std::vector<T> values(fs::file_size(path) / sizeof(T));
  BinaryReader reader(path);
  reader.read_array(std::span(values));
  return values;
}

// The largest |value| of record k (of `size` values).
double record_maximum(const std::vector<double>& values, std::size_t size, std::size_t k) {
  double largest = 0.0;
  for (std::size_t i = k * size; i < (k + 1) * size; ++i)
    largest = spodvr::testing::nan_aware_max(largest, std::abs(values[i]));
  return largest;
}

// The largest |a - b| of each record of `size` values (NaN stays NaN).
std::vector<double> record_differences(const std::vector<double>& a, const std::vector<double>& b, std::size_t size) {
  std::vector<double> worst;
  for (std::size_t start = 0; start + size <= a.size(); start += size) {
    double difference = 0.0;
    for (std::size_t i = start; i < start + size; ++i)
      difference = spodvr::testing::nan_aware_max(difference, std::abs(a[i] - b[i]));
    worst.push_back(difference);
  }
  return worst;
}

// How many records differ by more than their tolerance times their scale (NaN counts).
int failing_records(const std::vector<double>& differences, const std::vector<double>& scale) {
  int failing = 0;
  for (std::size_t k = 0; k < differences.size(); ++k)
    failing += (differences[k] <= record_tolerance(k) * scale[k]) ? 0 : 1;
  return failing;
}

}  // namespace

TEST(CudaDensities, ProgramRunWritesTheCpuRecords) {
  SPODVR_REQUIRE_GPU();
  TemporaryDirectory dir("cuda_densities");
  const nlohmann::json json = coupled_input();
  BackendCallCounts counts;
  ASSERT_EQ(run_counted(json, dir.path() / "cpu", make_cpu_backend, counts), 0);
  ASSERT_EQ(run_counted(json, dir.path() / "gpu", cuda::make_cuda_backend, counts), 0);
  const fs::path cpu = dir.path() / "cpu" / "densities";
  const fs::path gpu = dir.path() / "gpu" / "densities";

  ASSERT_EQ(read_values<double>(cpu / "times.bin").size(), records);
  EXPECT_EQ(read_values<double>(cpu / "times.bin"), read_values<double>(gpu / "times.bin"));
  EXPECT_EQ(read_values<std::int32_t>(cpu / "surface_index.bin"), read_values<std::int32_t>(gpu / "surface_index.bin"));

  const std::vector<double> scat_cpu = read_values<double>(cpu / "prob_scat_zr.bin");
  const std::vector<double> scat_gpu = read_values<double>(gpu / "prob_scat_zr.bin");
  const std::vector<double> spec_cpu = read_values<double>(cpu / "prob_spec_zr.bin");
  const std::vector<double> spec_gpu = read_values<double>(gpu / "prob_spec_zr.bin");
  ASSERT_EQ(scat_cpu.size(), scat_gpu.size());
  ASSERT_EQ(spec_cpu.size(), spec_gpu.size());
  const std::size_t scat_size = scat_cpu.size() / records;
  const std::size_t spec_size = spec_cpu.size() / records;

  // A record's scale: its largest cell of either wave. The scattering wave must carry weight in the last record,
  // else this compares rounding noise.
  std::vector<double> scale;
  for (std::size_t k = 0; k < records; ++k)
    scale.push_back(
        spodvr::testing::nan_aware_max(record_maximum(scat_cpu, scat_size, k), record_maximum(spec_cpu, spec_size, k)));
  ASSERT_GT(record_maximum(scat_cpu, scat_size, records - 1), 1e-8 * scale.back());

  EXPECT_EQ(failing_records(record_differences(scat_cpu, scat_gpu, scat_size), scale), 0);
  EXPECT_EQ(failing_records(record_differences(spec_cpu, spec_gpu, spec_size), scale), 0);

  // Negative: the largest scattering cell of the last record off by 10 times the tolerance, or NaN, fails that
  // record only.
  std::size_t peak = (records - 1) * scat_size;
  for (std::size_t i = peak; i < scat_cpu.size(); ++i)
    if (scat_cpu[i] > scat_cpu[peak]) peak = i;
  for (const double planted : {scat_gpu[peak] + 10 * record_tolerance(records - 1) * scale.back(), std::nan("")}) {
    std::vector<double> wrong = scat_gpu;
    wrong[peak] = planted;
    EXPECT_EQ(failing_records(record_differences(scat_cpu, wrong, scat_size), scale), 1) << planted;
  }
}
