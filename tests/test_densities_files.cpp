// Tests: densities/ (report/densities_files.*): written with writeDensities, the same bytes with 1 and 8
// threads; the projected plot-data size and the refusal of a run whose plot data exceed the free space.
#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "core/error.hpp"
#include "counting_backend.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "report/densities_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

std::string file_bytes(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream content;
  content << file.rdbuf();
  return content.str();
}

// Runs the program with writeDensities on `threads` OpenMP threads into `dir`.
void run_with_densities(const fs::path& dir, int threads) {
  nlohmann::json json = small_hardwall_input();
  json["writeDensities"] = true;
  fs::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir;
  cl.threads = threads;
  std::ostringstream out;
  std::ostringstream err;
  const int previous = omp_get_max_threads();
  ASSERT_EQ(run_program(cl, out, err, false), 0) << err.str();
  omp_set_num_threads(previous);
  EXPECT_NE(out.str().find("PLOT DATA"), std::string::npos);
}

}  // namespace

// Rule "results do not depend on the thread count": every file of densities/ is bytewise equal.
TEST(DensitiesFiles, OneAndEightThreadsWriteTheSameBytes) {
  TemporaryDirectory dir("densities_threads");
  run_with_densities(dir.path() / "one", 1);
  run_with_densities(dir.path() / "eight", 8);
  for (const char* name :
       {"manifest.json", "times.bin", "prob_scat_zr.bin", "prob_spec_zr.bin", "surface_index.bin", "norm.bin"}) {
    const std::string one = file_bytes(dir.path() / "one" / "densities" / name);
    EXPECT_FALSE(one.empty()) << name;
    EXPECT_EQ(one, file_bytes(dir.path() / "eight" / "densities" / name)) << name;
  }
}

TEST(DensitiesFiles, ProjectedSizeIsRecordsTimesRecordBytes) {
  nlohmann::json json = small_hardwall_input();
  json["writeDensities"] = true;
  const Config config = build_config(parse_input(json));
  const SpatialGrids grids = make_spatial_grids(config);
  const std::vector<PlotDataSize> sizes = plot_data_sizes(config, grids);
  ASSERT_EQ(sizes.size(), 1u);
  const std::uintmax_t record = 8u + 8u * (grids.z_scattering.size() + grids.z_specular.size()) * 56u + 4u + 8u;
  EXPECT_EQ(sizes[0].bytes, small_hardwall_averages_rows * record);
  EXPECT_EQ(sizes[0].key, "writeDensities");
  json.erase("writeDensities");
  EXPECT_TRUE(plot_data_sizes(build_config(parse_input(json)), grids).empty());
}

TEST(DensitiesFiles, RefusesPlotDataLargerThanTheFreeSpace) {
  const std::vector<PlotDataSize> sizes = {{"densities/", "writeDensities", 1000}, {"other/", "otherKey", 10}};
  EXPECT_NO_THROW(check_free_space(sizes, 1010));
  try {
    check_free_space(sizes, 1009);
    FAIL() << "no InputError";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("writeDensities"), std::string::npos) << e.what();
  }
  EXPECT_NO_THROW(check_free_space({}, 0));
}
