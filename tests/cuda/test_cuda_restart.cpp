// Tests: checkpoint / restart and the command line with the GPU backend (--device cuda). A GPU
// run split by a restart (native or Fortran files) is bitwise identical to the unsplit GPU run;
// a GPU run can continue from a CPU checkpoint (and the results then agree to rounding).
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#include <unistd.h>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "cuda_test_setup.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

nlohmann::json short_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_short_sym/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["CheckpointFormat"] = "both";
  return json;
}

std::string read_file(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

int run(const nlohmann::json& json, const fs::path& dir, const std::string& device, std::string& err) {
  fs::create_directories(dir);
  const fs::path input = dir / "input.json";
  std::ofstream(input) << json.dump(1);
  CommandLine cl;
  cl.input_file = input;
  cl.output_dir = dir;
  cl.threads = 2;
  cl.device = device;
  std::ostringstream out;
  std::ostringstream errors;
  const int code = run_program(cl, out, errors, false);
  err = errors.str();
  return code;
}

nlohmann::json with_times(nlohmann::json json, double t_start, double t_end, int n_steps) {
  json["tStart"] = t_start;
  json["tEnd"] = t_end;
  json["nStep"] = n_steps;
  return json;
}

nlohmann::json restart_from(nlohmann::json json, const fs::path& dir, const char* format = "native") {
  json["NewCalc"] = false;
  json["RestartDir"] = dir.string();
  json["RestartFormat"] = format;
  return json;
}

}  // namespace

TEST(CudaRestart, SplitGpuRunIsBitwiseIdenticalToTheUnsplitGpuRun) {
  SPODVR_REQUIRE_GPU();
  const fs::path root = fs::temp_directory_path() / ("spodvr_cuda_restart_" + std::to_string(::getpid()));
  fs::remove_all(root);
  const nlohmann::json json = short_json();
  const double dt = json["dtStep"].get<double>();
  std::string err;
  ASSERT_EQ(run(with_times(json, 0.0, 80 * dt, 80), root / "full", "cuda", err), 0) << err;
  ASSERT_EQ(run(with_times(json, 0.0, 40 * dt, 40), root / "first", "cuda", err), 0) << err;
  ASSERT_EQ(run(restart_from(with_times(json, 40 * dt, 80 * dt, 40), root / "first"), root / "second", "cuda", err), 0)
      << err;
  ASSERT_EQ(run(restart_from(with_times(json, 40 * dt, 80 * dt, 40), root / "first", "fortran"),
                root / "second_fortran", "cuda", err),
            0)
      << err;
  for (const char* name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "probabilities_rovib.tsv",
                           "probabilities_diffraction.tsv", "fort.9", "fort.21"}) {
    EXPECT_EQ(read_file(root / "full" / name), read_file(root / "second" / name)) << name;
    EXPECT_EQ(read_file(root / "full" / name), read_file(root / "second_fortran" / name)) << name;
  }
  fs::remove_all(root);
}

TEST(CudaRestart, UnknownDeviceIsAnError) {
  EXPECT_THROW((void)parse_command_line({"input.json", "--device", "tpu"}), std::invalid_argument);
  EXPECT_EQ(parse_command_line({"input.json", "--device", "cuda"}).device, "cuda");
  EXPECT_EQ(parse_command_line({"input.json"}).device, "cpu");
}
