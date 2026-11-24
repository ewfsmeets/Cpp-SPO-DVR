// Tests: command-line parsing, and the program's exit codes for inputs it refuses.
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include <nlohmann/json.hpp>

#include "app/backend_factory.hpp"
#include "app/cli.hpp"
#include "app/run.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/fortran_files.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

// Runs the program on `json` with output into `dir` (on `device`); returns the exit code (stderr in err).
int run(const nlohmann::json& json, const fs::path& dir, std::string& err, const std::string& device = "cpu") {
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

// hardwall_small with a flux analysis that passes every input check: r = 3.05 bohr is a grid
// point below ROptStart, and 0.2..0.6 eV is the wave packet's energy range.
nlohmann::json flux_input() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 3.05;
  json["EMinFluxAna"] = 0.2;
  json["EMaxFluxAna"] = 0.6;
  json["DeltaEFluxAna"] = 0.01;
  return json;
}

// The message of the std::invalid_argument that parse_command_line throws for `args` ("" if none).
std::string parse_error(const std::vector<std::string>& args) {
  try {
    (void)parse_command_line(args);
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  ADD_FAILURE() << "no std::invalid_argument for " << ::testing::PrintToString(args);
  return "";
}

// A scratch directory for one test, removed at the end of the test.
struct ScratchDir {
  fs::path path;
  explicit ScratchDir(const std::string& name)
      : path(fs::temp_directory_path() / (name + "_" + std::to_string(::getpid()))) {
    fs::remove_all(path);
  }
  ~ScratchDir() { fs::remove_all(path); }
};

}  // namespace

TEST(Cli, ParsesAllOptions) {
  const auto cl =
      parse_command_line({"in.json", "--threads", "4", "--schedule", "static", "--output-dir", "out", "--no-color"});
  EXPECT_EQ(cl.input_file, "in.json");
  EXPECT_EQ(cl.threads, 4);
  EXPECT_EQ(cl.schedule, "static");
  EXPECT_EQ(cl.output_dir, "out");
  EXPECT_TRUE(cl.no_color);
}

TEST(Cli, HelpAndVersionNeedNoInput) {
  EXPECT_TRUE(parse_command_line({"--help"}).help);
  EXPECT_TRUE(parse_command_line({"--version"}).version);
}

// Each error names its own problem (several share the exception type).
TEST(Cli, Errors) {
  EXPECT_EQ(parse_error({}), "no input file given");
  EXPECT_EQ(parse_error({"a.json", "b.json"}), "more than one input file given");
  EXPECT_EQ(parse_error({"a.json", "--threads"}), "--threads requires a value");
  EXPECT_EQ(parse_error({"a.json", "--threads", "0"}), "--threads requires a positive integer");
  EXPECT_EQ(parse_error({"a.json", "--threads", "4x"}), "--threads requires a positive integer");
  EXPECT_EQ(parse_error({"a.json", "--bogus"}), "unknown option '--bogus'");
  EXPECT_EQ(parse_error({"a.json", "--schedule"}), "--schedule requires a value");
  EXPECT_EQ(parse_error({"a.json", "--schedule", "guided"}), "--schedule must be 'dynamic' or 'static' (got 'guided')");
}

// --device: cpu (the default) or cuda; another value or none is an error. In a build without CUDA a
// propagating run on cuda stops with exit code 3 before its set-up; in a CUDA build the backend exists.
TEST(Cli, DeviceOption) {
  EXPECT_EQ(parse_command_line({"a.json"}).device, "cpu");
  EXPECT_EQ(parse_command_line({"a.json", "--device", "cpu"}).device, "cpu");
  EXPECT_EQ(parse_command_line({"a.json", "--device", "cuda"}).device, "cuda");
  EXPECT_EQ(parse_error({"a.json", "--device", "gpu"}), "--device must be cpu or cuda");
  EXPECT_EQ(parse_error({"a.json", "--device", "CUDA"}), "--device must be cpu or cuda");
  EXPECT_EQ(parse_error({"a.json", "--device"}), "--device requires a value");

  if (cuda_backend_available()) {
    EXPECT_NO_THROW((void)backend_factory(Device::Cuda));
    return;
  }
  const ScratchDir dir("spodvr_cli_device");
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  std::string err;
  EXPECT_EQ(run(nlohmann::json::parse(file), dir.path, err, "cuda"), 3);
  EXPECT_NE(err.find("--device cuda: this program was built without CUDA support"), std::string::npos) << err;
  EXPECT_FALSE(fs::exists(dir.path / "run_info.json"));
}

// The flux analysis is not ported: a propagating run with a valid flux input stops with exit code
// 3 (not supported), names doFluxAnalysis and writes no data file. An invalid flux input is still
// an input error (exit code 2): the input checks come first.
TEST(FluxAnalysisGuard, PropagatingRunIsRefusedWithExitCode3) {
  const ScratchDir dir("spodvr_flux_guard");
  std::string err;
  EXPECT_EQ(run(flux_input(), dir.path / "valid", err), 3);
  EXPECT_NE(err.find("doFluxAnalysis"), std::string::npos) << err;
  EXPECT_FALSE(fs::exists(dir.path / "valid" / "asymptotic_potential.dat"));

  nlohmann::json invalid = flux_input();
  invalid["EMaxFluxAna"] = 0.7;
  EXPECT_EQ(run(invalid, dir.path / "invalid", err), 2);
  EXPECT_NE(err.find("EMaxFluxAna"), std::string::npos) << err;
}

// An analysis-only run (onlyAnalysis) of the Fortran BK file of hardwall_small accepts the same flux
// keys and finishes normally. RestartDir holds that file as fort.21 and a Fortran wave-function file
// fort.9 (zero wave packet, analysis surface not fixed) for the analysis line it reads.
TEST(FluxAnalysisGuard, AnalysisOnlyRunAcceptsTheFluxKeys) {
  const ScratchDir dir("spodvr_flux_guard_analysis");
  const fs::path restart = dir.path / "restart";
  fs::create_directories(restart);
  fs::create_symlink(fs::path(SPODVR_REFERENCE_DIR) / "hardwall_small" / "fortran_fort21.bin", restart / "fort.21");

  nlohmann::json json = flux_input();
  json["onlyAnalysis"] = true;
  json["RestartDir"] = restart.string();
  json["RestartFormat"] = "fortran";
  const Config config = build_config(parse_input(json));
  const InputParameters& p = config.input;
  const SpatialGrids grids = make_spatial_grids(config);
  const RotationalBasis basis(make_angular_grid(config), p.j_in, p.mj_in);
  const WavePacket zero{ScatteringWave(basis.size(), p.nx * p.ny, p.nr * p.nz),
                        SpecularWave(p.nr, grids.z_specular.size())};
  const int before_absorber = grids.z_specular_optical->first - 1;
  write_fortran_wave_function(restart / "fort.9", zero, {before_absorber, before_absorber, false, 0.0});

  std::string err;
  EXPECT_EQ(run(json, dir.path / "out", err), 0) << err;
  EXPECT_EQ(err.find("doFluxAnalysis"), std::string::npos) << err;
  EXPECT_TRUE(fs::exists(dir.path / "out" / "fine_analysis" / "run_info.json"));
  EXPECT_TRUE(fs::exists(dir.path / "out" / "fine_analysis" / "s_matrix.tsv"));  // the fine analysis ran
}

// A potential-only run (Potential = true) accepts the same flux keys and finishes normally.
TEST(FluxAnalysisGuard, PotentialOnlyRunAcceptsTheFluxKeys) {
  const ScratchDir dir("spodvr_flux_guard_potential");
  nlohmann::json json = flux_input();
  json["Potential"] = true;
  std::string err;
  EXPECT_EQ(run(json, dir.path, err), 0) << err;
  EXPECT_TRUE(fs::exists(dir.path / "asymptotic_potential.dat"));
}
