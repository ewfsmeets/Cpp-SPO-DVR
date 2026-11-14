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

#include "app/cli.hpp"
#include "app/run.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

// Runs the program on `json` with output into `dir`; returns the exit code (stderr in err).
int run(const nlohmann::json& json, const fs::path& dir, std::string& err) {
  fs::create_directories(dir);
  const fs::path input = dir / "input.json";
  std::ofstream(input) << json.dump(1);
  CommandLine cl;
  cl.input_file = input;
  cl.output_dir = dir;
  cl.threads = 2;
  std::ostringstream out;
  std::ostringstream errors;
  const int code = run_program(cl, out, errors, false);
  err = errors.str();
  return code;
}

// hardwall_small with a flux analysis that passes every input check: r = 3.05 bohr is a grid
// point below ROptStart, and 0.2..0.6 eV is the wave packet's energy range; 60 steps to t = 120.
nlohmann::json flux_input() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 3.05;
  json["EMinFluxAna"] = 0.2;
  json["EMaxFluxAna"] = 0.6;
  json["DeltaEFluxAna"] = 0.01;
  json["tEnd"] = 120.0;
  json["nStep"] = 60;
  return json;
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

TEST(Cli, Errors) {
  EXPECT_THROW((void)parse_command_line({}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "b.json"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads", "0"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads", "4x"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--bogus"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--schedule"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--schedule", "guided"}), std::invalid_argument);
}

// A new run with the flux analysis runs through and writes its flux files and the flux coefficients
// for a restart; a restart with it runs through too and writes them again (that it continues the
// coefficients bit for bit: test_restart.cpp). An invalid flux input is an input error (exit code 2):
// the input checks come first.
TEST(FluxAnalysisRuns, NewRunAndItsRestartWriteTheFluxFilesAndAnInvalidInputIsRefused) {
  const ScratchDir dir("spodvr_flux_runs");
  std::string err;
  const std::vector<std::string> flux_files = {"flux_reaction.tsv", "flux_density.tsv", "flux_analysis.json",
                                               "flux_coefficients.bin"};
  EXPECT_EQ(run(flux_input(), dir.path / "new", err), 0) << err;
  for (const std::string& name : flux_files) EXPECT_TRUE(fs::exists(dir.path / "new" / name)) << name;

  nlohmann::json restart = flux_input();
  restart["NewCalc"] = false;
  restart["RestartDir"] = (dir.path / "new").string();
  restart["tStart"] = 120.0;
  restart["tEnd"] = 240.0;
  EXPECT_EQ(run(restart, dir.path / "restart", err), 0) << err;
  for (const std::string& name : flux_files) EXPECT_TRUE(fs::exists(dir.path / "restart" / name)) << name;

  nlohmann::json invalid = flux_input();
  invalid["EMaxFluxAna"] = 0.7;
  EXPECT_EQ(run(invalid, dir.path / "invalid", err), 2);
  EXPECT_NE(err.find("EMaxFluxAna"), std::string::npos) << err;
}

// A potential-only run (Potential = true) accepts the same flux keys and finishes normally.
TEST(FluxAnalysisRuns, PotentialOnlyRunAcceptsTheFluxKeys) {
  const ScratchDir dir("spodvr_flux_runs_potential");
  nlohmann::json json = flux_input();
  json["Potential"] = true;
  std::string err;
  EXPECT_EQ(run(json, dir.path, err), 0) << err;
  EXPECT_TRUE(fs::exists(dir.path / "asymptotic_potential.dat"));
}
