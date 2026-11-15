// Tests: checkpoint / restart. A run split by a restart must give bitwise the same results as
// the unsplit run; damaged or mismatching checkpoints are clear errors.
#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <unistd.h>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "io/fortran_unformatted.hpp"
#include "scm_example_input.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

nlohmann::json short_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_short/input.json");
  return nlohmann::json::parse(file);
}

std::string read_file(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

// Runs the program on `json` with output into `dir`; returns the exit code (stderr in err,
// stdout in out).
int run(const nlohmann::json& json, const fs::path& dir, std::string& err, std::string& out) {
  fs::create_directories(dir);
  const fs::path input = dir / "input.json";
  std::ofstream(input) << json.dump(1);
  CommandLine cl;
  cl.input_file = input;
  cl.output_dir = dir;
  cl.threads = 2;
  std::ostringstream output;
  std::ostringstream errors;
  const int code = run_program(cl, output, errors, false);
  err = errors.str();
  out = output.str();
  return code;
}

int run(const nlohmann::json& json, const fs::path& dir, std::string& err) {
  std::string out;
  return run(json, dir, err, out);
}

nlohmann::json with_times(nlohmann::json json, double t_start, double t_end, int n_steps) {
  json["tStart"] = t_start;
  json["tEnd"] = t_end;
  json["nStep"] = n_steps;
  return json;
}

nlohmann::json restart_from(nlohmann::json json, const fs::path& dir) {
  json["NewCalc"] = false;
  json["RestartDir"] = dir.string();
  return json;
}

// Scratch directory with the unsplit run (0 .. 120) and the split run (0 .. 60, 60 .. 120),
// the latter restarted once from the native and once from the Fortran files.
struct RestartRuns {
  fs::path root = fs::temp_directory_path() / ("spodvr_restart_test_" + std::to_string(::getpid()));
  int code_full;
  int code_first;
  int code_second;
  int code_second_fortran;
  std::string err;

  RestartRuns() {
    fs::remove_all(root);
    nlohmann::json json = short_json();
    json["CheckpointFormat"] = "both";
    code_full = run(with_times(json, 0.0, 120.0, 60), root / "full", err);
    code_first = run(with_times(json, 0.0, 60.0, 30), root / "first", err);
    code_second = run(restart_from(with_times(json, 60.0, 120.0, 30), root / "first"), root / "second", err);
    nlohmann::json fortran = restart_from(with_times(json, 60.0, 120.0, 30), root / "first");
    fortran["RestartFormat"] = "fortran";
    code_second_fortran = run(fortran, root / "second_fortran", err);
  }
  ~RestartRuns() { fs::remove_all(root); }
};

const RestartRuns& runs() {
  static const RestartRuns r;
  return r;
}

// Only the first part of the split run (for the error cases), in its own scratch directory.
struct FirstRun {
  fs::path root = fs::temp_directory_path() / ("spodvr_restart_errors_" + std::to_string(::getpid()));
  int code;
  std::string err;

  FirstRun() {
    fs::remove_all(root);
    code = run(with_times(short_json(), 0.0, 60.0, 30), root / "first", err);
  }
  ~FirstRun() { fs::remove_all(root); }
};

const FirstRun& first_run() {
  static const FirstRun r;
  return r;
}

// Fine analysis (alsoAnalysis / onlyAnalysis) on a short run with a large time step (dtStep =
// dtAna = 10, 0 .. 600 a.u.: 60 BK records, 3 grid points in the confidence interval; the
// values are not physical, the runs are only compared with each other):
//   full          0 .. 600 with alsoAnalysis
//   first/second  the same run split at 300 by a restart, alsoAnalysis in the second part
//   only_native   onlyAnalysis of full/bk_coefficients.bin
//   only_fortran  onlyAnalysis of full/fort.21
//   in_place      onlyAnalysis with the output directory = RestartDir (a copy of full)
struct FineAnalysisRuns {
  fs::path root = fs::temp_directory_path() / ("spodvr_fine_analysis_test_" + std::to_string(::getpid()));
  int code_full, code_first, code_second, code_native, code_fortran, code_in_place;
  std::string err;

  static nlohmann::json json(double t_start, double t_end) {
    nlohmann::json j = short_json();
    j["dtStep"] = 10.0;
    j["dtAna"] = 10.0;
    j["dtAve"] = 100.0;
    j["CheckpointFormat"] = "both";
    j["alsoAnalysis"] = true;
    return with_times(j, t_start, t_end, static_cast<int>((t_end - t_start) / 10.0 + 0.5));
  }
  static nlohmann::json only_analysis(const fs::path& dir, const std::string& format) {
    nlohmann::json j = json(0.0, 600.0);
    j["alsoAnalysis"] = false;
    j["onlyAnalysis"] = true;
    j["RestartDir"] = dir.string();
    j["RestartFormat"] = format;
    return j;
  }

  FineAnalysisRuns() {
    fs::remove_all(root);
    code_full = run(json(0.0, 600.0), root / "full", err);
    code_first = run(json(0.0, 300.0), root / "first", err);
    code_second = run(restart_from(json(300.0, 600.0), root / "first"), root / "second", err);
    code_native = run(only_analysis(root / "full", "native"), root / "only_native", err);
    code_fortran = run(only_analysis(root / "full", "fortran"), root / "only_fortran", err);
    fs::create_directories(root / "in_place");
    fs::copy_file(root / "full" / "bk_coefficients.bin", root / "in_place" / "bk_coefficients.bin");
    code_in_place = run(only_analysis(root / "in_place", "native"), root / "in_place", err);
  }
  ~FineAnalysisRuns() { fs::remove_all(root); }
};

const FineAnalysisRuns& fine_runs() {
  static const FineAnalysisRuns r;
  return r;
}

const std::vector<std::string> fine_analysis_files = {"fine_analysis/grid.json",
                                                      "fine_analysis/scattering_summary.tsv",
                                                      "fine_analysis/s_matrix.tsv",
                                                      "fine_analysis/probabilities_rovib.tsv",
                                                      "fine_analysis/probabilities_diffraction.tsv",
                                                      "fine_analysis/spectrum/scattering_summary.tsv",
                                                      "fine_analysis/spectrum/probabilities_rovib.tsv",
                                                      "fine_analysis/spectrum/probabilities_diffraction.tsv"};

}  // namespace

// alsoAnalysis is the same fine analysis as onlyAnalysis of the BK file the run wrote, from the
// native and from the Fortran file, and also after a restart.
TEST(FineAnalysisRuns, AllRunModesGiveBitwiseIdenticalResults) {
  const FineAnalysisRuns& r = fine_runs();
  for (const int code : {r.code_full, r.code_first, r.code_second, r.code_native, r.code_fortran, r.code_in_place})
    ASSERT_EQ(code, 0) << r.err;
  ASSERT_TRUE(fs::exists(r.root / "full" / "fine_analysis" / "s_matrix.tsv"));  // grid points in the interval
  for (const std::string& name : fine_analysis_files) {
    const std::string full = read_file(r.root / "full" / name);
    EXPECT_FALSE(full.empty()) << name;
    for (const std::string other : {"second", "only_native", "only_fortran", "in_place"})
      EXPECT_EQ(full, read_file(r.root / other / name)) << other << "/" << name;
  }
  // The analysis-only run in RestartDir leaves the BK file alone and writes no propagation files.
  EXPECT_EQ(read_file(r.root / "full" / "bk_coefficients.bin"), read_file(r.root / "in_place" / "bk_coefficients.bin"));
  EXPECT_FALSE(fs::exists(r.root / "in_place" / "checkpoint.bin"));
}

TEST(FineAnalysisRuns, Errors) {
  const FineAnalysisRuns& r = fine_runs();
  ASSERT_EQ(r.code_full, 0) << r.err;
  std::string err;
  // tEnd beyond the records in the file.
  nlohmann::json json = FineAnalysisRuns::only_analysis(r.root / "full", "native");
  json = with_times(json, 0.0, 700.0, 70);
  EXPECT_EQ(run(json, r.root / "too_long", err), 2);
  EXPECT_NE(err.find("tEnd"), std::string::npos) << err;
  // A BK file of a different input.
  json = FineAnalysisRuns::only_analysis(r.root / "full", "native");
  json["Z0"] = 11.6;
  EXPECT_EQ(run(json, r.root / "mismatch", err), 4);
  EXPECT_NE(err.find("Z0"), std::string::npos) << err;
  // No BK file: an input error, found before the set-up (RestartFiles below).
  json = FineAnalysisRuns::only_analysis(r.root / "missing", "fortran");
  EXPECT_EQ(run(json, r.root / "no_file", err), 2);
  EXPECT_NE(err.find("RestartDir"), std::string::npos) << err;
}

TEST(Restart, SplitRunIsBitwiseIdenticalToTheUnsplitRun) {
  const RestartRuns& r = runs();
  ASSERT_EQ(r.code_full, 0) << r.err;
  ASSERT_EQ(r.code_first, 0) << r.err;
  ASSERT_EQ(r.code_second, 0) << r.err;
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "scattering_summary.tsv",
                                 "probabilities_rovib.tsv", "probabilities_diffraction.tsv", "fort.9", "fort.21"})
    EXPECT_EQ(read_file(r.root / "full" / name), read_file(r.root / "second" / name)) << name;
  EXPECT_GT(fs::file_size(r.root / "second" / "bk_coefficients.bin"),
            fs::file_size(r.root / "first" / "bk_coefficients.bin"));
}

TEST(Restart, MismatchingInputIsAnError) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  nlohmann::json json = restart_from(with_times(short_json(), 60.0, 120.0, 30), r.root / "first");
  json["Z0"] = 11.6;
  std::string err;
  EXPECT_EQ(run(json, r.root / "mismatch", err), 4);
  EXPECT_NE(err.find("Z0"), std::string::npos) << err;

  json = restart_from(with_times(short_json(), 50.0, 120.0, 35), r.root / "first");
  EXPECT_EQ(run(json, r.root / "wrong_time", err), 4);
  EXPECT_NE(err.find("tStart"), std::string::npos) << err;
}

// UseSymmetry with nPhisOptimal may change the number of phi points (the propagation's
// quadrature): a restart must use the same number. It is checked itself, so a choice that keeps
// nphi is allowed (the hard wall's C2v needs an even nphi, which 10 already is), while the
// Dai-Light C6v case goes from 10 to 12 phi points.
TEST(Restart, NumberOfPhiPointsIsChecked) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  nlohmann::json json = restart_from(with_times(short_json(), 60.0, 120.0, 30), r.root / "first");
  json["UseSymmetry"] = true;
  json["nPhisOptimal"] = true;
  std::string err;
  EXPECT_EQ(run(json, r.root / "symmetric_restart", err), 0) << err;

  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_short/input.json");
  const nlohmann::json dai = with_times(nlohmann::json::parse(file), 0.0, 20.0, 10);
  ASSERT_EQ(run(dai, r.root / "dailight_first", err), 0) << err;
  nlohmann::json restart = restart_from(with_times(dai, 20.0, 40.0, 10), r.root / "dailight_first");
  restart["UseSymmetry"] = true;
  restart["nPhisOptimal"] = true;
  EXPECT_EQ(run(restart, r.root / "dailight_restart", err), 4);
  EXPECT_NE(err.find("number of phi points"), std::string::npos) << err;
}

TEST(Restart, TruncatedCheckpointIsAnError) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  const fs::path damaged = r.root / "damaged";
  fs::create_directories(damaged);
  fs::copy_file(r.root / "first" / "bk_coefficients.bin", damaged / "bk_coefficients.bin");
  const std::string checkpoint = read_file(r.root / "first" / "checkpoint.bin");
  std::ofstream(damaged / "checkpoint.bin", std::ios::binary) << checkpoint.substr(0, checkpoint.size() - 100);
  std::string err;
  EXPECT_EQ(run(restart_from(with_times(short_json(), 60.0, 120.0, 30), damaged), r.root / "after_damage", err), 4);
  EXPECT_NE(err.find("checkpoint"), std::string::npos) << err;
}

TEST(Restart, OutputDirectoryMustDifferFromRestartDir) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  const std::string before = read_file(r.root / "first" / "bk_coefficients.bin");
  std::string err;
  EXPECT_EQ(run(restart_from(with_times(short_json(), 60.0, 120.0, 30), r.root / "first"), r.root / "first", err), 2);
  EXPECT_NE(err.find("RestartDir"), std::string::npos) << err;
  EXPECT_EQ(read_file(r.root / "first" / "bk_coefficients.bin"), before);
}

TEST(Restart, RestartWithoutRestartDirIsAnInputError) {
  nlohmann::json json = with_times(short_json(), 60.0, 120.0, 30);
  json["NewCalc"] = false;
  std::string err;
  EXPECT_EQ(run(json, fs::temp_directory_path() / "spodvr_restart_no_dir", err), 2);
  EXPECT_NE(err.find("RestartDir"), std::string::npos) << err;
  fs::remove_all(fs::temp_directory_path() / "spodvr_restart_no_dir");
}

// A restart (run mode 1) or an analysis-only run (mode 2) whose RestartDir lacks a file it will
// read stops right after the input with an input error naming RestartDir and the file. Its stdout
// ends before the input summary ("Run mode"), so no set-up step ran (grids, PES, 6D potential);
// a valid restart prints both markers, so their absence is not vacuous. The files present are
// empty stand-ins: only their presence is checked this early.
TEST(RestartFiles, MissingFileStopsTheRunBeforeTheSetUp) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  std::string err;
  std::string out;
  const nlohmann::json restart = restart_from(with_times(short_json(), 60.0, 120.0, 30), r.root / "first");
  ASSERT_EQ(run(restart, r.root / "valid_restart", err, out), 0) << err;
  ASSERT_NE(out.find("Run mode"), std::string::npos) << out;
  ASSERT_NE(out.find("6D potential computed"), std::string::npos) << out;

  struct Case {
    std::string format;
    int run_mode;
    std::vector<std::string> present;
    std::string error;
  };
  const std::vector<Case> cases = {{"native", 1, {"bk_coefficients.bin"}, "no checkpoint.bin"},
                                   {"native", 1, {"checkpoint.bin"}, "no bk_coefficients.bin"},
                                   {"fortran", 1, {"fort.21"}, "neither fort.8 nor fort.9"},
                                   {"fortran", 1, {"fort.9"}, "neither fort.22 nor fort.21"},
                                   {"native", 2, {"checkpoint.bin", "fort.21"}, "no bk_coefficients.bin"},
                                   {"fortran", 2, {"bk_coefficients.bin"}, "neither fort.22 nor fort.21"}};
  const test_support::TemporaryDirectory dir("restart_files");
  for (std::size_t k = 0; k < cases.size(); ++k) {
    const Case& c = cases[k];
    const fs::path restart_dir = dir.path() / std::format("restart_dir_{}", k);
    fs::create_directories(restart_dir);
    for (const std::string& name : c.present) std::ofstream(restart_dir / name) << "";
    nlohmann::json json = c.run_mode == 1 ? restart_from(with_times(short_json(), 60.0, 120.0, 30), restart_dir)
                                          : FineAnalysisRuns::only_analysis(restart_dir, c.format);
    json["RestartFormat"] = c.format;
    EXPECT_EQ(run(json, dir.path() / std::format("output_{}", k), err, out), 2) << c.error << ": " << err;
    EXPECT_NE(err.find("RestartDir: " + c.error + " in '" + restart_dir.string() + "'"), std::string::npos) << err;
    EXPECT_EQ(out.find("Run mode"), std::string::npos) << c.error << ": " << out;
  }
}

// RestartDir must be a directory, and each file the run reads from it a readable regular file
// (here a directory named checkpoint.bin, and a directory fort.8 next to a regular fort.9: the
// name the reader picks is checked); both before the set-up.
TEST(RestartFiles, RestartDirAndItsFilesMustBeRegular) {
  const test_support::TemporaryDirectory dir("restart_dir");
  const nlohmann::json times = with_times(short_json(), 60.0, 120.0, 30);
  std::string err;
  std::string out;
  EXPECT_EQ(run(restart_from(times, dir.path() / "absent"), dir.path() / "output_1", err, out), 2) << err;
  EXPECT_NE(err.find("RestartDir: '" + (dir.path() / "absent").string() + "' is not a directory"), std::string::npos)
      << err;
  EXPECT_EQ(out.find("Run mode"), std::string::npos) << out;

  fs::create_directories(dir.path() / "restart" / "checkpoint.bin");
  std::ofstream(dir.path() / "restart" / "bk_coefficients.bin") << "";
  EXPECT_EQ(run(restart_from(times, dir.path() / "restart"), dir.path() / "output_2", err, out), 2) << err;
  EXPECT_NE(err.find("checkpoint.bin' (the wave packet) is not a readable regular file"), std::string::npos) << err;
  EXPECT_EQ(out.find("Run mode"), std::string::npos) << out;

  const fs::path fortran_dir = dir.path() / "fortran";
  fs::create_directories(fortran_dir / "fort.8");
  for (const std::string name : {"fort.9", "fort.21"}) std::ofstream(fortran_dir / name) << "";
  nlohmann::json fortran = restart_from(times, fortran_dir);
  fortran["RestartFormat"] = "fortran";
  EXPECT_EQ(run(fortran, dir.path() / "output_3", err, out), 2) << err;
  EXPECT_NE(err.find("fort.8' (the wave packet) is not a readable regular file"), std::string::npos) << err;
  EXPECT_EQ(out.find("Run mode"), std::string::npos) << out;
}

// Permissions: a file without read permission, and RestartDir without search permission
// (chmod 600: an InputError naming RestartDir, not a filesystem_error with exit code 1). Skipped
// where permissions do not restrict access (running as root).
TEST(RestartFiles, UnreadableFileOrRestartDirIsAnInputError) {
  const test_support::TemporaryDirectory dir("restart_permissions");
  const fs::path restart_dir = dir.path() / "restart";
  fs::create_directories(restart_dir);
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin"}) std::ofstream(restart_dir / name) << "";
  const fs::path checkpoint = restart_dir / "checkpoint.bin";
  fs::permissions(checkpoint, fs::perms::none);
  if (std::ifstream(checkpoint)) {
    fs::permissions(checkpoint, fs::perms::owner_read | fs::perms::owner_write);
    GTEST_SKIP() << "running as root: file permissions do not restrict access";
  }
  const nlohmann::json restart = restart_from(with_times(short_json(), 60.0, 120.0, 30), restart_dir);
  std::string err;
  std::string out;
  EXPECT_EQ(run(restart, dir.path() / "output_1", err, out), 2) << err;
  EXPECT_NE(err.find("checkpoint.bin' (the wave packet) is not a readable regular file"), std::string::npos) << err;
  EXPECT_EQ(out.find("Run mode"), std::string::npos) << out;
  fs::permissions(checkpoint, fs::perms::owner_read | fs::perms::owner_write);

  fs::permissions(restart_dir, fs::perms::owner_read | fs::perms::owner_write);
  EXPECT_EQ(run(restart, dir.path() / "output_2", err, out), 2) << err;
  EXPECT_NE(err.find("RestartDir: cannot access '" + checkpoint.string() + "': Permission denied"), std::string::npos)
      << err;
  EXPECT_EQ(out.find("Run mode"), std::string::npos) << out;
  fs::permissions(restart_dir, fs::perms::owner_all);  // so that the scratch directory can be removed
}

// Restarting from the Fortran-format files (fort.9 / fort.21) also continues the run exactly:
// the analysis surface state that the Fortran file lacks is reconstructed from the input.
TEST(Restart, RestartFromFortranFilesIsBitwiseIdenticalToTheUnsplitRun) {
  const RestartRuns& r = runs();
  ASSERT_EQ(r.code_second_fortran, 0) << r.err;
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "fort.9", "fort.21"})
    EXPECT_EQ(read_file(r.root / "full" / name), read_file(r.root / "second_fortran" / name)) << name;
}

namespace {

// Both low-Z absorbers on: the scattering one in the energy form, the specular one as a factor.
nlohmann::json with_low_z(nlohmann::json json) {
  json["ZOptLowStart"] = -2.0;
  json["ZOptLowEnd"] = 8.0;
  json["ZOptLowEnergy"] = 0.3;
  json["ZSpOptLowStart"] = -2.0;
  json["ZSpOptLowEnd"] = 6.0;
  json["ZSpOptLowFactor"] = 0.01;
  return json;
}

// With low-Z absorbers: the unsplit run (0 .. 120) and the split run (0 .. 60, 60 .. 120).
struct LowZRuns {
  fs::path root = fs::temp_directory_path() / ("spodvr_restart_low_z_" + std::to_string(::getpid()));
  int code_full;
  int code_first;
  int code_second;
  std::string err;

  LowZRuns() {
    fs::remove_all(root);
    const nlohmann::json json = with_low_z(short_json());
    code_full = run(with_times(json, 0.0, 120.0, 60), root / "full", err);
    code_first = run(with_times(json, 0.0, 60.0, 30), root / "first", err);
    code_second = run(restart_from(with_times(json, 60.0, 120.0, 30), root / "first"), root / "second", err);
  }
  ~LowZRuns() { fs::remove_all(root); }
};

// The low-Z runs, made once for every test that reads them.
const LowZRuns& low_z_runs() {
  static const LowZRuns r;
  return r;
}

}  // namespace

TEST(Restart, SplitRunWithLowZAbsorbersIsBitwiseIdenticalToTheUnsplitRun) {
  const LowZRuns& r = low_z_runs();
  ASSERT_EQ(r.code_full, 0) << r.err;
  ASSERT_EQ(r.code_first, 0) << r.err;
  ASSERT_EQ(r.code_second, 0) << r.err;
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "probabilities_rovib.tsv"})
    EXPECT_EQ(read_file(r.root / "full" / name), read_file(r.root / "second" / name)) << name;
}

// A restart or an analysis-only run must have the low-Z absorbers of the run that wrote the file:
// one switched on or off only in the input is named, and so is a changed low-Z parameter.
TEST(Restart, LowZAbsorberMismatchIsAnError) {
  const LowZRuns& r = low_z_runs();
  ASSERT_EQ(r.code_first, 0) << r.err;
  std::string err;
  const nlohmann::json second = with_times(short_json(), 60.0, 120.0, 30);
  EXPECT_EQ(run(restart_from(second, r.root / "first"), r.root / "without", err), 4);
  EXPECT_NE(err.find("the low-Z absorber on the scattering grid (ZOptLowStart, ZOptLowEnd) is off in the input but on "
                     "in checkpoint"),
            std::string::npos)
      << err;

  const FirstRun& plain = first_run();
  ASSERT_EQ(plain.code, 0) << plain.err;
  EXPECT_EQ(run(restart_from(with_low_z(second), plain.root / "first"), r.root / "with", err), 4);
  EXPECT_NE(err.find("is on in the input but off in checkpoint"), std::string::npos) << err;

  nlohmann::json changed = restart_from(with_low_z(second), r.root / "first");
  changed["ZSpOptLowEnd"] = 5.92;
  EXPECT_EQ(run(changed, r.root / "changed_end", err), 4);
  EXPECT_NE(err.find("ZSpOptLowEnd is"), std::string::npos) << err;
  changed = restart_from(with_low_z(second), r.root / "first");
  changed["ZOptLowEnergy"] = 0.25;
  EXPECT_EQ(run(changed, r.root / "changed_energy", err), 4);
  EXPECT_NE(err.find("ZOptLowFactor (effective) is"), std::string::npos) << err;

  // Mode 2 (onlyAnalysis) reads the same list from the BK file.
  nlohmann::json only = with_times(short_json(), 0.0, 120.0, 60);
  only["onlyAnalysis"] = true;
  only["RestartDir"] = (r.root / "full").string();
  EXPECT_EQ(run(only, r.root / "only_without", err), 4);
  EXPECT_NE(err.find("the low-Z absorber on the scattering grid (ZOptLowStart, ZOptLowEnd) is off in the input but on "
                     "in BK coefficient file"),
            std::string::npos)
      << err;
}

// Files written without low-Z absorbers and without SCM keep the parameter list of format version 2
// (35 entries, no low-Z or SCM names), so checkpoints and BK files of earlier versions of the program
// stay readable.
TEST(Restart, CheckpointWithoutLowZAbsorbersKeepsItsParameterList) {
  const FirstRun& r = first_run();
  ASSERT_EQ(r.code, 0) << r.err;
  const std::vector<std::string> format_2 = {"ZStart",
                                             "dZ",
                                             "nZ",
                                             "nZSpecular",
                                             "RStart",
                                             "dR",
                                             "nR",
                                             "nX",
                                             "nY",
                                             "JMax",
                                             "MJMax",
                                             "number of phi points (MJMax, nPhisOptimal)",
                                             "Mass1",
                                             "Mass2",
                                             "HomoNuc",
                                             "dtStep",
                                             "dtAna",
                                             "Zana",
                                             "Z0",
                                             "EMinWP",
                                             "EMaxWP",
                                             "EParallel",
                                             "PhiParallel",
                                             "VIn",
                                             "JIn",
                                             "MJIn",
                                             "ZOptStart",
                                             "ZOptEnd",
                                             "ROptStart",
                                             "ROptEnd",
                                             "ZSpOptStart",
                                             "ZSpOptEnd",
                                             "ZOptFactor (effective)",
                                             "ROptFactor (effective)",
                                             "ZSpOptFactor (effective)"};
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin"}) {
    BinaryReader file(r.root / "first" / name);
    (void)file.read<std::array<char, 16>>();  // magic
    EXPECT_EQ(file.read<std::uint32_t>(), 2u) << name;
    EXPECT_EQ(file.read_string(), "HardWallPES") << name;
    ASSERT_EQ(file.read<std::uint64_t>(), format_2.size()) << name;
    for (const std::string& expected : format_2) {
      EXPECT_EQ(file.read_string(), expected) << name;
      (void)file.read<double>();
    }
  }
}

namespace {

// The SCM example (SystemID CRPDistrtSurface) on 2 × 2 (x, y) points, from a position file with two
// identical snapshots: the first part of a run (0 .. 10, 4 steps) and the input of the second part.
struct ScmRuns {
  test_support::TemporaryDirectory dir{"restart_scm"};
  fs::path two_snapshots = dir.path() / "two_snapshots.dat";
  int code_first;
  std::string err;

  // The example position file with its snapshot twice (line 2: 2 snapshots).
  void write_two_snapshots() const {
    std::ifstream in(test_support::scm_example_directory / "examplePosfile.dat");
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) lines.push_back(line);
    ASSERT_EQ(lines.size(), 3u + 343u);
    std::ofstream out(two_snapshots);
    out << lines[0] << "\n 7 7 7 2 3 343\n" << lines[2] << '\n';
    for (int copy = 0; copy < 2; ++copy)
      for (std::size_t k = 3; k < lines.size(); ++k) out << lines[k] << '\n';
  }
  nlohmann::json json(double t_start, double t_end) const {
    nlohmann::json j = test_support::scm_example_json(two_snapshots);
    j.update({{"nX", 2}, {"nY", 2}, {"dtAna", 5.0}, {"dtAve", 5.0}});
    return with_times(j, t_start, t_end, static_cast<int>((t_end - t_start) / 2.5 + 0.5));
  }
  nlohmann::json second(const nlohmann::json& changes = nlohmann::json::object()) const {
    nlohmann::json j = restart_from(json(10.0, 20.0), dir.path() / "first");
    j.update(changes);
    return j;
  }

  ScmRuns() {
    write_two_snapshots();
    code_first = run(json(0.0, 10.0), dir.path() / "first", err);
  }
};

// The SCM runs, made once for every test that reads them.
const ScmRuns& scm_runs() {
  static const ScmRuns r;
  return r;
}

}  // namespace

// A restart of an SCM run with the same input continues it; one with another snapshot
// (PosFileSlabNum), cut-off, rescaling (EnableSurfExpan), parameter file (a changed P7 of layer 1) or
// position file (one atom moved by 0.01 bohr: the same a_T and atom count) is refused naming it (exit
// code 4). A Fortran-format restart has no parameters to check.
TEST(Restart, ScmIdentityMismatchIsAnError) {
  const ScmRuns& r = scm_runs();
  ASSERT_EQ(r.code_first, 0) << r.err;
  std::string err;
  EXPECT_EQ(run(r.second(), r.dir.path() / "same", err), 0) << err;
  EXPECT_EQ(run(r.second({{"PosFileSlabNum", 2}}), r.dir.path() / "slab_2", err), 4);
  EXPECT_NE(err.find("PosFileSlabNum is 2.0000000000000000e+00 in the input but 1.0000000000000000e+00 in checkpoint"),
            std::string::npos)
      << err;
  EXPECT_EQ(run(r.second({{"SCMCutoff", 15.0}}), r.dir.path() / "cutoff", err), 4);
  EXPECT_NE(err.find("SCMCutoff is"), std::string::npos) << err;
  EXPECT_EQ(run(r.second({{"EnableSurfExpan", false}}), r.dir.path() / "no_rescaling", err), 4);
  EXPECT_NE(err.find("EnableSurfExpan is 0.0000000000000000e+00 in the input"), std::string::npos) << err;

  std::ifstream original(test_support::scm_example_directory / "params.dat");
  const std::string text{std::istreambuf_iterator<char>(original), std::istreambuf_iterator<char>()};
  const std::string p7 = "7.44024869632";  // P7 of layer 1, the first of the two lines
  ASSERT_NE(text.find(p7), std::string::npos);
  std::string changed = text;
  changed.replace(text.find(p7), p7.size(), "7.44024869633");
  std::ofstream(r.dir.path() / "params.dat") << changed;
  EXPECT_EQ(run(r.second({{"SCMParFilepath", (r.dir.path() / "params.dat").string()}}), r.dir.path() / "params", err),
            4);
  EXPECT_NE(err.find("SCMParFilepath (layer 1, P7) is 7.4402486963300"), std::string::npos) << err;

  // The atom at the origin (ideal (0, 0, 0)) of snapshot 1 moved by 0.01 bohr in x.
  std::ifstream two(r.two_snapshots);
  std::vector<std::string> lines;
  for (std::string line; std::getline(two, line);) lines.push_back(line);
  for (std::size_t k = 3; k < lines.size(); ++k) {
    std::istringstream fields(lines[k]);
    double x, y, z, xi, yi, zi;
    std::string rest;
    if (!(fields >> x >> y >> z >> xi >> yi >> zi) || xi != 0.0 || yi != 0.0 || zi != 0.0) continue;
    std::getline(fields, rest);
    lines[k] = std::format("{} {} {} {} {} {}{}", x + 0.01, y, z, xi, yi, zi, rest);
    break;
  }
  std::ofstream moved(r.dir.path() / "moved.dat");
  for (const std::string& line : lines) moved << line << '\n';
  moved.close();
  EXPECT_EQ(run(r.second({{"PosFilepath", (r.dir.path() / "moved.dat").string()}}), r.dir.path() / "moved", err), 4);
  EXPECT_NE(err.find("PosFilepath (rms displacement) is"), std::string::npos) << err;
}

// The header of an SCM checkpoint: the 35 parameters of every file, then the SCM identity (snapshot,
// cut-off, rescaling, the slab's a_T, atom count and position statistics, P1–P18 of both layers).
TEST(Restart, ScmCheckpointHasTheScmIdentity) {
  const ScmRuns& r = scm_runs();
  ASSERT_EQ(r.code_first, 0) << r.err;
  BinaryReader file(r.dir.path() / "first" / "checkpoint.bin");
  (void)file.read<std::array<char, 16>>();  // magic
  EXPECT_EQ(file.read<std::uint32_t>(), 2u);
  EXPECT_EQ(file.read_string(), "CRPDistrtSurface");
  ASSERT_EQ(file.read<std::uint64_t>(), 35u + 10u + 36u);
  for (int k = 0; k < 35; ++k) {
    (void)file.read_string();
    (void)file.read<double>();
  }
  const std::pair<std::string, double> identity[] = {{"PosFileSlabNum", 1.0},
                                                     {"SCMCutoff", 16.0},
                                                     {"EnableSurfExpan", 1.0},
                                                     {"PosFilepath (a_T)", 4.8896663475},
                                                     {"PosFilepath (atoms used)", 73.0}};
  for (const auto& [name, value] : identity) {
    EXPECT_EQ(file.read_string(), name);
    EXPECT_EQ(file.read<double>(), value) << name;
  }
  // The position statistics: the values of the run's identity, bitwise.
  const ScmIdentity scm = *build_config(parse_input(r.json(0.0, 10.0))).derived.scm;
  const std::pair<std::string_view, double> statistics[] = {
      {"PosFilepath (rms displacement)", scm.rms_displacement},
      {"PosFilepath (layer 1 mean shift x)", scm.top_layer_shift[0]},
      {"PosFilepath (layer 1 mean shift y)", scm.top_layer_shift[1]},
      {"PosFilepath (layer 1 mean shift z)", scm.top_layer_shift[2]},
      {"PosFilepath (highest atom)", scm.highest_atom}};
  for (const auto& [name, value] : statistics) {
    EXPECT_EQ(file.read_string(), name);
    EXPECT_EQ(file.read<double>(), value) << name;
  }
  for (int layer = 1; layer <= 2; ++layer)
    for (int k = 1; k <= 18; ++k) {
      EXPECT_EQ(file.read_string(), "SCMParFilepath (layer " + std::to_string(layer) + ", P" + std::to_string(k) + ")");
      const double value = file.read<double>();
      if (k == 17) {
        EXPECT_EQ(value, 0.725);  // the r clamp range of params.dat
      }
      if (k == 18) {
        EXPECT_EQ(value, 1.847);
      }
    }
}

// analysis.json of an SCM run holds its whole identity, bitwise the values of the run's configuration;
// that of another PES has no "scm" entry.
TEST(AnalysisFiles, ScmIdentityRoundTrip) {
  const ScmRuns& r = scm_runs();
  ASSERT_EQ(r.code_first, 0) << r.err;
  std::ifstream file(r.dir.path() / "first" / "analysis.json");
  const nlohmann::json json = nlohmann::json::parse(file);
  ASSERT_TRUE(json.contains("scm"));
  const nlohmann::json& scm = json["scm"];
  const Config config = build_config(parse_input(r.json(0.0, 10.0)));
  ASSERT_TRUE(config.derived.scm.has_value());
  const ScmIdentity& identity = *config.derived.scm;
  EXPECT_EQ(scm["crp_data_set"], config.input.crp_dataset);
  EXPECT_EQ(scm["parameter_file"], config.input.scm.parameter_file);
  EXPECT_EQ(scm["position_file"], r.two_snapshots.string());
  EXPECT_EQ(scm["slab_tag"], "Multiseedv3_New");
  EXPECT_EQ(scm["slab_tag"], identity.slab_tag);
  EXPECT_EQ(scm["position_file_version"], 3);
  EXPECT_EQ(scm["slab_number"], 1);
  EXPECT_EQ(scm["cutoff"].get<double>(), 16.0);
  EXPECT_EQ(scm["rescaling"], true);
  EXPECT_EQ(scm["slab_lattice_constant"].get<double>(), identity.slab_lattice_constant);
  EXPECT_EQ(scm["slab_atoms"], 73);
  EXPECT_EQ(scm["rms_displacement"].get<double>(), identity.rms_displacement);
  ASSERT_EQ(scm["top_layer_shift"].size(), 3u);
  for (int k = 0; k < 3; ++k) EXPECT_EQ(scm["top_layer_shift"][k].get<double>(), identity.top_layer_shift[k]);
  EXPECT_EQ(scm["highest_atom"].get<double>(), identity.highest_atom);
  EXPECT_NEAR(identity.highest_atom, 0.936367, 1e-6);  // the report's value
  ASSERT_EQ(scm["parameters"].size(), 2u);
  for (int layer = 0; layer < 2; ++layer) {
    ASSERT_EQ(scm["parameters"][layer].size(), 18u);
    for (int k = 0; k < 18; ++k)
      EXPECT_EQ(scm["parameters"][layer][k].get<double>(), identity.parameters[layer][k]) << layer << ", " << k;
  }

  const FirstRun& plain = first_run();
  ASSERT_EQ(plain.code, 0) << plain.err;
  std::ifstream plain_file(plain.root / "first" / "analysis.json");
  EXPECT_FALSE(nlohmann::json::parse(plain_file).contains("scm"));
}

// A file written by gfortran with -fmax-subrecord-length=16 (tests/reference/fortran_unformatted).
TEST(FortranUnformatted, ReadsGfortranRecordsAndSubrecords) {
  const fs::path path = fs::path(SPODVR_REFERENCE_DIR) / "fortran_unformatted" / "subrecords.bin";
  FortranRecordReader file(path);
  FortranRecordParser first = file.read("record 1");
  EXPECT_EQ(first.read<std::int32_t>(), 7);
  first.expect_end();
  FortranRecordParser second = file.read("record 2");
  EXPECT_EQ(second.read<std::int32_t>(), 1);  // .TRUE.
  EXPECT_EQ(second.read<std::int32_t>(), 42);
  FortranRecordParser third = file.read("record 3");
  std::vector<double> a(5);
  third.read_array(std::span<double>(a));
  third.expect_end();
  EXPECT_EQ(a, (std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0}));
  EXPECT_TRUE(file.at_end());
}

TEST(FortranUnformatted, WritesTheSameBytesAsGfortran) {
  const fs::path reference = fs::path(SPODVR_REFERENCE_DIR) / "fortran_unformatted" / "subrecords.bin";
  const fs::path path = fs::temp_directory_path() / ("spodvr_records_" + std::to_string(::getpid()) + ".bin");
  {
    FortranRecordWriter file(path, 16);
    file.write(FortranRecord().add<std::int32_t>(7));
    file.write(FortranRecord().add<std::int32_t>(1).add<std::int32_t>(42));
    const std::vector<double> a = {1.0, 2.0, 3.0, 4.0, 5.0};
    file.write_array(std::span<const double>(a));
  }
  EXPECT_EQ(read_file(path), read_file(reference));
  fs::remove(path);
}

TEST(FortranUnformatted, BadMarkersAndShortRecordsAreErrors) {
  const fs::path path = fs::temp_directory_path() / ("spodvr_bad_records_" + std::to_string(::getpid()) + ".bin");
  {
    BinaryWriter file(path);  // head marker 8, 4 data bytes, tail marker 7
    for (const std::int32_t value : {8, 5, 7}) file.write(value);
  }
  FortranRecordReader file(path);
  EXPECT_THROW((void)file.read("record 1"), CheckpointError);
  fs::remove(path);

  FortranRecordParser record(std::vector<std::byte>(4), "test record");
  EXPECT_THROW((void)record.read<double>(), CheckpointError);
}
