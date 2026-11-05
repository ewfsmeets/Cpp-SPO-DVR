// Tests: checkpoint / restart. A run split by a restart must give bitwise the same results as
// the unsplit run; damaged or mismatching checkpoints are clear errors.
#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "io/fortran_unformatted.hpp"

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
  static nlohmann::json only_analysis(const fs::path& dir, const char* format) {
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
    for (const char* other : {"second", "only_native", "only_fortran", "in_place"})
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
  // No BK file.
  json = FineAnalysisRuns::only_analysis(r.root / "missing", "fortran");
  EXPECT_EQ(run(json, r.root / "no_file", err), 4);
}

TEST(Restart, SplitRunIsBitwiseIdenticalToTheUnsplitRun) {
  const RestartRuns& r = runs();
  ASSERT_EQ(r.code_full, 0) << r.err;
  ASSERT_EQ(r.code_first, 0) << r.err;
  ASSERT_EQ(r.code_second, 0) << r.err;
  for (const char* name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "scattering_summary.tsv",
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

// Restarting from the Fortran-format files (fort.9 / fort.21) also continues the run exactly:
// the analysis surface state that the Fortran file lacks is reconstructed from the input.
TEST(Restart, RestartFromFortranFilesIsBitwiseIdenticalToTheUnsplitRun) {
  const RestartRuns& r = runs();
  ASSERT_EQ(r.code_second_fortran, 0) << r.err;
  for (const char* name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "fort.9", "fort.21"})
    EXPECT_EQ(read_file(r.root / "full" / name), read_file(r.root / "second_fortran" / name)) << name;
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
