// Tests: non-finite results (exit code 5). The finder (analysis/non_finite.hpp) per group of
// quantities and its scan order, the wording of the error line (report/analysis_report.hpp), and
// two runs: a restart from a checkpoint whose wave packet holds a NaN, and an analysis-only run of
// a BK file whose last record holds one. Both write every file of a clean run and exit with 5.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/non_finite.hpp"
#include "app/cli.hpp"
#include "app/fine_analysis_run.hpp"
#include "app/run.hpp"
#include "fft/fftw_fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "report/analysis_report.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
namespace fs = std::filesystem;
using cplx = std::complex<double>;

namespace {

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();
constexpr double infinity = std::numeric_limits<double>::infinity();

// One energy's results with 2 states x 2 levels, 3 rotations, 2 orders, every number finite.
// Some values are odd but legitimate: P > 1, a negative reaction probability, S = -0, a zero phase.
EnergyResult finite_energy() {
  EnergyResult e{};
  e.total_energy = 0.5;
  e.collision_energy = 0.3;
  e.z_energy = 0.3;
  e.total_momentum = 7.0;
  e.s_matrix = {cplx(-0.0, 0.0), cplx(0.6, -0.2), cplx(1.1, 0.0), cplx(0.0, -0.0)};
  e.probability = {0.0, 0.4, 1.21, 0.0};
  e.phase = {0.0, -0.32, 0.0, 0.0};
  e.open = {1, 1, 1, 0};
  e.p_vibration = {0.4, 1.21};
  e.p_rotation = {0.4, 0.0, 0.0, 1.21, 0.0, 0.0};
  e.p_order = std::vector<double>(12, 0.1);
  e.e_rotation_v = {0.01, 0.0};
  e.e_diffraction_v = {0.0, 0.02};
  e.p_total = 1.61;
  e.p_reaction = -0.61;
  e.e_vibration = 0.1;
  e.e_rotation = 0.01;
  e.e_diffraction = 0.02;
  return e;
}

ScatteringResults finite_results(int n_energies) {
  ScatteringResults r;
  for (int i = 0; i < n_energies; ++i) r.energies.push_back(finite_energy());
  return r;
}

void expect_found(const ScatteringResults& r, int energy_index, ResultQuantity quantity, int index) {
  const std::optional<NonFiniteValue> found = find_non_finite(r);
  ASSERT_TRUE(found.has_value());
  EXPECT_EQ(found->energy_index, energy_index);
  EXPECT_EQ(found->quantity, quantity);
  EXPECT_EQ(found->index, index);
}

}  // namespace

TEST(NonFiniteResults, FiniteResultsHaveNone) {
  EXPECT_FALSE(find_non_finite(ScatteringResults{}).has_value());
  EXPECT_FALSE(find_non_finite(finite_results(3)).has_value());
  ScatteringResults empty_lists = finite_results(1);  // an energy without states
  EnergyResult& e = empty_lists.energies[0];
  e.s_matrix.clear();
  e.probability.clear();
  e.phase.clear();
  e.open.clear();
  e.p_vibration.clear();
  e.p_rotation.clear();
  e.p_order.clear();
  e.e_rotation_v.clear();
  e.e_diffraction_v.clear();
  EXPECT_FALSE(find_non_finite(empty_lists).has_value());
  // Extreme but finite numbers.
  ScatteringResults extreme = finite_results(1);
  extreme.energies[0].p_total = std::numeric_limits<double>::max();
  extreme.energies[0].p_reaction = -std::numeric_limits<double>::max();
  extreme.energies[0].probability[0] = std::numeric_limits<double>::denorm_min();
  EXPECT_FALSE(find_non_finite(extreme).has_value());
}

// Every quantity poisoned at once: the finder reports them one after the other in the documented
// order (removing each as it is found), energies before quantities.
TEST(NonFiniteResults, FirstNonFiniteFieldIsNamed) {
  using enum ResultQuantity;
  ScatteringResults r = finite_results(3);
  EnergyResult& e = r.energies[1];
  e.total_energy = not_a_number;
  e.collision_energy = infinity;
  e.z_energy = -infinity;
  e.total_momentum = not_a_number;
  e.s_matrix[2] = cplx(not_a_number, 0.0);
  e.probability[3] = not_a_number;
  e.phase[1] = not_a_number;
  e.p_vibration[1] = infinity;
  e.p_rotation[4] = not_a_number;
  e.p_order[7] = not_a_number;
  e.e_rotation_v[0] = not_a_number;
  e.e_diffraction_v[1] = -infinity;
  e.p_total = not_a_number;
  e.p_reaction = not_a_number;
  e.e_vibration = infinity;
  e.e_rotation = not_a_number;
  e.e_diffraction = not_a_number;
  r.energies[2].total_energy = not_a_number;  // a later energy comes last

  struct Expected {
    ResultQuantity quantity;
    int index;
  };
  const std::vector<Expected> order = {{total_energy, 0},   {collision_energy, 0}, {z_energy, 0},
                                       {total_momentum, 0}, {s_matrix, 2},         {probability, 3},
                                       {phase, 1},          {p_vibration, 1},      {p_rotation, 4},
                                       {p_order, 7},        {e_rotation_v, 0},     {e_diffraction_v, 1},
                                       {p_total, 0},        {p_reaction, 0},       {e_vibration, 0},
                                       {e_rotation, 0},     {e_diffraction, 0}};
  const EnergyResult clean = finite_energy();
  for (const Expected& x : order) {
    const std::optional<NonFiniteValue> found = find_non_finite(r);
    ASSERT_TRUE(found.has_value());
    EXPECT_EQ(found->energy_index, 1);
    EXPECT_EQ(found->quantity, x.quantity);
    EXPECT_EQ(found->index, x.index);
    // Remove what was found: the whole energy is replaced field by field.
    switch (x.quantity) {
      case total_energy: e.total_energy = clean.total_energy; break;
      case collision_energy: e.collision_energy = clean.collision_energy; break;
      case z_energy: e.z_energy = clean.z_energy; break;
      case total_momentum: e.total_momentum = clean.total_momentum; break;
      case s_matrix: e.s_matrix = clean.s_matrix; break;
      case probability: e.probability = clean.probability; break;
      case phase: e.phase = clean.phase; break;
      case p_vibration: e.p_vibration = clean.p_vibration; break;
      case p_rotation: e.p_rotation = clean.p_rotation; break;
      case p_order: e.p_order = clean.p_order; break;
      case e_rotation_v: e.e_rotation_v = clean.e_rotation_v; break;
      case e_diffraction_v: e.e_diffraction_v = clean.e_diffraction_v; break;
      case p_total: e.p_total = clean.p_total; break;
      case p_reaction: e.p_reaction = clean.p_reaction; break;
      case e_vibration: e.e_vibration = clean.e_vibration; break;
      case e_rotation: e.e_rotation = clean.e_rotation; break;
      case e_diffraction: e.e_diffraction = clean.e_diffraction; break;
    }
  }
  expect_found(r, 2, total_energy, 0);
  r.energies[2].total_energy = 1.0;
  EXPECT_FALSE(find_non_finite(r).has_value());
}

TEST(NonFiniteResults, EnergiesAndMomentum) {
  for (const double bad : {not_a_number, infinity, -infinity}) {
    ScatteringResults r = finite_results(2);
    r.energies[1].z_energy = bad;
    expect_found(r, 1, ResultQuantity::z_energy, 0);
    r.energies[0].total_momentum = bad;
    expect_found(r, 0, ResultQuantity::total_momentum, 0);
  }
}

// S counts when either part is not finite; the first state in (is, v) order is named.
TEST(NonFiniteResults, SMatrixEitherPart) {
  ScatteringResults r = finite_results(1);
  r.energies[0].s_matrix[3] = cplx(0.0, infinity);
  expect_found(r, 0, ResultQuantity::s_matrix, 3);
  r.energies[0].s_matrix[1] = cplx(not_a_number, 0.5);
  expect_found(r, 0, ResultQuantity::s_matrix, 1);
}

TEST(NonFiniteResults, ProbabilitiesAndPhases) {
  ScatteringResults r = finite_results(1);
  r.energies[0].phase[0] = not_a_number;
  expect_found(r, 0, ResultQuantity::phase, 0);
  r.energies[0].probability[2] = infinity;
  expect_found(r, 0, ResultQuantity::probability, 2);
}

TEST(NonFiniteResults, SummedProbabilities) {
  ScatteringResults r = finite_results(1);
  r.energies[0].p_total = not_a_number;
  expect_found(r, 0, ResultQuantity::p_total, 0);
  r.energies[0].p_reaction = not_a_number;  // after p_total
  expect_found(r, 0, ResultQuantity::p_total, 0);
  r.energies[0].p_order[11] = infinity;
  expect_found(r, 0, ResultQuantity::p_order, 11);
  r.energies[0].p_rotation[5] = not_a_number;
  expect_found(r, 0, ResultQuantity::p_rotation, 5);
  r.energies[0].p_vibration[0] = not_a_number;
  expect_found(r, 0, ResultQuantity::p_vibration, 0);
}

TEST(NonFiniteResults, MeanEnergies) {
  ScatteringResults r = finite_results(1);
  r.energies[0].e_diffraction = -infinity;
  expect_found(r, 0, ResultQuantity::e_diffraction, 0);
  r.energies[0].e_vibration = not_a_number;
  expect_found(r, 0, ResultQuantity::e_vibration, 0);
  r.energies[0].e_diffraction_v[1] = not_a_number;
  expect_found(r, 0, ResultQuantity::e_diffraction_v, 1);
  r.energies[0].e_rotation_v[1] = not_a_number;
  expect_found(r, 0, ResultQuantity::e_rotation_v, 1);
}

namespace {

nlohmann::json short_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_short/input.json");
  return nlohmann::json::parse(file);
}

// The analysis basis of hardwall_short (no 6D potential: the gas phase needs only V(r)).
struct ShortBasis {
  std::unique_ptr<Potential> pes = make_potential(parse_input(short_json()));
  Config config = build_config(parse_input(short_json()), *pes);
  AngularGrid angular = make_angular_grid(config);
  SpatialGrids grids = make_spatial_grids(config);
  FftwFft fft;
  GasPhase gas = compute_gas_phase(config, grids, angular, compute_asymptotic_potential(*pes, config, grids), fft,
                                   LapackEigensolver());
  AnalysisSetup setup{config, angular, gas};
};

// Results of the shape of `basis` (every number 0.25), at two energies.
ScatteringResults results_of(const AnalysisBasis& basis) {
  const int n_states = basis.n_states();
  const int n_v = basis.n_v();
  EnergyResult e{};
  e.collision_energy = 0.25;
  e.s_matrix.assign(n_states * n_v, cplx(0.25, 0.0));
  e.probability.assign(n_states * n_v, 0.25);
  e.phase.assign(n_states * n_v, 0.25);
  e.p_vibration.assign(n_v, 0.25);
  e.p_rotation.assign(basis.n_rotation() * n_v, 0.25);
  e.p_order.assign((basis.max_order() + 1) * basis.n_rotation() * n_v, 0.25);
  e.e_rotation_v.assign(n_v, 0.25);
  e.e_diffraction_v.assign(n_v, 0.25);
  ScatteringResults r{{e, e}};
  r.energies[1].collision_energy = 0.5;
  return r;
}

std::string labels(const AnalysisBasis& basis, int is, int v) {
  const DiffractionChannel& c = basis.channel(is);
  return "(v, j, mJ, n, m) = (" + std::to_string(v) + ", " + std::to_string(basis.j(is)) + ", " +
         std::to_string(basis.mj(is)) + ", " + std::to_string(c.n) + ", " + std::to_string(c.m) + ")";
}

}  // namespace

// The error line names the result set, the energy and the state labels decoded from the flat index
// (layouts of EnergyResult), on the basis of hardwall_short.
TEST(NonFiniteResults, WordingNamesTheStateLabels) {
  const ShortBasis b;
  const AnalysisBasis& basis = b.setup.basis;
  ASSERT_GT(basis.n_states(), 1);
  ASSERT_GT(basis.n_rotation(), 1);
  ASSERT_GT(basis.max_order(), 0);
  const int n_states = basis.n_states();

  ScatteringResults r = results_of(basis);
  const int k = (n_states - 1) + n_states * (basis.n_v() - 1);  // the last state of the last level
  r.energies[1].s_matrix[k] = cplx(0.5, not_a_number);
  std::optional<NonFiniteValue> found = find_non_finite(r);
  ASSERT_TRUE(found.has_value());
  EXPECT_EQ(describe_non_finite("the on-the-fly analysis", *found, basis, r),
            "the on-the-fly analysis, energy 2 of 2 (E_coll = 0.500000 eV): S of " +
                labels(basis, n_states - 1, basis.n_v() - 1) + " is (5.0000000000000000e-01, nan)");

  // P summed per diffraction order: (order, irotana, v) = order + (mdifa+1) * (irotana + nana * v).
  r = results_of(basis);
  const int order = basis.max_order();
  const int irotana = basis.n_rotation() - 1;
  const int v = basis.n_v() - 1;
  r.energies[0].p_order[order + (basis.max_order() + 1) * (irotana + basis.n_rotation() * v)] = infinity;
  found = find_non_finite(r);
  ASSERT_TRUE(found.has_value());
  const int is = irotana * basis.n_diffraction();
  EXPECT_EQ(describe_non_finite("the fine analysis (spectrum)", *found, basis, r),
            "the fine analysis (spectrum), energy 1 of 2 (E_coll = 0.250000 eV): P of (v, j, mJ) = (" +
                std::to_string(v) + ", " + std::to_string(basis.j(is)) + ", " + std::to_string(basis.mj(is)) +
                "), diffraction order " + std::to_string(order) + " is inf");

  // P summed over diffraction: (irotana, v) = irotana + nana * v.
  r = results_of(basis);
  r.energies[0].p_rotation[irotana + basis.n_rotation() * v] = not_a_number;
  found = find_non_finite(r);
  ASSERT_TRUE(found.has_value());
  EXPECT_EQ(describe_non_finite("x", *found, basis, r), "x, energy 1 of 2 (E_coll = 0.250000 eV): P of (v, j, mJ) = (" +
                                                            std::to_string(v) + ", " + std::to_string(basis.j(is)) +
                                                            ", " + std::to_string(basis.mj(is)) + ") is nan");

  r = results_of(basis);
  r.energies[0].p_reaction = -infinity;
  found = find_non_finite(r);
  ASSERT_TRUE(found.has_value());
  EXPECT_EQ(describe_non_finite("x", *found, basis, r), "x, energy 1 of 2 (E_coll = 0.250000 eV): P_reaction is -inf");
}

namespace {

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
  std::ostringstream output;
  std::ostringstream errors;
  const int code = run_program(cl, output, errors, false);
  err = errors.str();
  return code;
}

nlohmann::json with_times(nlohmann::json json, double t_start, double t_end, int n_steps) {
  json["tStart"] = t_start;
  json["tEnd"] = t_end;
  json["nStep"] = n_steps;
  return json;
}

// The files under dir, as paths relative to dir.
std::set<std::string> file_set(const fs::path& dir) {
  std::set<std::string> files;
  for (const fs::directory_entry& entry : fs::recursive_directory_iterator(dir))
    files.insert(fs::relative(entry.path(), dir).string());
  return files;
}

// Overwrites the double at byte `offset` of the file with a NaN.
void patch_nan(const fs::path& path, std::size_t offset) {
  std::string bytes = read_file(path);
  ASSERT_LE(offset + sizeof(double), bytes.size());
  const auto nan_bytes = std::bit_cast<std::array<char, sizeof(double)>>(not_a_number);
  std::copy(nan_bytes.begin(), nan_bytes.end(), bytes.begin() + static_cast<std::ptrdiff_t>(offset));
  std::ofstream(path, std::ios::binary) << bytes;
}

int count_lines_with(const std::string& text, const std::string& what) {
  int count = 0;
  std::istringstream lines(text);
  for (std::string line; std::getline(lines, line);)
    if (line.find(what) != std::string::npos) count = count + 1;
  return count;
}

}  // namespace

// A restart from a checkpoint whose wave packet holds one NaN (the last element of the scattering wave,
// just before the 8-byte end marker): the NaN spreads through the FFTs into the BK coefficients. The run
// writes the files of the clean restart, prints one "Non-finite result" line and exits with 5.
TEST(NonFiniteResults, RunWithANaNExitsWith5AndWritesItsFiles) {
  const test_support::TemporaryDirectory dir("non_finite_run");
  const nlohmann::json json = short_json();
  std::string err;
  ASSERT_EQ(run(with_times(json, 0.0, 60.0, 30), dir.path() / "first", err), 0) << err;
  fs::create_directories(dir.path() / "damaged");
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin"})
    fs::copy_file(dir.path() / "first" / name, dir.path() / "damaged" / name);
  const fs::path checkpoint = dir.path() / "damaged" / "checkpoint.bin";
  patch_nan(checkpoint, fs::file_size(checkpoint) - 8 - sizeof(cplx));

  nlohmann::json restart = with_times(json, 60.0, 120.0, 30);
  restart["NewCalc"] = false;
  restart["RestartDir"] = (dir.path() / "first").string();
  ASSERT_EQ(run(restart, dir.path() / "clean", err), 0) << err;
  EXPECT_EQ(err, "");
  restart["RestartDir"] = (dir.path() / "damaged").string();
  EXPECT_EQ(run(restart, dir.path() / "nan", err), 5) << err;

  EXPECT_EQ(std::count(err.begin(), err.end(), '\n'), 1) << err;  // one line, nothing else on stderr
  // Every S is NaN now; the first in (is, v) order is the first analysed rotation of the propagation basis
  // (j = 4, mJ = -4), open at the first energy.
  EXPECT_NE(err.find("Non-finite result: the on-the-fly analysis, energy 1 of 41 (E_coll = 0.200000 eV): "
                     "S of (v, j, mJ, n, m) = (0, 4, -4, 0, 0) is (nan, nan)"),
            std::string::npos)
      << err;
  EXPECT_EQ(file_set(dir.path() / "nan"), file_set(dir.path() / "clean"));
  EXPECT_NE(read_file(dir.path() / "nan" / "s_matrix.tsv").find("nan"), std::string::npos);
}

// A restart from such a checkpoint with alsoAnalysis: the fine analysis of the BK file the run has just written
// follows the on-the-fly analysis, and each result set holding a NaN gets its line (the on-the-fly one first); the
// files of the clean run, exit code 5. 60 records (dtStep = dtAna = 10, 0 .. 600 a.u.), so that the fine grid has
// energies in the confidence interval (nen95 = 3; with 12 records there are none, and no fine results).
TEST(NonFiniteResults, AlsoAnalysisRunNamesEachResultSet) {
  const test_support::TemporaryDirectory dir("non_finite_also");
  nlohmann::json json = short_json();
  json["dtStep"] = 10.0;
  json["dtAna"] = 10.0;
  json["dtAve"] = 100.0;
  std::string err;
  ASSERT_EQ(run(with_times(json, 0.0, 300.0, 30), dir.path() / "first", err), 0) << err;
  fs::create_directories(dir.path() / "damaged");
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin"})
    fs::copy_file(dir.path() / "first" / name, dir.path() / "damaged" / name);
  const fs::path checkpoint = dir.path() / "damaged" / "checkpoint.bin";
  patch_nan(checkpoint, fs::file_size(checkpoint) - 8 - sizeof(cplx));

  nlohmann::json restart = with_times(json, 300.0, 600.0, 30);
  restart["NewCalc"] = false;
  restart["alsoAnalysis"] = true;
  restart["RestartDir"] = (dir.path() / "first").string();
  ASSERT_EQ(run(restart, dir.path() / "clean", err), 0) << err;
  EXPECT_EQ(err, "");
  restart["RestartDir"] = (dir.path() / "damaged").string();
  EXPECT_EQ(run(restart, dir.path() / "nan", err), 5) << err;

  std::istringstream lines(err);
  std::vector<std::string> printed;
  for (std::string line; std::getline(lines, line);) printed.push_back(line);
  ASSERT_EQ(printed.size(), 3u) << err;  // the on-the-fly results, the fine spectrum, its interpolated energies
  EXPECT_TRUE(printed[0].starts_with("Non-finite result: the on-the-fly analysis, ")) << err;
  EXPECT_TRUE(printed[1].starts_with("Non-finite result: the fine analysis (spectrum), ")) << err;
  EXPECT_TRUE(printed[2].starts_with("Non-finite result: the fine analysis (interpolated), ")) << err;
  EXPECT_EQ(file_set(dir.path() / "nan"), file_set(dir.path() / "clean"));
}

// An analysis-only run of a BK file whose last record has a NaN in the coefficient of the initial
// state: the time-to-energy transform spreads it over the spectrum and the interpolated energies.
// One line per result set, the files of the clean run, exit code 5. The run of FineAnalysisRuns in
// test_restart.cpp (dtStep = dtAna = 10, 60 records: all used, 60 is an FFT length).
TEST(NonFiniteResults, OnlyAnalysisOfANaNBkFileExitsWith5) {
  const test_support::TemporaryDirectory dir("non_finite_analysis");
  nlohmann::json json = short_json();
  json["dtStep"] = 10.0;
  json["dtAna"] = 10.0;
  json["dtAve"] = 100.0;
  std::string err;
  ASSERT_EQ(run(with_times(json, 0.0, 600.0, 60), dir.path() / "run", err), 0) << err;

  fs::create_directories(dir.path() / "damaged");
  for (const std::string name : {"checkpoint.bin", "bk_coefficients.bin"})
    fs::copy_file(dir.path() / "run" / name, dir.path() / "damaged" / name);
  // A record: the time, then the coefficients (is, v); the initial state (v = VIn = 0) is k = is0.
  const nlohmann::json basis = nlohmann::json::parse(read_file(dir.path() / "run" / "analysis.json")).at("analysis");
  const std::size_t n_coefficients =
      basis.at("n_states").get<std::size_t>() * basis.at("n_vibrational_levels").get<std::size_t>();
  const std::size_t initial_state = basis.at("initial_state_index").get<std::size_t>();
  const std::size_t record_size = sizeof(double) + n_coefficients * sizeof(cplx);
  const fs::path bk = dir.path() / "damaged" / "bk_coefficients.bin";
  patch_nan(bk, fs::file_size(bk) - record_size + sizeof(double) + initial_state * sizeof(cplx));

  nlohmann::json only = with_times(json, 0.0, 600.0, 60);
  only["onlyAnalysis"] = true;
  only["RestartDir"] = (dir.path() / "run").string();
  ASSERT_EQ(run(only, dir.path() / "clean", err), 0) << err;
  EXPECT_EQ(err, "");
  only["RestartDir"] = (dir.path() / "damaged").string();
  EXPECT_EQ(run(only, dir.path() / "nan", err), 5) << err;

  EXPECT_EQ(std::count(err.begin(), err.end(), '\n'), 2) << err;
  EXPECT_EQ(count_lines_with(err, "Non-finite result"), 2) << err;
  EXPECT_NE(err.find("Non-finite result: the fine analysis (spectrum), energy 1 of "), std::string::npos) << err;
  EXPECT_NE(err.find("Non-finite result: the fine analysis (interpolated), energy 1 of "), std::string::npos) << err;
  EXPECT_EQ(count_lines_with(err, "S of (v, j, mJ, n, m) = (0, 0, 0, 0, 0) is (nan, nan)"), 2) << err;
  EXPECT_EQ(file_set(dir.path() / "nan"), file_set(dir.path() / "clean"));
}
