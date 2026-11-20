// Tests: the norm budget (NormBudget; analysis/norm_budget.*, report/norm_budget_files.*): one window by hand with a
// negative test, the absorber list, and on runs of a small hard-wall input the averages unchanged by the budget, the
// same bytes with 1 and 8 threads, the projected size and the memory item. The balance against the norm and the
// numpy recomputation run on the plot fixtures (tools/plots/tests/test_runs_norm_budget.py).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "analysis/norm_budget.hpp"
#include "counting_backend.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation_tables.hpp"
#include "report/norm_budget_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

// The configuration, grids and tables of an input (the asymptotic potential plays no part in the absorbers).
struct BudgetSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  RotationalBasis basis;
  PropagationTables tables;

  explicit BudgetSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        basis(angular, config.input.j_in, config.input.mj_in),
        tables(make_propagation_tables(config, grids, basis, std::vector<double>(config.input.nr, 0.01))) {}
};

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

// The first index whose factor is below 1 (inside an absorber), the strongest (the smallest factor), or the first
// whose factor is exactly 1 (outside all).
int first_absorbing(const std::vector<double>& factors) {
  for (int i = 0; i < static_cast<int>(factors.size()); ++i)
    if (factors[i] < 1.0) return i;
  return -1;
}
int strongest(const std::vector<double>& factors) {
  return static_cast<int>(std::min_element(factors.begin(), factors.end()) - factors.begin());
}
int first_free(const std::vector<double>& factors) {
  for (int i = 0; i < static_cast<int>(factors.size()); ++i)
    if (factors[i] == 1.0) return i;
  return -1;
}

std::string file_bytes(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream content;
  content << file.rdbuf();
  return content.str();
}

// Runs the program on `json` with `threads` OpenMP threads into `dir`; its printed output goes to dir/stdout.txt.
void run(const nlohmann::json& json, const fs::path& dir, int threads) {
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
  std::ofstream(dir / "stdout.txt") << out.str();
}

}  // namespace

// The absorbers of the hard-wall input: ZOpt, ROpt on the scattering grid, ZSpOpt, ROpt on the specular grid; with
// both low-Z absorbers on, each grid's low one after its high one.
TEST(NormBudget, AbsorberList) {
  const BudgetSetup plain(small_hardwall_input());
  std::vector<std::string> names;
  const NormBudget plain_budget(plain.config, plain.grids, plain.tables);
  for (const BudgetAbsorber& a : plain_budget.absorbers()) names.push_back(a.name + "/" + a.grid);
  EXPECT_EQ(names,
            (std::vector<std::string>{"ZOpt/scattering", "ROpt/scattering", "ZSpOpt/specular", "ROpt/specular"}));

  const BudgetSetup low(with_low_z(small_hardwall_input()));
  names.clear();
  const NormBudget low_budget(low.config, low.grids, low.tables);
  for (const BudgetAbsorber& a : low_budget.absorbers()) names.push_back(a.name + "/" + a.grid);
  EXPECT_EQ(names, (std::vector<std::string>{"ZOpt/scattering", "ZOptLow/scattering", "ROpt/scattering",
                                             "ZSpOpt/specular", "ZSpOptLow/specular", "ROpt/specular"}));
}

// One window by hand on the hard-wall input (n time steps per analysis step; a cell's density factor over the
// window F = f^2n, gamma = -ln F). Three non-zero cells:
//   A in the Z absorber alone, without motion (rho_end = F rho_start): it loses exactly rho_start (1 - F), and both
//     one-sided estimates agree (error 0);
//   B in the Z and the r absorber (f_z f_r), probability flowing out: loss by the constant-inflow formula, shared
//     in the ratio ln f_z : ln f_r;
//   C on the specular grid at its Z absorber's strongest point, steady (rho_end = rho_start, inflow balancing the
//     loss): it loses exactly rho gamma, while the mean of the one-sided estimates is off (the negative test: by far
//     more than the tolerance), as is the amplitude's factor f^n in place of the density's (a realistic slip).
// The first call starts the window and returns nothing; the norms are the cells' sums.
TEST(NormBudget, OneWindowByHand) {
  const BudgetSetup s(small_hardwall_input());
  const int nr = s.config.input.nr;
  const int n = s.config.derived.steps_per_analysis;
  const int iz = first_absorbing(s.tables.z_optical);
  const int ir_free = first_free(s.tables.r_optical);
  const int ir_r = first_absorbing(s.tables.r_optical);
  const int izsp = strongest(s.tables.z_specular_optical);
  ASSERT_GE(iz, 0);
  ASSERT_GE(ir_free, 0);
  ASSERT_GE(ir_r, 0);
  ASSERT_GE(izsp, s.grids.z_scattering.size());  // the specular grid's Z absorber lies beyond the scattering grid

  const double fz = s.tables.z_optical[iz];
  const double fr = s.tables.r_optical[ir_r];
  const double fsp = s.tables.z_specular_optical[izsp];
  const double big_fa = std::pow(fz, 2 * n);
  const double big_fb = std::pow(fz * fr, 2 * n);
  const double big_fc = std::pow(fsp, 2 * n);
  const double gamma_b = -std::log(big_fb);
  const double gamma_c = -std::log(big_fc);

  CellProbabilities start{std::vector<double>(std::size_t(s.grids.z_scattering.size()) * nr, 0.0),
                          std::vector<double>(std::size_t(s.grids.z_specular.size()) * nr, 0.0)};
  CellProbabilities end = start;
  const std::size_t a = ir_free + std::size_t(nr) * iz;
  const std::size_t b = ir_r + std::size_t(nr) * iz;
  const std::size_t c = ir_free + std::size_t(nr) * izsp;
  start.scattering[a] = 0.3;
  end.scattering[a] = 0.3 * big_fa;
  start.scattering[b] = 0.2;
  end.scattering[b] = 0.1;
  start.specular[c] = 0.4;
  end.specular[c] = 0.4;

  NormBudget budget(s.config, s.grids, s.tables);
  EXPECT_FALSE(budget.add(start).has_value());
  const std::optional<NormBudgetWindow> w = budget.add(end);
  ASSERT_TRUE(w.has_value());

  const double loss_a = 0.3 * (1.0 - big_fa);
  const double loss_b = 0.2 * (1.0 - big_fb) + (0.1 - big_fb * 0.2) * (gamma_b / (1.0 - big_fb) - 1.0);
  const double loss_c = 0.4 * gamma_c;
  const double z_share = std::log(fz) / (std::log(fz) + std::log(fr));
  const double tolerance = 1e-14;
  // Absorbers: 0 ZOpt, 1 ROpt (scattering), 2 ZSpOpt, 3 ROpt (specular).
  EXPECT_NEAR(w->loss[0], loss_a + z_share * loss_b, tolerance);
  EXPECT_NEAR(w->loss[1], (1.0 - z_share) * loss_b, tolerance);
  EXPECT_NEAR(w->loss[2], loss_c, tolerance);
  EXPECT_EQ(w->loss[3], 0.0);
  // The indicator: gamma^2 / 6 of the inflow rho_end - rho_start + loss (A, without motion, has none), the leading
  // distance between the loss and the mean of the one-sided estimates (within 10 % of it at C's gamma).
  EXPECT_NEAR(w->loss_error[2], gamma_c * gamma_c / 6.0 * loss_c, tolerance);
  EXPECT_NEAR(w->loss_error[0], z_share * gamma_b * gamma_b / 6.0 * std::abs(0.1 - 0.2 + loss_b), tolerance);
  const double mean_c = 0.5 * (0.4 * (1.0 - big_fc) + 0.4 * (1.0 / big_fc - 1.0));
  EXPECT_NEAR(w->loss_error[2], std::abs(loss_c - mean_c), 0.1 * std::abs(loss_c - mean_c));
  EXPECT_DOUBLE_EQ(w->norm_start, 0.9);
  EXPECT_DOUBLE_EQ(w->norm_end, 0.3 * big_fa + 0.1 + 0.4);

  EXPECT_GT(std::abs(mean_c - loss_c), 1e3 * tolerance);  // the mean of the one-sided estimates is biased
  const double amplitude_slip = 0.4 * (1.0 - std::pow(fsp, n));  // f^n: the amplitude's factor, not the density's
  EXPECT_GT(std::abs(w->loss[2] - amplitude_slip), 1e6 * tolerance);
}

// The low-Z absorbers (both on) take the Z part of their own cells: one steady cell inside ZOptLow (outside the r
// absorber) loses rho gamma = -2n rho ln f to ZOptLow (index 1), nothing to ZOpt; giving it to the high-end absorber
// (a slip) fails.
TEST(NormBudget, LowZAbsorberTakesItsOwnCells) {
  const BudgetSetup s(with_low_z(small_hardwall_input()));
  const int nr = s.config.input.nr;
  const int n = s.config.derived.steps_per_analysis;
  ASSERT_TRUE(s.grids.z_optical_low.has_value());
  const int iz = s.grids.z_optical_low->first;  // the low absorber's strongest point
  const int ir = first_free(s.tables.r_optical);
  ASSERT_LT(s.tables.z_optical[iz], 1.0);
  CellProbabilities cells{std::vector<double>(std::size_t(s.grids.z_scattering.size()) * nr, 0.0),
                          std::vector<double>(std::size_t(s.grids.z_specular.size()) * nr, 0.0)};
  cells.scattering[ir + std::size_t(nr) * iz] = 0.5;
  NormBudget budget(s.config, s.grids, s.tables);
  EXPECT_FALSE(budget.add(cells).has_value());
  const std::optional<NormBudgetWindow> w = budget.add(cells);
  ASSERT_TRUE(w.has_value());
  const double expected = -0.5 * 2 * n * std::log(s.tables.z_optical[iz]);
  ASSERT_GT(expected, 1e-6);
  EXPECT_NEAR(w->loss[1], expected, 1e-14);  // 1: ZOptLow (AbsorberList)
  EXPECT_EQ(w->loss[0], 0.0);                // 0: ZOpt
}

// The budget computes the averages' densities once and hands them to both: averages.dat is bytewise the run's
// without NormBudget. norm_budget/ does not depend on the thread count (1 and 8 threads write the same bytes).
TEST(NormBudgetFiles, AveragesUnchangedAndOneAndEightThreadsWriteTheSameBytes) {
  TemporaryDirectory dir("norm_budget_runs");
  nlohmann::json json = small_hardwall_input();
  run(json, dir.path() / "plain", 4);
  json["NormBudget"] = true;
  run(json, dir.path() / "one", 1);
  run(json, dir.path() / "eight", 8);
  EXPECT_EQ(file_bytes(dir.path() / "plain" / "averages.dat"), file_bytes(dir.path() / "one" / "averages.dat"));
  for (const char* name :
       {"manifest.json", "times.bin", "loss.bin", "loss_error.bin", "norm_start.bin", "norm_end.bin"}) {
    const std::string one = file_bytes(dir.path() / "one" / "norm_budget" / name);
    EXPECT_FALSE(one.empty()) << name;
    EXPECT_EQ(one, file_bytes(dir.path() / "eight" / "norm_budget" / name)) << name;
  }
}

// The projected size is the written size: 100 windows (one per analysis step, none at tStart) of the four datasets
// (4 absorbers) and the time. The memory item is what the budget keeps (NormBudget::memory_bytes).
TEST(NormBudgetFiles, ProjectedSizeAndMemoryItem) {
  TemporaryDirectory dir("norm_budget_size");
  nlohmann::json json = small_hardwall_input();
  json["NormBudget"] = true;
  run(json, dir.path(), 4);
  const BudgetSetup s(json);
  const std::optional<PlotDataSize> size = norm_budget_data_size(s.config);
  ASSERT_TRUE(size.has_value());
  std::uintmax_t written = 0;
  for (const char* name : {"times.bin", "loss.bin", "loss_error.bin", "norm_start.bin", "norm_end.bin"})
    written += fs::file_size(dir.path() / "norm_budget" / name);
  EXPECT_EQ(size->bytes, written);
  EXPECT_EQ(written, std::uintmax_t(small_hardwall_analysis_steps) * (4 + 4 + 1 + 1 + 1) * sizeof(double));
  // The run prints the measured cost: 100 windows, the seconds and their share of the propagation.
  EXPECT_NE(file_bytes(dir.path() / "stdout.txt").find("Norm budget (NormBudget): 100 windows written, "),
            std::string::npos);

  const MemoryEstimate estimate = estimate_memory(s.config, s.grids, s.angular, 1);
  const std::size_t nr = s.config.input.nr;
  const std::size_t scattering = s.grids.z_scattering.size() * nr;
  const std::size_t cells = scattering + s.grids.z_specular.size() * nr;
  int found = 0;
  for (const MemoryItem& item : estimate.items) {
    if (item.name != "Norm budget (NormBudget)") continue;
    ++found;
    EXPECT_EQ(item.bytes, cells * (4 * sizeof(double) + 2 * sizeof(int)) + 2 * scattering * sizeof(double));
  }
  EXPECT_EQ(found, 1);

  // With both low-Z absorbers the records hold 6 losses each (two per absorber dataset more).
  const BudgetSetup low(with_low_z(json));
  EXPECT_EQ(norm_budget_data_size(low.config)->bytes,
            std::uintmax_t(small_hardwall_analysis_steps) * (6 + 6 + 1 + 1 + 1) * sizeof(double));
  EXPECT_EQ(NormBudget(low.config, low.grids, low.tables).absorbers().size(), 6u);
}
