// The runs that the plot tool's tests read (tools/plots, ctest PlotRuns.*): real spodvr runs of small
// inputs, written by the current program into SPODVR_PLOT_FIXTURES_DIR. One ctest in one process
// (PlotRuns.Fixtures, FIXTURES_SETUP plot_fixtures, out of the quick tier); each case deletes its
// directories first. The fixtures are deliberately not converged (as the gate case hardwall_small is not).
// Every fixture writes the plot data: the propagating ones densities/ (writeDensities), frames/ (dtFrames
// 100 a.u., with every FrameExtras extra) and norm_budget/ (NormBudget), and all potential/ (writePotential, with
// PlotSites: hardwall_tiny "origin" and "centre", dailight_oblique "top" and "bridge").
//   hardwall_tiny             hardwall_small with nX = nY = 2, JMax = MJMax = 2, tEnd 10000, alsoAnalysis
//   hardwall_tiny_split       the same in two segments, split at t = 4000 (first/, second/; split_time.txt)
//   hardwall_budget           hardwall_tiny with dtAve = dtAna = 10 (a densities/ record at every window's end, for
//                             recomputing the norm budget) ...
//   hardwall_budget_dt4       ... and the same with dtAna = dtAve = 4 (the budget's error falls with dtAna^2)
//   dailight_oblique          cu111_dailight_small with nX = nY = 4, JMax = MJMax = 2, EParallel 0.1 eV at
//                             17 degrees, nZ 192, the Z absorber at 10.16-13.28 (its length, ending on the last
//                             point) and tEnd 14000 (JMax 2: a run takes about 65 s; with JMax 4 it took 183 s)
//   dailight_oblique_split    the same in two segments, split at the multiple of 100 a.u. nearest the time
//                             of minimum <Z> of dailight_oblique (split_time.txt)
//   potential_only            hardwall_tiny with Potential: true (potential/ only)
//   hardwall_long             hardwall_tiny to tEnd 20000 with only the energy-resolved sums at "origin" and the
//                             collision energies 0.3, 0.4, 0.5 eV (energy_resolved/; the hard-wall gate of rho_E)
//   dailight_rho_e            dailight_oblique with only the energy-resolved sums, at all 16 grid points (PlotSites
//                             g<ix><iy>) and the analysis energies 0.45 and 0.55 eV (the S-matrix gate of rho_E)
//   pes_cut                   dailight_oblique with write2DPES at the bridge (x = alat / 2, y = 0), at the last
//                             theta node of the run (in degrees) and phi 0, so that the plot tool's tests can
//                             compare it with dailight_oblique's potential/ at the bridge site: pes_cut_2d.dat,
//                             asymptotic_potential.dat and run_info.json only (no plot keys)
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "core/unit_conversion.hpp"
#include "grids/angular_grid.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

const fs::path fixtures_dir = SPODVR_PLOT_FIXTURES_DIR;

nlohmann::json reference_input(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return nlohmann::json::parse(file);
}

// Runs the program on `json` into `dir` (recreated; the input is kept there as input.json); threads:
// OMP_NUM_THREADS, as every test. Returns the exit code; the program's error output goes to the test log.
int run(const nlohmann::json& json, const fs::path& dir) {
  fs::remove_all(dir);
  fs::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir;
  std::ostringstream out;
  std::ostringstream err;
  const int code = run_program(cl, out, err, false);
  std::ofstream(dir / "stdout.txt") << out.str();
  if (code != 0) std::cerr << err.str();
  return code;
}

// The times of a run: tStart, tEnd and nStep for steps of dtStep.
nlohmann::json with_times(nlohmann::json json, double t_start, double t_end) {
  const double dt = json["dtStep"];
  json["tStart"] = t_start;
  json["tEnd"] = t_end;
  json["nStep"] = static_cast<int>(std::lround((t_end - t_start) / dt));
  return json;
}

// The second segment of a split run: a restart from `first`.
nlohmann::json restart_from(nlohmann::json json, const fs::path& first) {
  json["NewCalc"] = false;
  json["RestartDir"] = first.string();
  return json;
}

// Every extra of FrameExtras.
nlohmann::json every_frame_extra() {
  return nlohmann::json::array({"angular", "diffraction", "vibration", "rotation_moments", "momentum"});
}

nlohmann::json hardwall_tiny() {
  nlohmann::json json = reference_input("hardwall_small");
  json["nX"] = 2;
  json["nY"] = 2;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  json["ImposeVMax"] = false;
  json["alsoAnalysis"] = true;
  json["writeProbs"] = true;
  json["writeRDProbs"] = true;
  json["writeDensities"] = true;
  json["writePotential"] = true;
  json["dtFrames"] = 100.0;
  json["FrameExtras"] = every_frame_extra();
  json["NormBudget"] = true;
  json["PlotSites"] = {{{"name", "origin"}, {"x", 0.0}, {"y", 0.0}}, {{"name", "centre"}, {"x", 0.5}, {"y", 0.5}}};
  return with_times(json, 0.0, 10000.0);
}

const double dailight_t_end = 14000.0;

nlohmann::json dailight_oblique() {
  nlohmann::json json = reference_input("cu111_dailight_small");
  json["nX"] = 4;
  json["nY"] = 4;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  json["EParallel"] = 0.1;
  json["PhiParallel"] = 17.0;
  json["nZ"] = 192;  // last point -2 + 191 * 0.08 = 13.28
  json["ZOptStart"] = 10.16;
  json["ZOptEnd"] = 13.28;
  json["alsoAnalysis"] = true;
  json["writeDensities"] = true;
  json["dtFrames"] = 100.0;
  json["FrameExtras"] = every_frame_extra();
  json["NormBudget"] = true;
  json["writePotential"] = true;
  json["PlotSites"] = {{{"name", "top"}, {"x", 0.0}, {"y", 0.0}}, {{"name", "bridge"}, {"x", 0.5}, {"y", 0.0}}};
  return with_times(json, 0.0, dailight_t_end);
}

// The time of the smallest <Z> in averages.dat (columns t, norm, <Z>, <r>, rest norm).
double time_of_minimum_z(const fs::path& averages) {
  std::ifstream in(averages);
  std::string line;
  double best_time = NAN;
  double best_z = INFINITY;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    double t = NAN;
    double norm = NAN;
    double z = NAN;
    row >> t >> norm >> z;
    if (z < best_z) {
      best_z = z;
      best_time = t;
    }
  }
  return best_time;
}

// The files every run with alsoAnalysis writes. (A run that stopped early has only run_info.json, the
// gas-phase files and part of bk_coefficients.bin: the plot tool's tests make such a copy themselves.)
void expect_run_files(const fs::path& dir) {
  for (const std::string name : {"run_info.json", "averages.dat", "bk_coefficients.bin", "checkpoint.bin",
                                 "scattering_summary.tsv", "s_matrix.tsv", "gas_phase_states.dat",
                                 "asymptotic_potential.dat", "densities/times.bin", "potential/v_sites.bin"})
    EXPECT_TRUE(fs::exists(dir / name)) << dir / name;
  EXPECT_TRUE(fs::exists(dir / "fine_analysis")) << dir;
}

}  // namespace

TEST(PlotRuns, HardwallTiny) {
  const fs::path dir = fixtures_dir / "hardwall_tiny";
  ASSERT_EQ(run(hardwall_tiny(), dir), 0);
  expect_run_files(dir);
}

TEST(PlotRuns, HardwallTinySplit) {
  const fs::path root = fixtures_dir / "hardwall_tiny_split";
  const double split = 4000.0;  // the packet near the wall; on the dtAve grid
  fs::remove_all(root);
  fs::create_directories(root);
  std::ofstream(root / "split_time.txt") << split << '\n';
  ASSERT_EQ(run(with_times(hardwall_tiny(), 0.0, split), root / "first"), 0);
  ASSERT_EQ(run(restart_from(with_times(hardwall_tiny(), split, 10000.0), root / "first"), root / "second"), 0);
  expect_run_files(root / "first");
  expect_run_files(root / "second");
}

// hardwall_tiny with the averages (and so densities/) at every analysis step of dtAna: the budget's windows and the
// cells they start and end at.
nlohmann::json hardwall_budget(double dt_analysis) {
  nlohmann::json json = hardwall_tiny();
  json["dtAna"] = dt_analysis;
  json["dtAve"] = dt_analysis;
  return json;
}

TEST(PlotRuns, HardwallBudget) {
  for (const auto& [name, dt_analysis] : {std::pair{"hardwall_budget", 10.0}, std::pair{"hardwall_budget_dt4", 4.0}}) {
    const fs::path dir = fixtures_dir / name;
    ASSERT_EQ(run(hardwall_budget(dt_analysis), dir), 0) << name;
    EXPECT_TRUE(fs::exists(dir / "norm_budget" / "times.bin")) << name;
    EXPECT_TRUE(fs::exists(dir / "densities" / "times.bin")) << name;
  }
}

TEST(PlotRuns, DailightOblique) {
  const fs::path dir = fixtures_dir / "dailight_oblique";
  ASSERT_EQ(run(dailight_oblique(), dir), 0);
  expect_run_files(dir);
}

// After DailightOblique (declaration order; one process): the split time comes from its averages.
TEST(PlotRuns, DailightObliqueSplit) {
  const fs::path unsplit = fixtures_dir / "dailight_oblique" / "averages.dat";
  ASSERT_TRUE(fs::exists(unsplit)) << "run PlotRuns.DailightOblique first";
  const double split = 100.0 * std::round(time_of_minimum_z(unsplit) / 100.0);
  ASSERT_GT(split, 0.0);
  ASSERT_LT(split, dailight_t_end);
  const fs::path root = fixtures_dir / "dailight_oblique_split";
  fs::remove_all(root);
  fs::create_directories(root);
  std::ofstream(root / "split_time.txt") << split << '\n';
  ASSERT_EQ(run(with_times(dailight_oblique(), 0.0, split), root / "first"), 0);
  const nlohmann::json second = restart_from(with_times(dailight_oblique(), split, dailight_t_end), root / "first");
  ASSERT_EQ(run(second, root / "second"), 0);
  expect_run_files(root / "first");
  expect_run_files(root / "second");
}

// The plot keys of the fixtures above, left out where a fixture writes only the energy-resolved sums.
void without_plot_keys(nlohmann::json& json) {
  for (const std::string key :
       {"alsoAnalysis", "writeDensities", "writePotential", "dtFrames", "FrameExtras", "NormBudget", "PlotSites"})
    json.erase(key);
}

TEST(PlotRuns, HardwallLong) {
  nlohmann::json json = with_times(hardwall_tiny(), 0.0, 20000.0);
  without_plot_keys(json);
  json["PlotSites"] = {{{"name", "origin"}, {"x", 0.0}, {"y", 0.0}}};
  json["EnergyResolvedSites"] = {"origin"};
  json["EnergyResolvedEnergies"] = {0.3, 0.4, 0.5};
  const fs::path dir = fixtures_dir / "hardwall_long";
  ASSERT_EQ(run(json, dir), 0);
  EXPECT_TRUE(fs::exists(dir / "energy_resolved" / "steps_summed.bin"));
}

TEST(PlotRuns, DailightRhoE) {
  nlohmann::json json = dailight_oblique();
  without_plot_keys(json);
  json["PlotSites"] = nlohmann::json::array();
  json["EnergyResolvedSites"] = nlohmann::json::array();
  for (int ix = 0; ix < 4; ++ix)
    for (int iy = 0; iy < 4; ++iy) {
      const std::string name = "g" + std::to_string(ix) + std::to_string(iy);
      json["PlotSites"].push_back({{"name", name}, {"x", ix / 4.0}, {"y", iy / 4.0}});
      json["EnergyResolvedSites"].push_back(name);
    }
  json["EnergyResolvedEnergies"] = {0.45, 0.55};  // analysis energies: 0.3 + 0.01 k collision
  const fs::path dir = fixtures_dir / "dailight_rho_e";
  ASSERT_EQ(run(json, dir), 0);
  EXPECT_TRUE(fs::exists(dir / "energy_resolved" / "steps_summed.bin"));
}

TEST(PlotRuns, PotentialOnly) {
  nlohmann::json json = hardwall_tiny();
  json["Potential"] = true;
  json.erase("alsoAnalysis");
  json.erase("writeDensities");  // densities/, frames/ and norm_budget/ belong to a propagating run; potential/ is
  json.erase("dtFrames");        // written at set-up
  json.erase("FrameExtras");
  json.erase("NormBudget");
  const fs::path dir = fixtures_dir / "potential_only";
  ASSERT_EQ(run(json, dir), 0);
  EXPECT_TRUE(fs::exists(dir / "run_info.json"));
  EXPECT_TRUE(fs::exists(dir / "potential" / "v_min.bin"));
  EXPECT_FALSE(fs::exists(dir / "averages.dat"));
}

TEST(PlotRuns, PesCut) {
  nlohmann::json json = dailight_oblique();
  for (const std::string key :
       {"alsoAnalysis", "writeDensities", "writePotential", "dtFrames", "FrameExtras", "NormBudget", "PlotSites"})
    json.erase(key);
  json["write2DPES"] = true;
  json["x2DPES"] = 2.42;  // the bridge: alat / 2 (alat 4.84 bohr)
  json["y2DPES"] = 0.0;
  const AngularGrid angular = make_angular_grid(build_config(parse_input(dailight_oblique())));
  json["theta2DPES"] = rad_to_deg(angular.thetas().back());  // a node of the run's angular grid
  json["phi2DPES"] = 0.0;  // the first phi node
  const fs::path dir = fixtures_dir / "pes_cut";
  ASSERT_EQ(run(json, dir), 0);
  EXPECT_TRUE(fs::exists(dir / "pes_cut_2d.dat"));
  EXPECT_TRUE(fs::exists(dir / "asymptotic_potential.dat"));
  EXPECT_FALSE(fs::exists(dir / "gas_phase_states.dat"));
}
