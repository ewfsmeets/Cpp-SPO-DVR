// write2DPES: the cut of the PES (compute_pes_cut), its file (pes_cut_2d.dat) and the run that writes it
// and stops. Its values against the Fortran: gate 1c in test_regression_pes.cpp.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "core/unit_conversion.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/pes_cut.hpp"
#include "potential/potential_factory.hpp"
#include "report/pes_cut_file.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using Json = nlohmann::json;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

// The Dai-Light reference case (its PES depends on x, y, theta and phi) with a write2DPES cut.
Json cut_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_small/input.json");
  Json json = Json::parse(file);
  json["write2DPES"] = true;
  json["x2DPES"] = 1.452;
  json["y2DPES"] = 0.484;
  json["theta2DPES"] = 37.0;
  json["phi2DPES"] = 113.0;
  return json;
}

// The output of run_program on `json` (written to dir/input.json) into dir/out: exit status, stdout, stderr.
struct RunResult {
  int status;
  std::string out;
  std::string err;
};
RunResult run(const TemporaryDirectory& dir, const Json& json) {
  std::ofstream(dir.path() / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir.path() / "input.json";
  cl.output_dir = dir.path() / "out";
  std::ostringstream out;
  std::ostringstream err;
  const int status = run_program(cl, out, err, false);
  return {status, out.str(), err.str()};
}

// The names of the files in a directory.
std::set<std::string> files_in(const fs::path& dir) {
  std::set<std::string> names;
  for (const auto& entry : fs::directory_iterator(dir)) names.insert(entry.path().filename().string());
  return names;
}

}  // namespace

// The cut holds V of the PES at every (r, Z) of the r grid and the whole scattering Z grid (also Z >= zasymp,
// with no switching: the case switches from 6 to 8 bohr), at the angles in radians the Fortran passes
// (degrees x Constants_Deg2Rad), Z fastest.
TEST(PesCut, ValuesAreThePesItselfInTheFortranOrder) {
  const Config config = build_config(parse_input(cut_json()));
  ASSERT_TRUE(config.input.use_switching);
  const SpatialGrids grids = make_spatial_grids(config);
  const auto potential = make_potential(config.input);
  const PesCut cut = compute_pes_cut(*potential, config.input, grids);
  ASSERT_EQ(cut.r.size(), static_cast<std::size_t>(grids.r.size()));
  ASSERT_EQ(cut.z.size(), static_cast<std::size_t>(grids.z_scattering.size()));
  ASSERT_EQ(cut.v_ev.size(), cut.r.size() * cut.z.size());
  EXPECT_EQ(cut.theta, 37.0 * Constants::deg_to_rad);
  EXPECT_EQ(cut.phi, 113.0 * Constants::deg_to_rad);
  int differing = 0;
  for (int ir = 0; ir < grids.r.size(); ++ir)
    for (int iz = 0; iz < grids.z_scattering.size(); ++iz) {
      const double v = potential->value(1.452, 0.484, grids.z_scattering[iz], grids.r[ir], cut.theta, cut.phi, iz, ir);
      differing += (cut.v(ir, iz) == v && cut.v_ev[iz + grids.z_scattering.size() * ir] == v) ? 0 : 1;
    }
  EXPECT_EQ(differing, 0);
  EXPECT_GT(cut.z.back(), config.derived.z_asymptotic);  // the cut also covers the gas phase
}

// pes_cut_2d.dat: four '#' lines, then per r a block of nZ lines "r Z V(hartree) V(eV)" (%.16e: they read back
// exactly) and a blank line; V(hartree) = V(eV) / 27.2114.
TEST(PesCut, FileHoldsEveryPointInBlocksOfR) {
  TemporaryDirectory dir("pes_cut_file");
  PesCut cut{0.5, 0.25, 90.0, 45.0, Constants::pi / 2.0, Constants::pi / 4.0, {1.0, 1.5}, {-1.0, 0.0, 2.0}, {}};
  for (int i = 0; i < 6; ++i) cut.v_ev.push_back(0.1 * i - 0.25);
  write_pes_cut_file(dir.path() / "cut.dat", cut);
  std::ifstream file(dir.path() / "cut.dat");
  std::string line;
  int header = 0;
  int blank = 0;
  std::vector<std::vector<double>> rows;
  while (std::getline(file, line)) {
    if (line.starts_with("#")) {
      ++header;
      continue;
    }
    if (line.empty()) {
      ++blank;
      EXPECT_EQ(rows.size() % 3, 0u) << "a blank line inside an r block";
      continue;
    }
    std::istringstream fields(line);
    std::vector<double> row(4);
    fields >> row[0] >> row[1] >> row[2] >> row[3];
    rows.push_back(row);
  }
  EXPECT_EQ(header, 4);
  EXPECT_EQ(blank, 2);
  ASSERT_EQ(rows.size(), 6u);
  for (int ir = 0; ir < 2; ++ir)
    for (int iz = 0; iz < 3; ++iz) {
      const std::vector<double>& row = rows[iz + 3 * ir];
      EXPECT_EQ(row[0], cut.r[ir]);
      EXPECT_EQ(row[1], cut.z[iz]);
      EXPECT_EQ(row[2], cut.v(ir, iz) / 27.2114);
      EXPECT_EQ(row[3], cut.v(ir, iz));
    }
}

// A write2DPES run writes the cut, V_asy and run_info.json (with the cut, without the gas phase) and stops
// with exit status 0 before the 6D potential and the gas phase, as the Fortran does.
TEST(PesCut, RunWritesTheCutAndStops) {
  TemporaryDirectory dir("pes_cut_run");
  const RunResult result = run(dir, cut_json());
  ASSERT_EQ(result.status, 0) << result.err;
  EXPECT_NE(result.out.find("2D CUT OF THE PES"), std::string::npos);
  EXPECT_EQ(result.out.find("6D potential computed"), std::string::npos);
  const fs::path out = dir.path() / "out";
  EXPECT_EQ(files_in(out), (std::set<std::string>{"asymptotic_potential.dat", "pes_cut_2d.dat", "run_info.json"}));
  std::ifstream info_file(out / "run_info.json");
  const Json info = Json::parse(info_file);
  EXPECT_EQ(info["pes_cut"]["file"], "pes_cut_2d.dat");
  EXPECT_EQ(info["pes_cut"]["theta_deg"], 37.0);
  EXPECT_TRUE(info["analysis_energies"]["collision_ev"].is_null());  // no gas phase
}

// The Fortran writes the cut in every run mode before it reads a restart file: so does a restart whose
// RestartDir does not exist. Plot data keys are refused: the run stops before anything writes them.
TEST(PesCut, RunNeedsNoRestartFilesAndRefusesPlotData) {
  TemporaryDirectory dir("pes_cut_restart");
  Json json = cut_json();
  json["NewCalc"] = false;
  // A restart time and its end (the input checks of a restart still apply: tStart > 0, tEnd = tStart + nStep dt).
  json["tStart"] = 1000.0;
  json["tEnd"] = 1000.0 + json["nStep"].get<double>() * json["dtStep"].get<double>();
  json["RestartDir"] = (dir.path() / "no_such_directory").string();
  EXPECT_EQ(run(dir, json).status, 0);
  Json plot = cut_json();
  plot["writePotential"] = true;
  const RunResult refused = run(dir, plot);
  EXPECT_EQ(refused.status, 2);
  EXPECT_NE(refused.err.find("writePotential"), std::string::npos) << refused.err;
}
