// Tests: run_info.json (report/run_info_file.*), written by every run that reaches the output
// directory: its fields against values computed here from the input, in every run mode.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "report/run_info_file.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

nlohmann::json reference_input(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return nlohmann::json::parse(file);
}

// One program run: exit code, printed output and the run_info.json it wrote (null if none).
struct ProgramRun {
  int code = -1;
  std::string out;
  std::string err;
  nlohmann::json info;
};

// Runs the program on `json` with output into `dir` (threads: OMP_NUM_THREADS, as every test) and
// reads `dir`/`info_file` (an analysis-only run writes it below fine_analysis/).
ProgramRun run(const nlohmann::json& json, const fs::path& dir, const fs::path& info_file = "run_info.json") {
  fs::create_directories(dir);
  const fs::path input = dir / "input.json";
  std::ofstream(input) << json.dump(1);
  CommandLine cl;
  cl.input_file = input;
  cl.output_dir = dir;
  std::ostringstream out;
  std::ostringstream err;
  ProgramRun r;
  r.code = run_program(cl, out, err, false);
  r.out = out.str();
  r.err = err.str();
  if (fs::exists(dir / info_file)) r.info = nlohmann::json::parse(std::ifstream(dir / info_file));
  return r;
}

// E(v, j) in hartree from gas_phase_states.dat (columns v, j, E, ψ…; `#` lines are comments): an
// independent reading of the program's own file.
double gas_file_energy(const fs::path& file, int v, int j) {
  std::ifstream in(file);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    int row_v = -1;
    int row_j = -1;
    double energy = NAN;
    row >> row_v >> row_j >> energy;
    if (row_v == v && row_j == j) return energy;
  }
  return NAN;
}

// The documented factors (docs/method.md, units::legacy in core/unit_conversion.hpp: the set the program uses).
constexpr double hartree_to_ev = 27.2114;
constexpr double amu_to_me = 1822.89;

// Fields of run_info.json fixed by the input and a documented unit, with the value expected here
// (hartree, electron masses, bohr, a.u. of time; eV where the name says so).
using Expected = std::vector<std::pair<std::string, double>>;

Expected expected_fields(const InputParameters& p, double initial_energy) {
  const double total_mass = (p.mass1 + p.mass2) * amu_to_me;
  return {{"/potential/v_max", p.v_max / hartree_to_ev},
          {"/potential/v_max_ev", p.v_max},
          {"/potential/v_min", p.v_min / hartree_to_ev},
          {"/potential/v_min_ev", p.v_min},
          {"/molecule/total_mass", total_mass},
          {"/molecule/reduced_mass", p.mass1 * p.mass2 / (p.mass1 + p.mass2) * amu_to_me},
          {"/molecule/mass1_amu", p.mass1},
          {"/initial_state/energy", initial_energy},
          {"/initial_state/energy_ev", initial_energy * hartree_to_ev},
          {"/wave_packet/e_min_ev", p.e_min_wave_packet},
          {"/wave_packet/e_max_ev", p.e_max_wave_packet},
          {"/wave_packet/e_parallel_ev", p.e_parallel},
          {"/wave_packet/phi_parallel_deg", p.phi_parallel},
          {"/wave_packet/z0", p.z0},
          {"/wave_packet/k_parallel", std::sqrt(2.0 * total_mass * p.e_parallel / hartree_to_ev)},
          {"/times/t_start", p.t_start},
          {"/times/t_end", p.t_end},
          {"/times/dt_step", p.time_step},
          {"/times/dt_analysis", p.dt_analysis},
          {"/times/dt_averages", p.dt_averages},
          {"/grids/z_start", p.z_start},
          {"/grids/dz", p.dz},
          {"/grids/r_start", p.r_start},
          {"/grids/dr", p.dr},
          {"/unit_factors/hartree_to_ev", hartree_to_ev},
          {"/unit_factors/amu_to_me", amu_to_me},
          {"/unit_factors/bohr_to_angstrom", 0.529177},
          {"/unit_factors/hartree_to_wavenumber", 219475.0}};
}

// Fields that follow from the input and the lattice the PES reports (a, b, skew), computed here with
// independent formulas: the (x, y) grid and cell, the packet's mean energy ((sqrt(Emin) + sqrt(Emax))/2)^2
// and alpha = (k0 - kmin)^2 / 4, k_total, kz = sqrt(k_total^2 - k_parallel^2), the incidence angles and
// the analysis line on the Z grid.
Expected expected_derived(const InputParameters& p, const nlohmann::json& info) {
  const double pi = std::acos(-1.0);
  const double a = info["lattice"]["a"];
  const double b = info["lattice"]["b"];
  const double skew = info["lattice"]["skew_angle_rad"];
  const double total_mass = (p.mass1 + p.mass2) * amu_to_me;
  const auto momentum = [&](double energy_ev) { return std::sqrt(2.0 * total_mass * energy_ev / hartree_to_ev); };
  const double mean_normal = std::pow((std::sqrt(p.e_min_wave_packet) + std::sqrt(p.e_max_wave_packet)) / 2.0, 2);
  const double translational = p.e_parallel + mean_normal;
  const double k_min = momentum(p.e_min_wave_packet);
  const double k_mean = (k_min + momentum(p.e_max_wave_packet)) / 2.0;
  const double k_total = momentum(translational);
  const double k_parallel = momentum(p.e_parallel);
  const int iz_ana = info["analysis"]["iz_ana"];
  return {{"/grids/dx", a / p.nx},
          {"/grids/dy", b / p.ny},
          {"/grids/xy_volume_element", a / p.nx * b / p.ny * std::sin(skew)},
          {"/grids/unit_cell_area", a * b * std::sin(skew)},
          {"/grids/unit_cell_norm", 1.0 / std::sqrt(a * b * std::sin(skew))},
          {"/lattice/delta_kx", 2.0 * pi / (a * std::sin(skew))},
          {"/lattice/delta_ky", 2.0 * pi / (b * std::sin(skew))},
          {"/wave_packet/translational_energy_ev", translational},
          {"/wave_packet/alpha", std::pow(k_mean - k_min, 2) / 4.0},
          {"/wave_packet/k_total", k_total},
          {"/wave_packet/kz", std::sqrt(k_total * k_total - k_parallel * k_parallel)},
          {"/wave_packet/theta_incidence_rad", std::asin(std::sqrt(p.e_parallel / translational))},
          {"/wave_packet/phi_incidence_rad", p.phi_parallel * pi / 180.0},
          {"/analysis/z_ana", p.z_start + iz_ana * p.dz}};
}

// The fields of `info` that differ from `expected` by more than 1e-12 relative (or 1e-14
// absolute near zero). A missing field, a non-number or NaN counts as different.
std::vector<std::string> mismatches(const nlohmann::json& info, const Expected& expected) {
  std::vector<std::string> wrong;
  for (const auto& [pointer, value] : expected) {
    const nlohmann::json::json_pointer at(pointer);
    if (!info.contains(at) || !info.at(at).is_number()) {
      wrong.push_back(pointer);
      continue;
    }
    const double got = info.at(at).get<double>();
    const bool close = std::abs(got - value) <= std::max(1e-12 * std::abs(value), 1e-14);  // NaN: false
    if (!close) wrong.push_back(pointer);
  }
  return wrong;
}

// The fields every run_info.json must agree with: its input, its gas-phase file, the program.
void expect_matches_input(const ProgramRun& r, const nlohmann::json& input, const fs::path& dir) {
  ASSERT_EQ(r.code, 0) << r.err;
  ASSERT_FALSE(r.info.is_null()) << "no run_info.json in " << dir;
  const InputParameters p = parse_input(input);
  const double initial_energy = gas_file_energy(dir / "gas_phase_states.dat", p.v_in, p.j_in);
  EXPECT_EQ(mismatches(r.info, expected_fields(p, initial_energy)), std::vector<std::string>{});
  EXPECT_EQ(mismatches(r.info, expected_derived(p, r.info)), std::vector<std::string>{});
  EXPECT_LE(std::abs(r.info["analysis"]["z_ana"].get<double>() - p.z_analysis),
            p.dz / 2.0);  // Zana moved onto the grid
  EXPECT_EQ(r.info["program"]["name"], "spodvr");
  EXPECT_EQ(r.info["program"]["version"], SPODVR_VERSION);
  EXPECT_EQ(r.info["program"]["git"], "unknown");  // the tests pass the default ProgramInfo
  EXPECT_EQ(r.info["system_id"], p.system_id);
  EXPECT_EQ(r.info["grids"]["nz"], p.nz);
  EXPECT_EQ(r.info["grids"]["nr"], p.nr);
  EXPECT_EQ(r.info["grids"]["nx"], p.nx);
  EXPECT_EQ(r.info["grids"]["ny"], p.ny);
}

// The rotational basis lists each (J, mJ) once, within JMax and the effective MJMax, with J of the
// parity of JIn for a homonuclear molecule, and its initial index points at (JIn, MJIn).
void expect_valid_basis(const nlohmann::json& info, const InputParameters& p, int mj_max_effective) {
  const nlohmann::json& states = info["rotational_basis"]["states"];
  ASSERT_GT(states.size(), 0u);
  std::set<std::pair<int, int>> seen;
  for (const nlohmann::json& s : states) {
    const int j = s[0];
    const int mj = s[1];
    EXPECT_LE(j, p.j_max);
    EXPECT_LE(std::abs(mj), std::min(j, mj_max_effective));
    if (p.homonuclear) {
      EXPECT_EQ(j % 2, p.j_in % 2);
    }
    EXPECT_TRUE(seen.insert({j, mj}).second) << "duplicate state " << j << " " << mj;
  }
  const int initial = info["rotational_basis"]["initial_index"];
  EXPECT_EQ(states[initial][0], p.j_in);
  EXPECT_EQ(states[initial][1], p.mj_in);
}

// The grid points an absorber lists lie inside its range, and the neighbouring points outside it
// (points x_i = start + i * spacing, 0 <= i < size; 1e-9 bohr for the program's tolerant comparisons).
void expect_points_inside(const nlohmann::json& absorber, double grid_start, double spacing, int size) {
  const int first = absorber["points"]["first"];
  const int last = absorber["points"]["last"];
  const double start = absorber["start"];
  const double end = absorber["end"];
  const auto x = [&](int i) { return grid_start + i * spacing; };
  EXPECT_GE(x(first), start - 1e-9) << absorber["name"];
  EXPECT_LE(x(last), end + 1e-9) << absorber["name"];
  if (first > 0) {
    EXPECT_LT(x(first - 1), start + 1e-9) << absorber["name"];
  }
  if (last + 1 < size) {
    EXPECT_GT(x(last + 1), end - 1e-9) << absorber["name"];
  }
}

nlohmann::json potential_only(nlohmann::json json) {
  json["Potential"] = true;
  return json;
}

}  // namespace

// A potential-only run of hardwall_small: the unit-carrying and derived fields against the input, the
// basis, the analysis energies (EMinWP + EParallel … EMaxWP + EParallel; total − collision = E(vin,
// jin)), the packet width and the three absorbers as given, with their grid points.
TEST(RunInfo, PotentialOnlyRunMatchesItsInput) {
  const TemporaryDirectory dir("spodvr_run_info_potential");
  const nlohmann::json input = potential_only(reference_input("hardwall_small"));
  const ProgramRun r = run(input, dir.path());
  expect_matches_input(r, input, dir.path());
  const InputParameters p = parse_input(input);
  EXPECT_EQ(r.info["run_mode"]["code"], 6);
  EXPECT_EQ(r.info["run_mode"]["name"], "potential only");
  expect_valid_basis(r.info, p, p.mj_max < 1 ? p.j_max : p.mj_max);

  const std::vector<double> collision = r.info["analysis_energies"]["collision_ev"];
  const std::vector<double> total = r.info["analysis_energies"]["total_ev"];
  ASSERT_EQ(collision.size(), static_cast<std::size_t>(p.n_energies));
  EXPECT_NEAR(collision.front(), p.e_min_wave_packet + p.e_parallel, 1e-12);
  EXPECT_NEAR(collision.back(), p.e_max_wave_packet + p.e_parallel, 1e-12);
  const double initial_ev = r.info["initial_state"]["energy_ev"];
  for (std::size_t i = 0; i < total.size(); ++i) EXPECT_NEAR(total[i] - collision[i], initial_ev, 1e-12);

  const double alpha = r.info["wave_packet"]["alpha"];
  EXPECT_NEAR(r.info["wave_packet"]["sigma_z"].get<double>() * 2.0 * std::sqrt(alpha), 1.0, 1e-15);
  const nlohmann::json& absorbers = r.info["absorbers"]["list"];
  ASSERT_EQ(absorbers.size(), 3u);
  EXPECT_EQ(absorbers[0]["name"], "ZOpt");
  EXPECT_EQ(absorbers[0]["start"], p.z_opt_start);
  EXPECT_EQ(absorbers[0]["end"], p.z_opt_end);
  EXPECT_EQ(absorbers[1]["name"], "ROpt");
  EXPECT_EQ(absorbers[1]["start"], p.r_opt_start);
  EXPECT_EQ(absorbers[2]["name"], "ZSpOpt");
  EXPECT_EQ(absorbers[2]["end"], p.zsp_opt_end);
  for (const nlohmann::json& a : absorbers) EXPECT_GT(a["strength"].get<double>(), 0.0) << a["name"];
  expect_points_inside(absorbers[0], p.z_start, p.dz, p.nz);
  expect_points_inside(absorbers[1], p.r_start, p.dr, p.nr);
  expect_points_inside(absorbers[2], p.z_start, p.dz, p.nz_specular);
}

// Both low-Z absorbers, MJMax 0 (which means JMax) and the absorbers given as factors: an absorber given
// by its factor keeps it as its effective strength (hartree), and each lists its grid points.
TEST(RunInfo, LowZAbsorbersFactorsAndMJMaxZero) {
  const TemporaryDirectory dir("spodvr_run_info_low_z");
  nlohmann::json input = potential_only(reference_input("hardwall_short"));
  input["MJMax"] = 0;
  input.erase("ZOptEnergy");
  input["ZOptFactor"] = 0.02;
  input.erase("ROptEnergy");
  input["ROptFactor"] = 0.03;
  input["ZOptLowStart"] = -2.0;
  input["ZOptLowEnd"] = 8.0;
  input["ZOptLowEnergy"] = 0.3;
  input["ZSpOptLowStart"] = -2.0;
  input["ZSpOptLowEnd"] = 6.0;
  input["ZSpOptLowFactor"] = 0.01;
  const ProgramRun r = run(input, dir.path());
  expect_matches_input(r, input, dir.path());
  const InputParameters p = parse_input(input);
  EXPECT_EQ(r.info["grids"]["mj_max"], p.j_max);
  expect_valid_basis(r.info, p, p.j_max);

  const nlohmann::json& absorbers = r.info["absorbers"]["list"];
  ASSERT_EQ(absorbers.size(), 5u);
  EXPECT_EQ(absorbers[0]["strength"], 0.02);
  EXPECT_EQ(absorbers[1]["strength"], 0.03);
  EXPECT_EQ(absorbers[3]["name"], "ZOptLow");
  EXPECT_EQ(absorbers[3]["start"], -2.0);
  EXPECT_EQ(absorbers[3]["end"], 8.0);
  EXPECT_EQ(absorbers[4]["name"], "ZSpOptLow");
  EXPECT_EQ(absorbers[4]["end"], 6.0);
  EXPECT_EQ(absorbers[4]["strength"], 0.01);
  expect_points_inside(absorbers[0], p.z_start, p.dz, p.nz);
  expect_points_inside(absorbers[1], p.r_start, p.dr, p.nr);
  expect_points_inside(absorbers[3], p.z_start, p.dz, p.nz);
  expect_points_inside(absorbers[4], p.z_start, p.dz, p.nz_specular);
}

// Oblique incidence on the 60-degree Cu(111) cell (EParallel 0.1 eV at 17 degrees, off every mirror line):
// kx and ky are the components along the reciprocal axes, which differ from the Cartesian ones here; and
// with switching z_asymptotic is ZmaxSwitch, not the PES's own value.
TEST(RunInfo, ObliqueIncidenceOnTheSkewedCell) {
  const TemporaryDirectory dir("spodvr_run_info_oblique");
  nlohmann::json input = potential_only(reference_input("cu111_dailight_small"));
  input["EParallel"] = 0.1;
  input["PhiParallel"] = 17.0;
  const ProgramRun r = run(input, dir.path());
  expect_matches_input(r, input, dir.path());
  const double skew = r.info["lattice"]["skew_angle_rad"];
  ASSERT_NEAR(skew, std::acos(-1.0) / 3.0, 1e-12);
  const double k = r.info["wave_packet"]["k_parallel"];
  const double phi = 17.0 * std::acos(-1.0) / 180.0;
  const double kx_cartesian = k * std::cos(phi);
  const double ky_cartesian = k * std::sin(phi);
  EXPECT_NEAR(r.info["wave_packet"]["kx"].get<double>(), kx_cartesian / std::sin(skew), 1e-12);
  EXPECT_NEAR(r.info["wave_packet"]["ky"].get<double>(),
              (std::cos(skew) * kx_cartesian + std::sin(skew) * ky_cartesian) / std::sin(skew), 1e-12);
  EXPECT_GT(std::abs(r.info["wave_packet"]["kx"].get<double>() - kx_cartesian), 1e-3);  // the test can tell them apart
  EXPECT_EQ(r.info["z_asymptotic"], 8.0);  // ZmaxSwitch
  EXPECT_NE(r.info["z_asymptotic_unswitched"], r.info["z_asymptotic"]);
}

// A new run, its restart and an analysis-only run of the first run's BK file each write run_info.json
// with their run mode; the restart and the analysis name their RestartDir.
TEST(RunInfo, EveryRunModeWritesIt) {
  const TemporaryDirectory dir("spodvr_run_info_modes");
  nlohmann::json first = reference_input("hardwall_short");
  first["tStart"] = 0.0;
  first["tEnd"] = 60.0;
  first["nStep"] = 30;
  const ProgramRun new_run = run(first, dir.path() / "first");
  expect_matches_input(new_run, first, dir.path() / "first");
  EXPECT_EQ(new_run.info["run_mode"]["name"], "new");
  EXPECT_TRUE(new_run.info["restart_dir"].is_null());
  EXPECT_EQ(new_run.info["new_calculation"], true);

  nlohmann::json second = reference_input("hardwall_short");
  second["tStart"] = 60.0;
  second["tEnd"] = 120.0;
  second["nStep"] = 30;
  second["NewCalc"] = false;
  second["RestartDir"] = (dir.path() / "first").string();
  const ProgramRun restart = run(second, dir.path() / "second");
  expect_matches_input(restart, second, dir.path() / "second");
  EXPECT_EQ(restart.info["run_mode"]["code"], 1);
  EXPECT_EQ(restart.info["run_mode"]["name"], "restart");
  EXPECT_EQ(restart.info["restart_dir"], (dir.path() / "first").string());
  EXPECT_EQ(restart.info["new_calculation"], false);

  nlohmann::json analysis = first;
  analysis["onlyAnalysis"] = true;
  analysis["RestartDir"] = (dir.path() / "first").string();
  const ProgramRun only = run(analysis, dir.path() / "analysis", "fine_analysis/run_info.json");
  expect_matches_input(only, analysis, dir.path() / "analysis");
  EXPECT_EQ(only.info["run_mode"]["code"], 2);
  EXPECT_EQ(only.info["run_mode"]["name"], "analysis only");
  EXPECT_EQ(only.info["restart_dir"], (dir.path() / "first").string());
  EXPECT_FALSE(fs::exists(dir.path() / "analysis" / "run_info.json"));

  // In the propagation's own directory, the analysis leaves the propagation's run_info.json alone.
  const ProgramRun in_place = run(analysis, dir.path() / "first", "fine_analysis/run_info.json");
  EXPECT_EQ(in_place.info["run_mode"]["name"], "analysis only");
  const nlohmann::json propagation = nlohmann::json::parse(std::ifstream(dir.path() / "first" / "run_info.json"));
  EXPECT_EQ(propagation["run_mode"]["name"], "new");
}

// Without a gas phase (a write2DPES run, src/app/run.cpp write_pes_cut_run) the fields that need it are null;
// the program's version and commit are written as given.
TEST(RunInfo, WithoutGasPhaseItsFieldsAreNull) {
  const TemporaryDirectory dir("spodvr_run_info_no_gas");
  const Config config = build_config(parse_input(reference_input("hardwall_small")));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  write_run_info(dir.path() / "run_info.json", config, grids, angular, ProgramInfo{"9.9.9", "abc1234-dirty"},
                 "input.json", std::nullopt);
  const nlohmann::json info = nlohmann::json::parse(std::ifstream(dir.path() / "run_info.json"));
  EXPECT_TRUE(info["initial_state"]["energy"].is_null());
  EXPECT_TRUE(info["initial_state"]["energy_ev"].is_null());
  EXPECT_TRUE(info["analysis_energies"]["collision_ev"].is_null());
  EXPECT_EQ(info["program"]["version"], "9.9.9");
  EXPECT_EQ(info["program"]["git"], "abc1234-dirty");
  EXPECT_EQ(info["grids"]["nz"], config.input.nz);
}

// Negative test of the field check: an eV value in a hartree field, NaN, and a missing field are
// found, and exactly those.
TEST(RunInfo, FieldCheckFindsUnitMistakes) {
  const TemporaryDirectory dir("spodvr_run_info_negative");
  const nlohmann::json input = potential_only(reference_input("hardwall_small"));
  const ProgramRun r = run(input, dir.path());
  ASSERT_EQ(r.code, 0) << r.err;
  const InputParameters p = parse_input(input);
  const Expected expected = expected_fields(p, gas_file_energy(dir.path() / "gas_phase_states.dat", p.v_in, p.j_in));
  ASSERT_EQ(mismatches(r.info, expected), std::vector<std::string>{});

  nlohmann::json wrong = r.info;
  wrong["potential"]["v_max"] = p.v_max;  // written in eV instead of hartree
  wrong["initial_state"]["energy"] = std::nan("");  // NaN
  wrong["molecule"].erase("reduced_mass");  // missing
  EXPECT_EQ(mismatches(wrong, expected),
            (std::vector<std::string>{"/potential/v_max", "/molecule/reduced_mass", "/initial_state/energy"}));
}
