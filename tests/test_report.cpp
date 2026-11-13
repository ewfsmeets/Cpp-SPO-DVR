// Unit tests of the screen formatting (report/format, report/propagation_report, report/sections,
// report/flux_report): durations, the live progress line (only with colour, so that logs stay plain
// text), the absorber lines of the set-up printout, and the flux analysis's set-up block, result table
// and data files (report/flux_files).
#include <gtest/gtest.h>
#include <unistd.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "core/unit_conversion.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "report/flux_files.hpp"
#include "report/flux_report.hpp"
#include "report/format.hpp"
#include "report/propagation_report.hpp"
#include "report/sections.hpp"
#include "symmetry/surface_lattice.hpp"

namespace spodvr {
namespace {

bool has_control_codes(const std::string& text) {
  return text.find('\r') != std::string::npos || text.find('\033') != std::string::npos;
}

// The hard wall reference case with changed keys.
Config config_with(const nlohmann::json& changes) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json.update(changes);
  return build_config(parse_input(json));
}

// The set-up printout of the input summary and the grids, without colour.
std::string setup_printout(const Config& config) {
  std::ostringstream out;
  const Style style{false};
  print_input_summary(out, style, config);
  print_grids(out, style, config, make_spatial_grids(config), make_angular_grid(config));
  return out.str();
}

// The printed line that contains `label` ("" if there is none).
std::string line_with(const std::string& text, const std::string& label) {
  std::istringstream lines(text);
  for (std::string line; std::getline(lines, line);)
    if (line.find(label) != std::string::npos) return line;
  return "";
}

TEST(Format, DurationUnderAnHourIsMinutesAndSeconds) {
  EXPECT_EQ(format_duration(0.0), "0:00");
  EXPECT_EQ(format_duration(59.4), "0:59");
  EXPECT_EQ(format_duration(61.0), "1:01");
}

TEST(Format, DurationFromAnHourHasHours) {
  EXPECT_EQ(format_duration(3600.0), "1:00:00");
  EXPECT_EQ(format_duration(5532.1), "1:32:12");
}

TEST(AveragesTable, WithoutColourThereIsNoProgressLineAndNoControlCodes) {
  std::ostringstream out;
  const Style style{false};
  AveragesTable table(out, style, 8, 2.5, 0.0);
  table.print_header();
  table.print_row(Averages{0.0, 1.0, 17.0, 1.4, 1.0});
  table.print_status(10.0);
  table.print_note("analysis surface fixed", 10.0);
  table.finish();
  const std::string text = out.str();
  EXPECT_FALSE(has_control_codes(text));
  EXPECT_EQ(text.find("step 4/8"), std::string::npos);
  EXPECT_NE(text.find("Propagation time: "), std::string::npos);
}

TEST(AveragesTable, WithColourTheProgressLineIsDrawnAndClearedBeforeTheNextRow) {
  std::ostringstream out;
  const Style style{true};
  AveragesTable table(out, style, 8, 2.5, 0.0);
  table.print_status(10.0);
  const std::string status = out.str();
  EXPECT_EQ(status.rfind("\r\033[K", 0), 0u);
  EXPECT_NE(status.find("step 4/8"), std::string::npos);
  EXPECT_NE(status.find(" 50%"), std::string::npos);
  table.print_row(Averages{10.0, 1.0, 17.0, 1.4, 1.0});
  EXPECT_EQ(out.str().substr(status.size(), 4), "\r\033[K");
}

// low-Z D4: the set-up shows each low-Z absorber's own range, strength and points (the Fortran's
// specular block prints the scattering grid's values instead), and nothing for an absorber that is
// off. The grid-point trap gives a warning; an input without low-Z keys gets neither.
TEST(Report, LowZD4PrintsEachAbsorbersOwnValues) {
  const Config specular = config_with({{"ZSpOptLowStart", -2.0}, {"ZSpOptLowEnd", 6.0}, {"ZSpOptLowEnergy", 0.3}});
  const std::string text = setup_printout(specular);
  const std::string cap = line_with(text, "Low-Z CAP in Z (specular)");
  EXPECT_NE(cap.find("-2.000000 .. 6.000000 bohr, strength " + scientific(specular.derived.zsp_opt_low_factor, 8) +
                     " hartree (energy 0.300000 eV)"),
            std::string::npos)
      << cap;
  const std::string points = line_with(text, "Low-Z optical potential Z specular (indices)");
  EXPECT_NE(points.find(" 0 .. 100"), std::string::npos) << points;  // Z = -2.0 .. 6.0
  EXPECT_EQ(line_with(text, "Low-Z CAP in Z (scattering)"), "");
  EXPECT_EQ(line_with(text, "Low-Z optical potential Z (indices)"), "");
  EXPECT_EQ(line_with(text, "WARNING"), "");

  const std::string trap =
      setup_printout(config_with({{"ZOptLowStart", 0.4}, {"ZOptLowEnd", 8.0}, {"ZOptLowEnergy", 0.3}}));
  EXPECT_NE(line_with(trap,
                      "WARNING: ZOptLowStart = 0.4 bohr: grid point 30 is 0.39999999999999991 here and "
                      "0.40000000000000002 in the Fortran (quad precision), so it is outside the absorber here "
                      "and inside it in the Fortran"),
            "");

  const std::string plain = setup_printout(config_with(nlohmann::json::object()));
  EXPECT_EQ(line_with(plain, "Low-Z"), "");
  EXPECT_EQ(line_with(plain, "WARNING"), "");
}

// A reference case's input with the flux analysis at r_c = 2.0 bohr (r index 8), 0.2 .. 0.6 eV by
// 0.15 eV (3 energies, flux Q3).
Config flux_config(const std::string& case_name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + case_name + "/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json.update({{"doFluxAnalysis", true},
               {"rFluxAnalysis", 2.0},
               {"EMinFluxAna", 0.2},
               {"EMaxFluxAna", 0.6},
               {"DeltaEFluxAna", 0.15}});
  return build_config(parse_input(json));
}

// Flux results of n_energies energies on an nx x ny grid with values that need all 17 digits (thirds,
// negative, tiny), distinct per (ix, iy, n), so that a wrong order or a lost digit shows.
FluxResults made_up_flux_results(int n_energies, int nx, int ny) {
  FluxResults r;
  r.density_index = Index3D(nx, ny, n_energies);
  for (int n = 0; n < n_energies; ++n) {
    r.energy.push_back(0.2 + n / 3.0);
    r.probability.push_back((n % 2 == 0 ? 1.0 : -1e-17) / (n + 3.0));
  }
  r.density.resize(r.density_index.size());
  for (std::size_t k = 0; k < r.density.size(); ++k) r.density[k] = (static_cast<double>(k) + 0.1) / 7.0;
  return r;
}

// The lines of a file (or a printout).
std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> lines;
  std::istringstream in(text);
  for (std::string line; std::getline(in, line);) lines.push_back(line);
  return lines;
}

std::string file_text(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::ostringstream text;
  text << file.rdbuf();
  return text.str();
}

// The tab-separated cells of a line.
std::vector<std::string> cells_of(const std::string& line) {
  std::vector<std::string> cells;
  std::istringstream in(line);
  for (std::string cell; std::getline(in, cell, '\t');) cells.push_back(cell);
  return cells;
}

// The rows of a TSV file after its '#' comments: the header first.
std::vector<std::vector<std::string>> table_of(const std::filesystem::path& path) {
  std::vector<std::vector<std::string>> rows;
  for (const std::string& line : lines_of(file_text(path)))
    if (!line.starts_with("#")) rows.push_back(cells_of(line));
  return rows;
}

double number(const std::string& cell) { return std::strtod(cell.c_str(), nullptr); }

// A scratch directory for one test, removed with it.
struct ScratchDirectory {
  std::filesystem::path path;
  explicit ScratchDirectory(const std::string& name)
      : path(std::filesystem::temp_directory_path() / (name + "_" + std::to_string(::getpid()))) {
    std::filesystem::remove_all(path);
    std::filesystem::create_directories(path);
  }
  ~ScratchDirectory() { std::filesystem::remove_all(path); }
};

// The three flux files of a skewed (60 deg) grid: the header lines, one row per energy and one per
// (energy, ix, iy) in the Fortran's order (energy outer, iy fastest), every number read back bitwise
// (%.16e), X and Y the Cartesian position of (x, y), and the JSON metadata.
TEST(FluxFiles, TablesHoldEveryValueBitwiseInTheFortranOrder) {
  const Config config = flux_config("cu111_dailight_small");
  const SpatialGrids grids = make_spatial_grids(config);
  const int nx = grids.x.size();
  const int ny = grids.y.size();
  const FluxResults results = made_up_flux_results(3, nx, ny);
  const ScratchDirectory dir("spodvr_flux_files");
  const std::vector<std::filesystem::path> files = write_flux_files(dir.path, config, grids, results);
  ASSERT_EQ(files.size(), 3U);

  const auto reaction = table_of(dir.path / "flux_reaction.tsv");
  ASSERT_EQ(reaction.size(), 4U);
  EXPECT_EQ(reaction[0], (std::vector<std::string>{"E_z", "P_flux"}));
  for (int n = 0; n < 3; ++n) {
    ASSERT_EQ(reaction[n + 1].size(), 2U);
    EXPECT_EQ(number(reaction[n + 1][0]), results.energy[n]);
    EXPECT_EQ(number(reaction[n + 1][1]), results.probability[n]);
  }

  const auto density = table_of(dir.path / "flux_density.tsv");
  ASSERT_EQ(density.size(), std::size_t(1 + 3 * nx * ny));
  EXPECT_EQ(density[0], (std::vector<std::string>{"E_z", "ix", "iy", "x", "y", "X", "Y", "density"}));
  const SurfaceLattice lattice(config.derived.system);
  std::size_t row = 1;
  for (int n = 0; n < 3; ++n)
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy, ++row) {
        const std::vector<std::string>& c = density[row];
        ASSERT_EQ(c.size(), 8U) << "row " << row;
        const CartesianVector cartesian = lattice.to_cartesian(grids.x[ix], grids.y[iy]);
        EXPECT_EQ(number(c[0]), results.energy[n]);
        EXPECT_EQ(c[1], std::to_string(ix));
        EXPECT_EQ(c[2], std::to_string(iy));
        EXPECT_EQ(number(c[3]), grids.x[ix]);
        EXPECT_EQ(number(c[4]), grids.y[iy]);
        EXPECT_EQ(number(c[5]), cartesian.x);
        EXPECT_EQ(number(c[6]), cartesian.y);
        EXPECT_EQ(number(c[7]), results.density[results.density_index(ix, iy, n)]) << "row " << row;
      }
  // The skewed grid: X differs from x once y > 0 (cos 60 deg).
  EXPECT_NE(number(density[2][5]), number(density[2][3]));

  const nlohmann::json json = nlohmann::json::parse(file_text(dir.path / "flux_analysis.json"));
  EXPECT_EQ(json["files"], (nlohmann::json{"flux_reaction.tsv", "flux_density.tsv", "flux_analysis.json"}));
  EXPECT_EQ(json["flux_surface"]["ir"], 8);
  EXPECT_EQ(json["flux_surface"]["r"].get<double>(), grids.r[8]);
  EXPECT_EQ(json["energies"]["count"], 3);
  EXPECT_EQ(json["surface_grid"]["nx"], nx);
  EXPECT_EQ(json["surface_grid"]["point_group"], config.derived.system.point_group);
}

// The writer refuses results of another surface grid instead of writing a misread table.
TEST(FluxFiles, ResultsOfAnotherGridAreRefused) {
  const Config config = flux_config("hardwall_small");
  const SpatialGrids grids = make_spatial_grids(config);
  const ScratchDirectory dir("spodvr_flux_files_refused");
  EXPECT_THROW((void)write_flux_files(dir.path, config, grids, made_up_flux_results(3, grids.x.size() + 1, 4)),
               std::invalid_argument);
  EXPECT_FALSE(std::filesystem::exists(dir.path / "flux_reaction.tsv"));
}

// The result rows of a printed P(E) table (the lines below its dashed rule).
int printed_rows(const std::string& text) {
  const std::vector<std::string> lines = lines_of(text);
  int rows = 0;
  bool below_rule = false;
  for (const std::string& line : lines) {
    if (below_rule && !line.empty()) ++rows;
    if (line.starts_with("---")) below_rule = true;
  }
  return rows;
}

// Like the analysis tables, P(E) shows every energy up to 50 and every k-th above (k = ceil(n / 50)),
// with a note; the first energy is always shown.
TEST(FluxReport, ResultTableIsThinnedAbove50Energies) {
  for (const auto& [n_energies, rows, stride] : {std::tuple{3, 3, 1}, {50, 50, 1}, {51, 26, 2}, {120, 40, 3}}) {
    std::ostringstream out;
    print_flux_results(out, Style{false}, made_up_flux_results(n_energies, 1, 1));
    const std::string text = out.str();
    EXPECT_EQ(printed_rows(text), rows) << n_energies << " energies";
    const std::string note = "Every " + std::to_string(stride) + "th energy is shown";
    EXPECT_EQ(text.find(note) != std::string::npos, stride > 1) << n_energies << " energies";
    EXPECT_NE(text.find(fixed(0.2, 4)), std::string::npos);
  }
}

// The FLUX ANALYSIS set-up block: r_c and its index, the energies, the memory numbers of the rule, and
// the flux Q7 warning only above 1e-6 (also for NaN); the flux Q3 warning of the energies.
TEST(FluxReport, SetupBlockShowsTheMemoryRuleAndWarnsAboveTheFluxQ7Limit) {
  const Config config = flux_config("hardwall_small");
  const SpatialGrids grids = make_spatial_grids(config);
  FluxEnergies energies;
  energies.z_energy = {units::precise.ev_to_hartree(0.2), units::precise.ev_to_hartree(0.35)};
  energies.warnings = {"the flux energies end at 0.35 eV"};
  FluxMemoryCheck memory{Device::Cpu, {100, 200}, {1536, 2, 512}, 1u << 20, true, {5, true}, 680};
  auto printout = [&](double ratio) {
    std::ostringstream out;
    print_flux_analysis_setup(out, Style{false}, config, grids, energies, ratio, memory);
    return out.str();
  };
  const std::string text = printout(1e-6);
  EXPECT_NE(line_with(text, "Flux surface r_c (index)").find("2.000000 (8)"), std::string::npos) << text;
  EXPECT_NE(line_with(text, "Flux energies (E_z)").find("0.200000 .. 0.350000"), std::string::npos) << text;
  EXPECT_NE(line_with(text, "Slots (analysis steps per flush)").find(" 5"), std::string::npos) << text;
  EXPECT_NE(line_with(text, "Memory per flux energy").find("1.5 KiB"), std::string::npos) << text;
  EXPECT_NE(line_with(text, "Most flux energies that fit").find("680 (with 1 slot)"), std::string::npos) << text;
  EXPECT_NE(line_with(text, "Warning: the flux energies end at 0.35 eV"), "");
  EXPECT_EQ(line_with(text, "flux Q7"), "");
  EXPECT_NE(line_with(printout(1.01e-6), "flux Q7"), "");
  EXPECT_NE(line_with(printout(std::numeric_limits<double>::quiet_NaN()), "flux Q7"), "");

  memory.max_energies.reset();
  memory.physical_bytes.reset();
  EXPECT_NE(line_with(printout(0.0), "Most flux energies that fit").find("unknown"), std::string::npos);
}

}  // namespace
}  // namespace spodvr
