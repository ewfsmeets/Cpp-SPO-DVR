// Unit tests of the screen formatting (report/format, report/propagation_report, report/sections):
// durations, the live progress line (only with colour, so that logs stay plain text) and the
// absorber lines of the set-up printout.
#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "report/format.hpp"
#include "report/propagation_report.hpp"
#include "report/sections.hpp"

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

}  // namespace
}  // namespace spodvr
