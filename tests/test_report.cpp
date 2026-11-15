// Unit tests of the screen formatting (report/format, report/propagation_report, report/sections,
// report/potential_report): durations, the live progress line (only with colour, so that logs stay
// plain text), the absorber lines and warnings of the set-up printout and the SCM block.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "report/format.hpp"
#include "report/potential_report.hpp"
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

// An absorber whose L/lambda lies outside the Vibok & Balint-Kurti table gets a warning naming its
// energy key, L/lambda and the table end (the Fortran's getACoeffOptPotential warns without naming
// the absorber). ROptEnergy 0.01 eV: L/lambda = 0.58847; ZSpOptEnergy 3 eV: 15.9458.
TEST(Report, AbsorberOutsideTheStrengthTableIsWarnedAbout) {
  const std::string text = setup_printout(config_with({{"ROptEnergy", 0.01}, {"ZSpOptEnergy", 3.0}}));
  EXPECT_NE(line_with(text,
                      "WARNING: ROptEnergy: at this energy the absorber is 0.58847 wavelengths long (L/lambda), "
                      "below the Vibok-Balint-Kurti table (L/lambda = 1 .. 15), so its strength is the "
                      "table's value at L/lambda = 1 (inputs/README.md, section \"Absorbers\")"),
            "")
      << text;
  EXPECT_NE(line_with(text,
                      "WARNING: ZSpOptEnergy: at this energy the absorber is 15.9458 wavelengths long "
                      "(L/lambda), above the Vibok-Balint-Kurti table (L/lambda = 1 .. 15), so its strength is "
                      "the table's value at L/lambda = 15 (inputs/README.md, section \"Absorbers\")"),
            "")
      << text;
  EXPECT_EQ(line_with(text, "ZOptEnergy:"), "");
}

// absorber D2: at L/lambda = 1 exactly the strength is 0 and the warning says the absorber is off.
// A value just above 15 is printed with all digits, so that it does not read "15" (and one just
// below 1 so that it does not read "1"). The records are set by hand: an input that hits 1 exactly
// is impractical.
TEST(Report, AbsorberWarningAtTheFirstTablePointAndNearTheEnds) {
  Config config = config_with(nlohmann::json::object());
  config.derived.absorber_table_warnings = {{"ZOptEnergy", 1.0, AbsorberTableRange::at_first_point},
                                            {"ROptEnergy", std::nextafter(1.0, 0.0), AbsorberTableRange::below},
                                            {"ZSpOptEnergy", std::nextafter(15.0, 16.0), AbsorberTableRange::above}};
  const std::string text = setup_printout(config);
  EXPECT_NE(line_with(text,
                      "WARNING: ZOptEnergy: at this energy the absorber is exactly 1 wavelength long "
                      "(L/lambda = 1), the first point of the Vibok-Balint-Kurti table, which no interval of "
                      "the table includes: its strength is 0 and the absorber is off"),
            "")
      << text;
  EXPECT_NE(line_with(text, "WARNING: ROptEnergy: at this energy the absorber is 0.9999999999999999 wavelengths long"),
            "")
      << text;
  EXPECT_NE(
      line_with(text, "WARNING: ZSpOptEnergy: at this energy the absorber is 15.000000000000002 wavelengths long"), "")
      << text;
}

// An SCM slab of 3 + 2 atoms and its diagnostics, as diagnose_scm_potential would return them.
struct ScmReportData {
  ScmInput input;
  ScmSlab slab;
  ScmDiagnostics diagnostics;

  ScmReportData() {
    input.parameter_file = "params.dat";
    input.position_file = "slab.dat";
    slab.header = {"MD 300 K", "FCC111", 7, 7, 7, 20, 3, 0, 4.8896663475};
    slab.configuration = 4;
    slab.cutoff = 16.0;
    ScmDiagnostics& d = diagnostics;
    d.crp_lattice_constant = 4.9161534512607137;
    d.rescaling_factor = 4.8896663475 / 4.9161534512607137;
    d.slab = {3, 2, 0.25, 0.5, 0.75, {0.01, -0.02, 0.03}};
    d.lattice = {3, 2, 1.2e-6, 4e-7, 0};
    d.asymptotic.z_asymptotic = 11.9;
    d.asymptotic.point_heights = {11.9, 12.283219};
    d.asymptotic.largest_at_points.points = 96;
    d.asymptotic.largest_on_plane = ScmExtreme{4.8e-5, {1.0, 2.0, 11.9, 4.0, 0.4363323129985824, 0.0}, 20520};
    d.asymptotic.plane_r_limit = 4.2;
    d.edge_jumps = ScmEdgeJumps{-3.7, 12, 12, {{2.25, 1.95, 1.1, 194, 0.61}, {4.0, 1.4, 0.13, 0, 0.0}}};
    d.holes = ScmHoleCheck{0.472431, 0.75589, 4.34637, 0.5, 0, {}};
    d.outside_r_range = ScmOutsideRRangeCheck{
        0.75589, 4.34637, 3, 1, {-231.7, {2.4448, 0.0, 1.55, 5.8, 1.0016, 0.0}, 246240}, 29, 4.2, 9.8, true};
  }
  [[nodiscard]] std::string printout() const {
    std::ostringstream out;
    print_scm_potential(out, Style{false}, input, slab, diagnostics);
    return out.str();
  }
};

// The SCM block shows the files, the snapshot, the lattice constants with s, the displacements and
// one line per check (SCM D9, D3, D15) plus the cell-edge table (SCM D10); without a finding it has
// no warning.
TEST(Report, ScmBlockShowsTheSlabAndTheChecks) {
  const std::string text = ScmReportData().printout();
  EXPECT_NE(line_with(text, "Position file (PosFilepath)").find("slab.dat"), std::string::npos);
  EXPECT_NE(line_with(text, "Slab tag").find("MD 300 K"), std::string::npos);
  EXPECT_NE(line_with(text, "Structure, file version, snapshots").find("FCC111, 3, 20"), std::string::npos);
  EXPECT_NE(line_with(text, "Snapshot used (PosFileSlabNum)").find(" 4"), std::string::npos);
  EXPECT_NE(line_with(text, "Atoms used (layer 1 + layer 2)").find("3 + 2"), std::string::npos);
  EXPECT_NE(line_with(text, "Rescaling (EnableSurfExpan)").find("on, s = a_T/a0 = 9.946122300649"), std::string::npos);
  EXPECT_NE(line_with(text, "Highest atom (distorted z)").find("0.750000 bohr"), std::string::npos);
  EXPECT_NE(line_with(text, "(SCM D9)").find("1.20e-06 bohr at most (3 + 2 atoms checked)"), std::string::npos);
  EXPECT_NE(
      line_with(text, "Gas-phase points with V != V_CRP (SCM D3)").find("0 of 96 (Z = 11.900000, 12.283219 bohr)"),
      std::string::npos);
  EXPECT_NE(line_with(text, "Largest |V - V_CRP| on Z = 11.9000, r < 4.2000")
                .find("4.800e-05 eV at (u, v) = (1.0000, 2.0000) bohr, r = 4.0000 bohr, theta = 25.0, phi = 0.0 deg "
                      "(20520 points)"),
            std::string::npos);
  EXPECT_NE(line_with(text, "(SCM D15)").find("none (ZStart 0.500000 >= Zmin 0.472431 bohr)"), std::string::npos);
  EXPECT_NE(line_with(text, "2.25      1.95").find("1.100000        194 of 288                  0.610000"),
            std::string::npos);
  EXPECT_NE(line_with(text, "4.00      1.40").find("0 of 288                         -"), std::string::npos);
  EXPECT_NE(line_with(text, "r outside Rmin .. Rmax, every Z plane (SCM D15)")
                .find("Rmin .. Rmax = 0.755890 .. 4.346370 bohr: lowest V - V_CRP -2.317e+02 eV at Z = 1.5500 bohr, "
                      "(u, v) = (2.4448, 0.0000) bohr, r = 5.8000 bohr, theta = 57.4, phi = 0.0 deg, beyond Rmax, "
                      "inside the r absorber (r = 4.2000 .. 9.8000 bohr); 246240 points: 3 theta x 1 phi of the grid"),
            std::string::npos);
  EXPECT_NE(line_with(text, "(min V_asy + EMaxWP + 0.5 eV)"), "");  // scm_accessible_margin
  EXPECT_EQ(line_with(text, "Warning"), "");
  EXPECT_FALSE(has_control_codes(text));
}

// Each finding gives its warning: an ideal atom more than 1e-5 bohr off the lattice (SCM D9), |V -
// V_CRP| above 1 meV on the plane (SCM D3), a Z grid below Zmin with its lowest V - V_CRP (SCM D15).
// An analysis-only run shows no cell-edge table.
TEST(Report, ScmBlockWarnings) {
  ScmReportData data;
  data.diagnostics.lattice.top_layer_deviation = 2e-5;
  data.diagnostics.asymptotic.largest_on_plane->value = 2e-3;
  data.diagnostics.holes =
      ScmHoleCheck{0.472431, 0.75589, 4.34637, -0.55, 7, {-257.5, {0.0, 0.0, -0.55, 3.0, 1.0, 0.5}, 9000}};
  const std::string text = data.printout();
  EXPECT_NE(line_with(text, "Warning: an ideal atom lies 2.00e-05 bohr from its site of the CRP lattice"), "");
  EXPECT_NE(line_with(text, "Warning: |V - V_CRP| reaches 2.000e-03 eV on the plane Z = z_asymptotic"), "");
  EXPECT_NE(line_with(text, "Warning: the Z grid starts at ZStart = -0.550000 bohr, below the CRP's Zmin = 0.472431"),
            "");
  EXPECT_NE(line_with(text, "the SCM term at the grid's Z: V - V_CRP reaches -2.575e+02 eV at Rmin <= r <= Rmax"), "");
  EXPECT_NE(line_with(text, "(SCM D15). Use ZStart >= Zmin, or ImposeVMin."), "");
  EXPECT_NE(
      line_with(text, "(SCM D15)").find("7, lowest V - V_CRP at Rmin <= r <= Rmax: -2.575e+02 eV at Z = -0.5500 bohr"),
      std::string::npos);

  // Outside the r range only a line, never a warning: beyond Rmax outside the r absorber, with no r
  // grid point in the absorber, below Rmin, and with no r outside the range.
  ScmOutsideRRangeCheck& outside_range = *data.diagnostics.outside_r_range;
  outside_range.in_r_absorber = false;
  EXPECT_NE(
      line_with(data.printout(), "every Z plane (SCM D15)").find(", beyond Rmax, outside the r absorber (r = 4.2"),
      std::string::npos);
  const std::string outside = data.printout();
  std::size_t warnings = 0;
  for (std::size_t at = outside.find("Warning:"); at != std::string::npos; at = outside.find("Warning:", at + 1))
    ++warnings;
  EXPECT_EQ(warnings, 3u);  // SCM D9, D3 and D15 below Zmin; none for the hole outside the r range
  outside_range.r_absorber_points = 0;
  EXPECT_NE(line_with(data.printout(), "every Z plane (SCM D15)")
                .find(", beyond Rmax; no r grid point lies in the r absorber;"),
            std::string::npos);
  outside_range.lowest.at.r = 0.4;
  EXPECT_NE(line_with(data.printout(), "every Z plane (SCM D15)")
                .find("r = 0.4000 bohr, theta = 57.4, phi = 0.0 deg, "
                      "below Rmin; 246240 points"),
            std::string::npos);
  outside_range.lowest.points = 0;
  EXPECT_NE(line_with(data.printout(), "every Z plane (SCM D15)")
                .find("Rmin .. Rmax = 0.755890 .. 4.346370 bohr: no r grid point outside"),
            std::string::npos);

  data.diagnostics.edge_jumps.reset();
  data.diagnostics.holes.reset();
  data.diagnostics.outside_r_range.reset();
  data.diagnostics.asymptotic.largest_on_plane.reset();
  const std::string analysis_only = data.printout();
  EXPECT_NE(line_with(analysis_only, "Analysis-only run: no 6D potential, so only the gas-phase points are checked."),
            "");
  EXPECT_EQ(line_with(analysis_only, "max accessible jump"), "");
  EXPECT_EQ(line_with(analysis_only, "(SCM D15)"), "");
  EXPECT_EQ(line_with(analysis_only, "Rmax"), "");
  EXPECT_EQ(line_with(analysis_only, "Largest |V - V_CRP|"), "");
}

}  // namespace
}  // namespace spodvr
