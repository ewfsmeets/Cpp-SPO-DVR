// Gate 5: the fine analysis (run mode 2) against the Fortran code. The C++ analyses the Fortran's
// own BK coefficients (hardwall_small/fortran_fort21.bin), so the fine analysis is tested in
// isolation. Reference data: tools/fortran_reference/run_fine_analysis_reference.sh (unmodified
// Fortran binary, onlyAnalysis runs), in hardwall_small/fine_<variant>/.
// The Fortran interpolated S-matrix elements (fort.25) are compared with an emulation of the
// Fortran interpolation including its bugs (docs/fortran_fine_analysis_bugs.md); everything
// else with the C++ results directly. Every value of fil10, unit 3 and fort.25 is compared, or
// skipped with a named, counted reason; the counts are asserted against the Fortran's write rules.
#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "app/fine_analysis_run.hpp"
#include "comparison.hpp"
#include "core/constants.hpp"
#include "core/real_compare.hpp"
#include "fft/fftw_fft.hpp"
#include "fine_analysis_reference.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"

namespace fs = std::filesystem;
using namespace spodvr;
using namespace spodvr::testing;
using cplx = std::complex<double>;

namespace {

fs::path case_dir(const std::string& case_name) { return fs::path(SPODVR_REFERENCE_DIR) / case_name; }

// The input of a Fortran reference variant (run_fine_analysis_reference.sh) of a case.
nlohmann::json variant_input(const std::string& case_name, const std::string& variant) {
  std::ifstream file(case_dir(case_name) / "input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["onlyAnalysis"] = true;
  json["writeSMatrix"] = true;
  json["RestartDir"] = case_dir(case_name).string();
  json["RestartFormat"] = "fortran";
  if (variant == "truncated") {
    json["tEnd"] = 13860.0;
    json["nStep"] = 6930;
  }
  if (variant == "small") {
    json["tEnd"] = 800.0;
    json["nStep"] = 400;
  }
  if (variant == "empty") {
    json["tEnd"] = 100.0;
    json["nStep"] = 50;
  }
  return json;
}

// Setup of an analysis-only run (no 6D potential) and its fine analysis.
struct FineRun {
  explicit FineRun(const nlohmann::json& json)
      : pes(make_potential(parse_input(json))),
        config(build_config(parse_input(json), *pes)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        gas(compute_gas_phase(config, grids, angular, compute_asymptotic_potential(*pes, config, grids), fft,
                              LapackEigensolver())),
        setup(config, angular, gas) {}

  [[nodiscard]] FineAnalysisResults run() const {
    const fs::path bk_file = fs::path(config.input.restart_dir) / "fortran_fort21.bin";
    // No checkpoint: every phase refers to Zana, as in the Fortran's fort.25 (docs/fortran_elastic_phase_bug.md).
    return run_fine_analysis(config, setup, gas, BkInput{bk_file, true}, fft);
  }
  // State index is of the labels (order, n, m, j, mJ); -1 if not in the basis.
  [[nodiscard]] int state(int order, int n, int m, int j, int mj) const {
    const AnalysisBasis& b = setup.basis;
    for (int is = 0; is < b.n_states(); ++is) {
      const DiffractionChannel& c = b.channel(is);
      if (c.order == order && c.n == n && c.m == m && b.j(is) == j && b.mj(is) == mj) return is;
    }
    return -1;
  }
  [[nodiscard]] int rotation(int j, int mj) const {
    const AnalysisBasis& b = setup.basis;
    for (int irot = 0; irot < b.n_rotation(); ++irot)
      if (b.j(irot * b.n_diffraction()) == j && b.mj(irot * b.n_diffraction()) == mj) return irot;
    return -1;
  }

  FftwFft fft;
  std::unique_ptr<Potential> pes;  // built once, as the program does
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  GasPhase gas;
  AnalysisSetup setup;
};

// Counts and worst deviation of the printed values compared.
struct Comparison {
  int n = 0;
  double worst = 0.0;

  void check(double value, const std::string& text, const std::string& what, double floor = 1e-20) {
    const double tolerance = printed_text_tolerance(text, floor);
    const double deviation = std::abs(value - parse_fortran_real(text)) / tolerance;
    EXPECT_LE(deviation, 1.0) << what << ": C++ " << value << ", Fortran " << text;
    worst = nan_aware_max(worst, deviation);
    n = n + 1;
  }
  // Phases: modulo 2 pi.
  void check_phase(double value, const std::string& text, const std::string& what) {
    const double f = parse_fortran_real(text);
    double d = value - f;
    d = d - 2.0 * Constants::pi * std::round(d / (2.0 * Constants::pi));
    check(f + d, text, what);
  }
  void report(const std::string& what) const {
    std::cout << "  " << what << ": " << n << " values, largest deviation " << worst << " x (half printed unit)\n";
  }
};

// The Fortran interpolation of S between two grid points, literally (Analyzor.f90 smatrix_calc,
// nread = 2), as in the reference (5.24 since svn r1713; the same up to r1727), including its remaining bugs: ATAN
// instead of ATAN2 for the phases, and the signs of Re S and Im S taken separately from the
// linearly interpolated S. r1713 fixed the third one, s1_real for s1_imag in the interpolated Im S
// (5.19; the unit tests in test_fine_analysis.cpp keep that documented 5.19 model). Test use only.
cplx fortran_interpolation_with_bugs(cplx s1, cplx s2, double f) {
  const double s1_real = s1.real();
  const double s1_imag = s1.imag();
  const double s2_real = s2.real();
  const double s2_imag = s2.imag();
  const double prob1 = s1_real * s1_real + s1_imag * s1_imag;
  const double prob2 = s2_real * s2_real + s2_imag * s2_imag;
  const double prob_interpolated = prob1 + (prob2 - prob1) * f;
  const double phase1 = fortran_compare::rne(s1_real, 0.0) ? std::atan(s1_imag / s1_real) : 0.0;
  const double phase2 = fortran_compare::rne(s2_real, 0.0) ? std::atan(s2_imag / s2_real) : 0.0;
  const double phase_interpolated = phase1 + (phase2 - phase1) * f;
  const double sreal_interpolated = s1_real + (s2_real - s1_real) * f;
  const double simag_interpolated = s1_imag + (s2_imag - s1_imag) * f;  // 5.24 r1713; 5.19 had s1_real
  const double sreal = std::sqrt(prob_interpolated) * std::cos(phase_interpolated);
  const double simag = std::sqrt(prob_interpolated) * std::sin(phase_interpolated);
  // Fortran SIGN(a, b): |a| with the sign of b (b = 0 counts as positive).
  auto sign = [](double a, double b) { return b >= 0.0 ? std::abs(a) : -std::abs(a); };
  return {sign(sreal, sreal_interpolated), sign(simag, simag_interpolated)};
}

// The Fortran files gate 5 reads for a variant: the case's BK file and the variant's outputs (an empty fine grid
// writes no unit 3 and no fort.25, and the gate stops after the AsympInfo numbers and the analysis limits).
std::vector<fs::path> reference_files(const std::string& case_name, const std::string& variant) {
  const fs::path dir = case_dir(case_name) / ("fine_" + variant);
  std::vector<fs::path> files = {case_dir(case_name) / "fortran_fort21.bin", dir / "stdout.txt"};
  if (variant != "empty")
    for (const std::string name : {"fort.10", "fort.3", "fort.25"}) files.push_back(dir / name);
  return files;
}

// The fort.25 rows whose S and phase gate 5 skips, per reason (measured; the probability of every row is compared):
// both grid elements around the energy below 1e-8 of the largest |S| (negligible), or the linear interpolation of
// Re S or Im S within 1e-10 of the largest |S| of zero (its sign, which the Fortran copies, is undecidable).
struct Fort25Skips {
  int negligible;
  int undecidable;
};
const std::map<std::string, Fort25Skips>& fort25_skips() {
  static const std::map<std::string, Fort25Skips> skips = {
      {"hardwall_small/full", {9856, 0}},
      {"hardwall_small/truncated", {9856, 0}},
      {"hardwall_small/small", {7629, 0}},
      {"cu111_dailight_small/full", {4131, 0}},
      // The undecidable rows: elements below 9e-8 (2e-5 of the largest |S|) on both grid points (2026-11-24).
      {"cu111_dailight_oblique/full", {8066, 22}},
  };
  return skips;
}

// Gate 5 for a variant of a case; references: whether a missing Fortran file fails the test or skips it.
void run_gate(const std::string& case_name, const std::string& variant, References references) {
  std::string missing;
  for (const fs::path& file : reference_files(case_name, variant))
    if (!fs::exists(file)) missing = missing + " " + file.string();
  if (!missing.empty()) {
    if (references == References::required) {
      ADD_FAILURE() << "reference files missing:" << missing;
      return;
    }
    GTEST_SKIP() << "no Fortran reference (" << missing
                 << "); run tools/fortran_reference/run_fine_analysis_reference.sh " << case_name;
  }
  const fs::path dir = case_dir(case_name) / ("fine_" + variant);
  const FineRun run(variant_input(case_name, variant));
  const FineAnalysisResults r = run.run();
  const FineAnalysisGrid& g = r.grid;
  const AnalysisBasis& basis = run.setup.basis;
  const int n_states = basis.n_states();
  const int n_rot = basis.n_rotation();
  const int n_orders = basis.max_order() + 1;

  // 1. Asymptotic info and the analysis limits.
  const FortranAsympInfo info = read_fortran_asymp_info(dir / "stdout.txt");
  EXPECT_EQ(g.n_records_available, info.n_records);
  EXPECT_EQ(g.n_records_used, info.n_used);
  EXPECT_EQ(g.size(), info.nen95);
  Comparison c;
  c.check(g.spectral_resolution * 1000.0, info.spres_mev, "spres");
  c.check(g.meaningful_min, info.emin, "emintr");
  c.check(g.meaningful_max, info.emax, "emaxtr");
  expect_analysis_limits(compare_analysis_limits(analysis_limits(basis),
                                                 read_fortran_analysis_limits(dir / "stdout.txt"), limits_of_a_run),
                         limits_of_a_run);
  if (g.size() == 0) {
    EXPECT_TRUE(r.interpolated.energies.empty());
    return;
  }

  // 2. probaf table: reaction probability and phase of the initial channel.
  const int initial = basis.initial_state() + n_states * (run.config.derived.v_in_index - 1);
  Comparison table;
  for (const auto& [idat, row] : info.table) {
    const EnergyResult& e = r.spectrum.energies.at(idat - 1);
    const auto& [ekin, preact, phase] = row;
    table.check(e.z_energy, ekin, "probaf ekin");
    table.check(e.p_reaction, preact, "probaf preact");
    table.check_phase(e.phase[initial], phase, "probaf phase");
  }
  EXPECT_EQ(static_cast<int>(info.table.size()), g.size());
  table.report("probaf table");

  // 3. fil10.
  const FortranFil10 fil10 = read_fortran_fil10(dir / "fort.10", run.config.input.write_rd_probabilities);
  ASSERT_EQ(fil10.nen95, g.size());
  // Every label of fil10 names a C++ level, rotational state and diffraction order.
  const auto level = [&](int v) { return v >= 0 && v < basis.n_v(); };
  for (const auto& [v, values] : fil10.p_v) ASSERT_TRUE(level(v)) << "fil10: no C++ level v = " << v;
  for (const auto& [v, values] : fil10.e_dif_v) ASSERT_TRUE(level(v)) << "fil10: no C++ level v = " << v;
  for (const auto& [key, values] : fil10.p_rovib) {
    const auto& [v, j, mj] = key;
    ASSERT_TRUE(level(v) && run.rotation(j, mj) >= 0)
        << "fil10: no C++ state (v, j, mJ) = (" << v << ", " << j << ", " << mj << ")";
  }
  for (const auto& [key, values] : fil10.p_order) {
    const auto& [order, v, j, mj] = key;
    ASSERT_TRUE(level(v) && run.rotation(j, mj) >= 0 && order >= 0 && order < n_orders)
        << "fil10: no C++ state (order, v, j, mJ) = (" << order << ", " << v << ", " << j << ", " << mj << ")";
  }
  Comparison f10;
  for (int i = 0; i < g.size(); ++i) {
    const EnergyResult& e = r.spectrum.energies[i];
    const std::string at = "E_z = " + fil10.energies[i];
    f10.check(e.z_energy, fil10.energies[i], at + " energy");
    f10.check(e.probability[initial], fil10.p_initial[i], at + " P(initial)");
    f10.check_phase(e.phase[initial], fil10.phase_initial[i], at + " phase(initial)");
    f10.check(e.p_reaction, fil10.p_reaction[i], at + " P_reaction");
    f10.check(e.e_vibration, fil10.e_vib[i], at + " <E_vib>");
    f10.check(e.e_rotation, fil10.e_rot[i], at + " <E_rot>");
    f10.check(e.e_diffraction, fil10.e_dif[i], at + " <E_dif>");
    for (const auto& [v, values] : fil10.p_v) f10.check(e.p_vibration[v], values[i], at + " P(v)");
    for (const auto& [v, values] : fil10.e_dif_v) f10.check(e.e_diffraction_v[v], values[i], at + " <E_dif>(v)");
    for (const auto& [v, values] : fil10.e_rot_v) f10.check(e.e_rotation_v[v], values[i], at + " <E_rot>(v)");
    for (const auto& [key, values] : fil10.p_rovib) {
      const auto& [v, j, mj] = key;
      f10.check(e.p_rotation[run.rotation(j, mj) + n_rot * v], values[i], at + " P(v, j, mJ)");
    }
    for (const auto& [key, values] : fil10.p_order) {
      const auto& [order, v, j, mj] = key;
      f10.check(e.p_order[order + n_orders * (run.rotation(j, mj) + n_rot * v)], values[i], at + " P(order)");
    }
  }
  f10.report("fil10");

  // 4. Unit 3: every state on the fine grid.
  Comparison u3;
  for (const auto& [key, values] : read_fortran_unit3(dir / "fort.3")) {
    const auto& [order, n, m, v, j, mj] = key;
    const int is = run.state(order, n, m, j, mj);
    ASSERT_TRUE(is >= 0 && level(v)) << "unit 3: no C++ state (order, n, m, v, j, mJ) = (" << order << ", " << n << ", "
                                     << m << ", " << v << ", " << j << ", " << mj << ")";
    for (int i = 0; i < g.size(); ++i)
      u3.check(r.spectrum.energies[i].probability[is + n_states * v], values[i], "unit 3 P");
  }
  u3.report("unit 3 (all states, fine grid)");

  // 5. The interpolated per-energy blocks of the stdout (probabilities and average energies: unaffected by the bugs).
  PropagationRun interpolated;
  interpolated.scattering = r.interpolated;
  set_analysis_labels(interpolated, run.config, basis, run.gas);
  compare_energy_blocks_with_fortran(interpolated, read_fortran_energy_blocks(dir / "stdout.txt"));
  EXPECT_EQ(static_cast<int>(r.skipped.size()), info.n_skip_messages);

  // 6. fort.25: every row's probability against the C++; S and phase against the bug emulation, or skipped with
  // a reason (fort25_skips).
  const Fort25File fort25 = read_fortran_fort25_rows(dir / "fort.25");
  for (const std::string& problem : fort25.problems) ADD_FAILURE() << "fort.25 " << problem;
  ASSERT_EQ(fort25.total_energy.size(), r.interpolated.energies.size());
  EXPECT_EQ(static_cast<int>(fort25.rows.size()), fort25.row_like_lines)
      << "fort.25 rows read vs an independent count of the file";
  double s_max = 0.0;
  for (const EnergyResult& e : r.spectrum.energies)
    for (const cplx s : e.s_matrix) s_max = nan_aware_max(s_max, std::abs(s));
  Comparison p25;
  Comparison s25;
  Fort25Skips skipped{0, 0};
  for (const Fort25Row& row : fort25.rows) {
    const EnergyResult& e = r.interpolated.energies[row.block];
    const int is = run.state(row.order, row.n, row.m, row.j, row.mj);
    if (is < 0 || row.v < 0 || row.v >= basis.n_v()) {
      ADD_FAILURE() << "fort.25 line " << row.line << ": no C++ state";
      continue;
    }
    const int k = is + n_states * row.v;
    p25.check(e.probability[k], row.p, "fort.25 P");

    const auto [upper, f] = locate_on_grid(g, e.z_energy);
    const cplx s1 = r.spectrum.energies[upper - 1].s_matrix[k];
    const cplx s2 = r.spectrum.energies[upper].s_matrix[k];
    if (std::abs(s1) < 1e-8 * s_max && std::abs(s2) < 1e-8 * s_max) {
      skipped.negligible = skipped.negligible + 1;
      continue;
    }
    // Signs decided by a linear interpolation within rounding of zero cannot be compared.
    const double re_lin = s1.real() + (s2.real() - s1.real()) * f;
    // As the reference since 5.24 r1713 (5.19 had s1.real()).
    const double im_lin = s1.imag() + (s2.imag() - s1.imag()) * f;
    if (std::abs(re_lin) < 1e-10 * s_max || std::abs(im_lin) < 1e-10 * s_max) {
      skipped.undecidable = skipped.undecidable + 1;
      continue;
    }
    // The interval located by the C++ is the Fortran's on every gate-5 case: a fallback to the neighbouring interval
    // for an energy on a grid point (the Fortran's E_z is rounded from quad precision) matched better 0 times
    // (measured 2026-11-24) and was removed; NrOfEs = -95 (energies on the grid) is an input error since r1719.
    const cplx emulated = fortran_interpolation_with_bugs(s1, s2, f);
    s25.check(emulated.real(), row.re, "fort.25 Re S (emulated)");
    s25.check(emulated.imag(), row.im, "fort.25 Im S (emulated)");
    const double phase = emulated.real() == 0.0 ? 0.0 : std::atan2(emulated.imag(), emulated.real());
    s25.check(phase, row.phase, "fort.25 phase (emulated)");
  }
  p25.report("fort.25 probabilities");
  s25.report("fort.25 S-matrix and phases (Fortran interpolation emulated)");
  std::cout << "  fort.25: " << fort25.rows.size() << " rows; S and phase skipped: " << skipped.negligible
            << " negligible, " << skipped.undecidable << " with an undecidable sign\n";

  // 7. The counts: every value of every output compared, every skip pinned.
  const auto pin = fort25_skips().find(case_name + "/" + variant);
  ASSERT_NE(pin, fort25_skips().end()) << "no pinned fort.25 skips for " << case_name << "/" << variant;
  EXPECT_EQ(skipped.negligible, pin->second.negligible) << "fort.25 rows skipped as negligible vs the pin";
  EXPECT_EQ(skipped.undecidable, pin->second.undecidable) << "fort.25 rows skipped as undecidable vs the pin";
  expect_fine_counts({f10.n, u3.n, p25.n, s25.n},
                     expected_fine_counts(interpolated, g.size(), skipped.negligible + skipped.undecidable));
}

}  // namespace

// The committed variants: a missing reference file fails them.
TEST(RegressionFineAnalysis, Gate5Full) { run_gate("hardwall_small", "full", References::required); }
TEST(RegressionFineAnalysis, Gate5Truncated) { run_gate("hardwall_small", "truncated", References::required); }
TEST(RegressionFineAnalysis, Gate5Small) { run_gate("hardwall_small", "small", References::required); }
TEST(RegressionFineAnalysis, Gate5Empty) { run_gate("hardwall_small", "empty", References::required); }
// The corrugated C6v case (60 deg cell): on the BK file of the full Dai-Light Fortran run.
TEST(RegressionFineAnalysis, Gate5DaiLight) { run_gate("cu111_dailight_small", "full", References::required); }
// Oblique incidence with J > 0 (cu111_dailight_oblique): on the BK file of its Fortran run.
TEST(RegressionFineAnalysis, Gate5DaiLightOblique) { run_gate("cu111_dailight_oblique", "full", References::required); }

// A required case whose reference files are missing fails (once, naming them) instead of skipping.
TEST(RegressionFineAnalysis, MissingRequiredReferenceFails) {
  EXPECT_NONFATAL_FAILURE(run_gate("no_such_case", "full", References::required), "reference files missing");
}

// A NaN on either side of a printed value fails its check, and the reported largest deviation is NaN.
TEST(RegressionFineAnalysis, PrintedValueCheckFailsOnNaN) {
  const double not_a_number = std::numeric_limits<double>::quiet_NaN();
  Comparison cpp_nan;
  EXPECT_NONFATAL_FAILURE(cpp_nan.check(not_a_number, "0.1234567E+00", "C++ value"), "C++ value");
  EXPECT_TRUE(std::isnan(cpp_nan.worst));
  Comparison fortran_nan;
  EXPECT_NONFATAL_FAILURE(fortran_nan.check(0.1234567, "NaN", "Fortran value"), "Fortran value");
  EXPECT_TRUE(std::isnan(fortran_nan.worst));
}
