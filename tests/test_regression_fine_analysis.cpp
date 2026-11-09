// Gate 5: the fine analysis (run mode 2) against the Fortran code. The C++ analyses the Fortran's
// own BK coefficients (hardwall_small/fortran_fort21.bin), so the fine analysis is tested in
// isolation. Reference data: tools/fortran_reference/run_fine_analysis_reference.sh (unmodified
// Fortran binary, onlyAnalysis runs), in hardwall_small/fine_<variant>/.
// The Fortran interpolated S-matrix elements (fort.25) are compared with an emulation of the
// Fortran interpolation including its bugs (docs/fortran_fine_analysis_bugs.md); everything
// else with the C++ results directly.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "app/fine_analysis_run.hpp"
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
#include "symmetry/surface_lattice.hpp"

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
  if (variant == "truncated" || variant == "nen95") {
    json["tEnd"] = 13860.0;
    json["nStep"] = 6930;
  }
  if (variant == "nen95") json["NrOfEs"] = -95;
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
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        gas(compute_gas_phase(config, grids, angular,
                              compute_asymptotic_potential(*make_potential(config.input), config, grids), fft,
                              LapackEigensolver())),
        setup(config, angular, gas) {}

  [[nodiscard]] FineAnalysisResults run() const {
    const fs::path bk_file = fs::path(config.input.restart_dir) / "fortran_fort21.bin";
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
    worst = std::max(worst, deviation);
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
// nread = 2, non-NAG, lines 739-787), including its bugs: ATAN instead of ATAN2 for the phases,
// s1_real instead of s1_imag in the interpolated Im S, and the signs of Re S and Im S taken
// separately from the linearly interpolated S. Test use only.
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
  const double simag_interpolated = s1_imag + (s2_imag - s1_real) * f;  // the Fortran typo
  const double sreal = std::sqrt(prob_interpolated) * std::cos(phase_interpolated);
  const double simag = std::sqrt(prob_interpolated) * std::sin(phase_interpolated);
  // Fortran SIGN(a, b): |a| with the sign of b (b = 0 counts as positive).
  auto sign = [](double a, double b) { return b >= 0.0 ? std::abs(a) : -std::abs(a); };
  return {sign(sreal, sreal_interpolated), sign(simag, simag_interpolated)};
}

bool have_reference(const std::string& case_name, const std::string& variant) {
  return fs::exists(case_dir(case_name) / ("fine_" + variant) / "stdout.txt");
}

void run_gate(const std::string& case_name, const std::string& variant) {
  if (!have_reference(case_name, variant))
    GTEST_SKIP() << "no Fortran reference; run tools/fortran_reference/run_fine_analysis_reference.sh " << case_name;
  const fs::path dir = case_dir(case_name) / ("fine_" + variant);
  const FineRun run(variant_input(case_name, variant));
  const FineAnalysisResults r = run.run();
  const FineAnalysisGrid& g = r.grid;
  const AnalysisBasis& basis = run.setup.basis;
  const int n_states = basis.n_states();
  const int n_rot = basis.n_rotation();
  const int n_orders = basis.max_order() + 1;

  // 1. Asymptotic info.
  const FortranAsympInfo info = read_fortran_asymp_info(dir / "stdout.txt");
  EXPECT_EQ(g.n_records_available, info.n_records);
  EXPECT_EQ(g.n_records_used, info.n_used);
  EXPECT_EQ(g.size(), info.nen95);
  Comparison c;
  c.check(g.spectral_resolution * 1000.0, info.spres_mev, "spres");
  c.check(g.meaningful_min, info.emin, "emintr");
  c.check(g.meaningful_max, info.emax, "emaxtr");
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
    ASSERT_GE(is, 0);
    for (int i = 0; i < g.size(); ++i)
      u3.check(r.spectrum.energies[i].probability[is + n_states * v], values[i], "unit 3 P");
  }
  u3.report("unit 3 (all states, fine grid)");

  // 5. The interpolated per-energy blocks of the stdout (probabilities: unaffected by the bugs).
  PropagationRun interpolated;
  interpolated.scattering = r.interpolated;
  for (int irot = 0; irot < n_rot; ++irot)
    interpolated.rotation_labels.emplace_back(basis.j(irot * basis.n_diffraction()),
                                              basis.mj(irot * basis.n_diffraction()));
  compare_probabilities_with_fortran(interpolated, read_fortran_probabilities(dir / "stdout.txt"));
  EXPECT_EQ(static_cast<int>(r.skipped.size()), info.n_skip_messages);

  // 6. fort.25: probabilities against the C++; S and phase against the bug emulation.
  const std::vector<Fort25Block> fort25 = read_fortran_fort25(dir / "fort.25");
  ASSERT_EQ(fort25.size(), r.interpolated.energies.size());
  double s_max = 0.0;
  for (const EnergyResult& e : r.spectrum.energies)
    for (const cplx s : e.s_matrix) s_max = std::max(s_max, std::abs(s));
  Comparison p25;
  Comparison s25;
  int n_undecidable = 0;
  int n_other_interval = 0;
  for (std::size_t ie = 0; ie < fort25.size(); ++ie) {
    const EnergyResult& e = r.interpolated.energies[ie];
    const auto [upper, f] = locate_on_grid(g, e.z_energy);
    for (const auto& [key, entry] : fort25[ie].states) {
      const auto& [v, j, mj, n, m] = key;
      const int order = SurfaceLattice(run.config.derived.system).diffraction_order(n, m);
      const int k = run.state(order, n, m, j, mj) + n_states * v;
      p25.check(e.probability[k], entry.p, "fort.25 P");

      const cplx s1 = r.spectrum.energies[upper - 1].s_matrix[k];
      const cplx s2 = r.spectrum.energies[upper].s_matrix[k];
      if (std::abs(s1) < 1e-8 * s_max && std::abs(s2) < 1e-8 * s_max) continue;
      // Signs decided by a linear interpolation within rounding of zero cannot be compared.
      const double re_lin = s1.real() + (s2.real() - s1.real()) * f;
      const double im_lin = s1.imag() + (s2.imag() - s1.real()) * f;
      if (std::abs(re_lin) < 1e-10 * s_max || std::abs(im_lin) < 1e-10 * s_max) {
        n_undecidable = n_undecidable + 1;
        continue;
      }
      cplx emulated = fortran_interpolation_with_bugs(s1, s2, f);
      // On a grid point (f = 0 or 1 to rounding) the Fortran may have chosen the neighbouring
      // interval: its E_z (ensctr, quad precision) differs from the C++ one in the last bits.
      // The emulation of that interval is used if it matches the Fortran value better.
      const cplx printed(parse_fortran_real(entry.re), parse_fortran_real(entry.im));
      const int n_grid = g.size();
      std::optional<cplx> other;
      if (f > 1.0 - 1e-9 && upper + 1 < n_grid)
        other = fortran_interpolation_with_bugs(s2, r.spectrum.energies[upper + 1].s_matrix[k], f - 1.0);
      if (f < 1e-9 && upper >= 2)
        other = fortran_interpolation_with_bugs(r.spectrum.energies[upper - 2].s_matrix[k], s1, f + 1.0);
      if (other && std::abs(*other - printed) < std::abs(emulated - printed)) {
        emulated = *other;
        n_other_interval = n_other_interval + 1;
      }
      s25.check(emulated.real(), entry.re, "fort.25 Re S (emulated)");
      s25.check(emulated.imag(), entry.im, "fort.25 Im S (emulated)");
      const double phase = emulated.real() == 0.0 ? 0.0 : std::atan2(emulated.imag(), emulated.real());
      s25.check(phase, entry.phase, "fort.25 phase (emulated)");
    }
  }
  p25.report("fort.25 probabilities");
  s25.report("fort.25 S-matrix and phases (Fortran interpolation emulated)");
  std::cout << "  fort.25: " << n_undecidable << " values with an undecidable sign skipped, " << n_other_interval
            << " on a grid point matched with the neighbouring interval\n";
}

}  // namespace

TEST(RegressionFineAnalysis, Gate5Full) { run_gate("hardwall_small", "full"); }
TEST(RegressionFineAnalysis, Gate5Truncated) { run_gate("hardwall_small", "truncated"); }
TEST(RegressionFineAnalysis, Gate5Nen95) { run_gate("hardwall_small", "nen95"); }
TEST(RegressionFineAnalysis, Gate5Small) { run_gate("hardwall_small", "small"); }
TEST(RegressionFineAnalysis, Gate5Empty) { run_gate("hardwall_small", "empty"); }
// The corrugated C6v case (60 deg cell): on the BK file of the full Dai-Light Fortran run.
TEST(RegressionFineAnalysis, Gate5DaiLight) { run_gate("cu111_dailight_small", "full"); }
