// Tests: fine analysis (run mode 2) without Fortran output: the time-to-energy transform (direct
// sum, analytic transform, the two-FFT form of Fortran hlftrf), the probaf prefactors against
// the on-the-fly analysis, the interpolation, and the model of the Fortran interpolation bugs
// (docs/fortran_fine_analysis_bugs.md) against the exact S at the analysis energies.
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "app/fine_analysis_run.hpp"
#include "comparison.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "io/fortran_files.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"

namespace fs = std::filesystem;
using namespace spodvr;
using spodvr::testing::nan_aware_max;
using cplx = std::complex<double>;

namespace {

const fs::path case_dir = fs::path(SPODVR_REFERENCE_DIR) / "hardwall_small";
const fs::path fortran_bk = case_dir / "fortran_fort21.bin";

nlohmann::json hardwall_small_input() {
  std::ifstream file(case_dir / "input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["onlyAnalysis"] = true;
  json["RestartDir"] = case_dir.string();
  json["RestartFormat"] = "fortran";
  return json;
}

// Setup of an analysis-only run (gas phase and analysis basis, no 6D potential).
struct RunSetup {
  explicit RunSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        gas(compute_gas_phase(config, grids, angular,
                              compute_asymptotic_potential(*make_potential(config.input), config, grids), fft,
                              LapackEigensolver())),
        analysis(config, angular, gas) {}

  [[nodiscard]] const AnalysisBasis& basis() const { return analysis.basis; }
  [[nodiscard]] double initial_channel_energy() const {
    return basis().channel_energy(basis().initial_state(), config.derived.v_in_index - 1);
  }
  [[nodiscard]] int n_series() const { return basis().n_states() * basis().n_v(); }

  FftwFft fft;
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  GasPhase gas;
  AnalysisSetup analysis;
};

const RunSetup& hardwall_small() {
  static const RunSetup setup(hardwall_small_input());
  return setup;
}

// Random BK records (deterministic seed).
std::vector<BkCoefficients> random_records(const RunSetup& s, int n_records, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  std::vector<BkCoefficients> records;
  for (int j = 0; j < n_records; ++j) {
    BkCoefficients bk(s.basis().n_states(), s.basis().n_v());
    for (cplx& c : bk.values) c = cplx(u(rng), u(rng));
    records.push_back(bk);
  }
  return records;
}

BkSeries make_series(const RunSetup& s, const FineAnalysisGrid& grid, const std::vector<BkCoefficients>& records) {
  BkSeries series(grid, s.basis(), s.initial_channel_energy());
  for (int j = 0; j < static_cast<int>(records.size()); ++j) series.add_record((j + 1) * grid.time_step, records[j]);
  series.check_complete();
  return series;
}

BkSeries fortran_series(const RunSetup& s, const FineAnalysisGrid& grid) {
  BkSeries series(grid, s.basis(), s.initial_channel_energy());
  (void)read_fortran_bk_file(fortran_bk, s.config, s.analysis.rotations, s.basis(),
                             [&](double t, const BkCoefficients& bk) { series.add_record(t, bk); });
  series.check_complete();
  return series;
}

// The transform by its definition (test oracle): A_k(E) = deltat sum_j w_j c_k(j) exp(i E t_j),
// E = E_z (eV) / chaev, at arbitrary energies. Layout k + n_series * i.
std::vector<cplx> direct_amplitudes(const BkSeries& series, const FineAnalysisGrid& grid,
                                    const std::vector<double>& z_energies) {
  const double chaev = units::legacy.energy_factor();
  const int n = grid.n_records_used;
  std::vector<cplx> a(static_cast<std::size_t>(series.n_series()) * z_energies.size());
  for (std::size_t i = 0; i < z_energies.size(); ++i) {
    const double e = z_energies[i] / chaev;
    for (int k = 0; k < series.n_series(); ++k) {
      const std::span<const cplx> c = series.series(k);
      cplx sum(0.0, 0.0);
      for (int j = 0; j <= n; ++j) {
        const double w = (j == 0 || j == n) ? 0.5 : 1.0;
        sum = sum + w * c[j] * std::exp(cplx(0.0, e * j * grid.time_step));
      }
      a[k + static_cast<std::size_t>(series.n_series()) * i] = grid.time_step * sum;
    }
  }
  return a;
}

double max_abs(const std::vector<cplx>& v) {
  double m = 0.0;
  for (const cplx x : v) m = nan_aware_max(m, std::abs(x));
  return m;
}

// Fortran hlftrf literally: two FFTs of length notraf (sign +1) of the Hermitian and the
// anti-Hermitian extension of c; fwr = deltat/2 Re, fwi = deltat/2 Im. Then shbkar.
std::vector<cplx> fortran_hlftrf(const BkSeries& series, const FineAnalysisGrid& grid, const Fft& fft) {
  const int n = grid.n_records_used;
  const int length = grid.transform_length;
  const auto plan = fft.plan_complex(FftLayout{length}, FftSign::Plus, 1.0);
  std::vector<cplx> a(static_cast<std::size_t>(series.n_series()) * grid.size());
  std::vector<cplx> ftraf(length);
  for (int k = 0; k < series.n_series(); ++k) {
    const std::span<const cplx> c = series.series(k);
    for (const double sign : {1.0, -1.0}) {
      ftraf[0] = c[0];
      ftraf[n] = c[n];
      for (int j = 1; j < n; ++j) {
        ftraf[j] = c[j];
        ftraf[length - j] = sign * std::conj(c[j]);
      }
      plan->execute(ftraf);
      for (int i = 0; i < grid.size(); ++i) {
        cplx& out = a[k + static_cast<std::size_t>(series.n_series()) * i];
        const cplx f = ftraf[grid.first_index + i];
        if (sign > 0.0)
          out = cplx(0.5 * grid.time_step * f.real(), out.imag());
        else
          out = cplx(out.real(), 0.5 * grid.time_step * f.imag());
      }
    }
  }
  return a;
}

}  // namespace

// --- grid ------------------------------------------------------------------------------------

TEST(FineAnalysisGrid, HardwallSmall) {
  const RunSetup& s = hardwall_small();
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  EXPECT_EQ(g.n_records_available, 1400);
  EXPECT_EQ(g.n_records_used, 1386);
  EXPECT_EQ(g.transform_length, 2772);
  EXPECT_DOUBLE_EQ(g.effective_time, 13860.0);
  EXPECT_DOUBLE_EQ(g.spectral_resolution, Constants::pi / 13860.0 * units::legacy.energy_factor());
  EXPECT_EQ(g.size(), 65);  // Fortran: nen95 = 65
  // All kept points lie in the confidence interval, their neighbours outside.
  EXPECT_GE(g.z_energy(0), s.config.derived.confidence_e_min);
  EXPECT_LT((g.first_index - 1) * g.spectral_resolution, s.config.derived.confidence_e_min);
  EXPECT_LE(g.z_energy(g.size() - 1), s.config.derived.confidence_e_max);
  EXPECT_GT((g.last_index + 1) * g.spectral_resolution, s.config.derived.confidence_e_max);
}

TEST(FineAnalysisGrid, RecordLimits) {
  const RunSetup& s = hardwall_small();
  EXPECT_THROW((void)make_fine_analysis_grid(s.config, 1399), InputError);  // tEnd needs 1400 records
  EXPECT_EQ(make_fine_analysis_grid(s.config, 2000).n_records_available, 1400);  // only up to tEnd
  nlohmann::json json = hardwall_small_input();
  json["tEnd"] = 10.0;
  json["nStep"] = 5;
  const Config one_record = build_config(parse_input(json));
  EXPECT_THROW((void)make_fine_analysis_grid(one_record, 1400), InputError);  // fewer than 2 records
}

// --- series ----------------------------------------------------------------------------------

TEST(BkSeries, FrequencyShiftAndChecks) {
  const RunSetup& s = hardwall_small();
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  const std::vector<BkCoefficients> records = random_records(s, 1400, 1);
  const BkSeries series = make_series(s, g, records);
  EXPECT_EQ(series.length(), 1387);
  const double e0 = s.initial_channel_energy();
  for (const int k : {0, 7, s.n_series() - 1}) {
    EXPECT_EQ(series.series(k)[0], cplx(0.0, 0.0));
    for (const int j : {1, 500, 1386}) {
      const cplx expected = records[j - 1].values[k] * std::exp(cplx(0.0, e0 * j * g.time_step));
      EXPECT_EQ(series.series(k)[j], expected);
    }
  }

  BkSeries incomplete(g, s.basis(), e0);
  incomplete.add_record(g.time_step, records[0]);
  EXPECT_THROW(incomplete.check_complete(), CheckpointError);
  BkSeries wrong_time(g, s.basis(), e0);
  EXPECT_THROW(wrong_time.add_record(2.0 * g.time_step, records[0]), CheckpointError);
  // Native record times (t_end + (step - n_steps) dt) are j deltat only up to rounding.
  BkSeries rounded(g, s.basis(), e0);
  EXPECT_NO_THROW(rounded.add_record(g.time_step * (1.0 + 1e-14), records[0]));
}

// --- transform -------------------------------------------------------------------------------

TEST(HalfFourierTransform, EqualsDirectSum) {
  const RunSetup& s = hardwall_small();
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  const BkSeries series = make_series(s, g, random_records(s, 1400, 2));
  const std::vector<cplx> fft = half_fourier_transform(series, g, s.fft);
  const std::vector<cplx> direct = direct_amplitudes(series, g, g.z_energies());
  ASSERT_EQ(fft.size(), direct.size());
  const double scale = max_abs(direct);
  for (std::size_t i = 0; i < fft.size(); ++i) ASSERT_LE(std::abs(fft[i] - direct[i]), 1e-12 * scale) << i;
}

TEST(HalfFourierTransform, EqualsFortranTwoFftForm) {
  const RunSetup& s = hardwall_small();
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  const BkSeries series = fortran_series(s, g);
  const std::vector<cplx> one = half_fourier_transform(series, g, s.fft);
  const std::vector<cplx> two = fortran_hlftrf(series, g, s.fft);
  const double scale = max_abs(two);
  for (std::size_t i = 0; i < one.size(); ++i) ASSERT_LE(std::abs(one[i] - two[i]), 1e-13 * scale) << i;
}

// A Gaussian pulse well inside (0, T): the trapezoid sum equals the continuous transform
// sigma sqrt(2 pi) exp(-sigma^2 w^2 / 2) exp(i w t0), w = E + ev0jl - E0, to ~1e-12.
TEST(HalfFourierTransform, GaussianPulseAnalytic) {
  const RunSetup& s = hardwall_small();
  const double chaev = units::legacy.energy_factor();
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  const double t0 = 7000.0;
  const double sigma = 800.0;
  const double e0 = 0.4 / chaev + s.initial_channel_energy();  // peak at E_z = 0.4 eV
  std::vector<BkCoefficients> records;
  for (int j = 1; j <= 1400; ++j) {
    const double t = j * g.time_step;
    BkCoefficients bk(s.basis().n_states(), s.basis().n_v());
    for (cplx& c : bk.values) c = std::exp(-(t - t0) * (t - t0) / (2.0 * sigma * sigma)) * std::exp(cplx(0.0, -e0 * t));
    records.push_back(bk);
  }
  const BkSeries series = make_series(s, g, records);
  const std::vector<cplx> a = half_fourier_transform(series, g, s.fft);
  const double peak = sigma * std::sqrt(2.0 * Constants::pi);
  for (int i = 0; i < g.size(); ++i) {
    const double w = g.z_energy(i) / chaev + s.initial_channel_energy() - e0;
    const cplx exact = peak * std::exp(-sigma * sigma * w * w / 2.0) * std::exp(cplx(0.0, w * t0));
    EXPECT_LE(std::abs(a[static_cast<std::size_t>(s.n_series()) * i] - exact), 1e-12 * peak) << "i = " << i;
  }
}

// --- prefactors: the probaf formulas against those of the on-the-fly analysis ----------------

namespace {

// The same trapezoid sum evaluated by the on-the-fly accumulator (the last record weighted 1/2)
// and by fine_spectrum at the NrOfEs analysis energies.
void compare_with_on_the_fly(const RunSetup& s) {
  const double chaev = units::legacy.energy_factor();
  const int n_records = static_cast<int>(std::lround(s.config.input.t_end / s.config.derived.analysis_time_step));
  const FineAnalysisGrid g = make_fine_analysis_grid(s.config, n_records);
  std::vector<BkCoefficients> records = random_records(s, g.n_records_used, 3);

  SMatrixAccumulator on_the_fly(s.config, s.basis(), s.gas);
  for (int j = 1; j <= g.n_records_used; ++j) {
    BkCoefficients bk = records[j - 1];
    if (j == g.n_records_used)
      for (cplx& c : bk.values) c = 0.5 * c;
    on_the_fly.add(bk, j * g.time_step);
  }
  const ScatteringResults expected = on_the_fly.results();

  std::vector<double> z_energies;
  for (const double e : on_the_fly.energies().collision)
    z_energies.push_back(e - s.basis().diffraction_energy(0) * chaev);
  const BkSeries series = make_series(s, g, records);
  const ScatteringResults fine =
      fine_spectrum(s.config, z_energies, s.basis(), s.gas, direct_amplitudes(series, g, z_energies));

  ASSERT_EQ(fine.energies.size(), expected.energies.size());
  int n_compared = 0;
  for (std::size_t i = 0; i < fine.energies.size(); ++i) {
    const EnergyResult& a = fine.energies[i];
    const EnergyResult& b = expected.energies[i];
    EXPECT_NEAR(a.total_energy, b.total_energy, 1e-12);
    EXPECT_NEAR(a.collision_energy, b.collision_energy, 1e-12);
    double scale = 0.0;
    for (const cplx x : b.s_matrix) scale = nan_aware_max(scale, std::abs(x));
    for (std::size_t k = 0; k < a.s_matrix.size(); ++k) {
      if (!a.open[k] || !b.open[k]) continue;
      EXPECT_LE(std::abs(a.s_matrix[k] - b.s_matrix[k]), 1e-10 * scale) << "energy " << i << " state " << k;
      n_compared = n_compared + 1;
    }
  }
  EXPECT_GT(n_compared, 0);
}

}  // namespace

TEST(FineSpectrum, PrefactorsMatchOnTheFlyAnalysis) { compare_with_on_the_fly(hardwall_small()); }

// With parallel momentum (EParallel != 0) the two analyses use different expressions for the
// momentum of the elastic channel (sk0 from E_coll minus k_parallel^2, or from E_z).
TEST(FineSpectrum, PrefactorsMatchOnTheFlyAnalysisWithParallelMomentum) {
  nlohmann::json json = hardwall_small_input();
  json["EParallel"] = 0.05;
  json["PhiParallel"] = 30.0;
  const RunSetup s(json);
  ASSERT_GT(s.basis().diffraction_energy(0), 0.0);
  compare_with_on_the_fly(s);
}

// --- interpolation ---------------------------------------------------------------------------

TEST(InterpolateSMatrix, GridPointsAndLinearProbability) {
  const cplx s1 = std::polar(0.8, 3.0);
  const cplx s2 = std::polar(0.6, -3.0);  // the phase passes through pi
  EXPECT_NEAR(std::abs(interpolate_s_matrix(s1, s2, 0.0) - s1), 0.0, 1e-15);
  EXPECT_NEAR(std::abs(interpolate_s_matrix(s1, s2, 1.0) - s2), 0.0, 1e-15);
  for (const double f : {0.25, 0.5, 0.75}) {
    const cplx s = interpolate_s_matrix(s1, s2, f);
    EXPECT_NEAR(std::norm(s), 0.64 + (0.36 - 0.64) * f, 1e-15);       // |S|^2 linear (as Fortran)
    EXPECT_GT(std::abs(std::arg(s)), 3.0 - 1e-12);                     // across pi, not through 0
    EXPECT_NEAR(std::arg(s), std::arg(s1 + (s2 - s1) * f), 1e-15);     // phase of the chord
  }
}

TEST(InterpolateSMatrix, SignChangesAndZeros) {
  // Re S changes sign (the case the Fortran ATAN gets wrong): the phase moves smoothly.
  const cplx s1 = std::polar(1.0, Constants::pi / 2.0 - 0.1);
  const cplx s2 = std::polar(1.0, Constants::pi / 2.0 + 0.1);
  EXPECT_NEAR(std::arg(interpolate_s_matrix(s1, s2, 0.5)), Constants::pi / 2.0, 1e-12);
  // A zero end point (closed channel): the phase of the other one.
  EXPECT_NEAR(std::arg(interpolate_s_matrix(cplx(0.0, 0.0), s2, 0.3)), std::arg(s2), 1e-15);
  EXPECT_NEAR(std::norm(interpolate_s_matrix(cplx(0.0, 0.0), s2, 0.3)), 0.3, 1e-15);
  EXPECT_EQ(interpolate_s_matrix(cplx(0.0, 0.0), cplx(0.0, 0.0), 0.5), cplx(0.0, 0.0));
}

TEST(LocateOnGrid, NeighboursAndClamps) {
  FineAnalysisGrid g{};
  g.first_index = 10;
  g.last_index = 14;
  g.spectral_resolution = 0.01;  // points 0.10 .. 0.14
  EXPECT_EQ(locate_on_grid(g, 0.115).upper, 2);
  EXPECT_NEAR(locate_on_grid(g, 0.115).f, 0.5, 1e-12);
  EXPECT_EQ(locate_on_grid(g, 0.12).upper, 2);  // on a point: first point with E >= E_z, f = 1
  EXPECT_NEAR(locate_on_grid(g, 0.12).f, 1.0, 1e-12);
  EXPECT_EQ(locate_on_grid(g, 0.10).upper, 1);  // clamped: f = 0
  EXPECT_NEAR(locate_on_grid(g, 0.10).f, 0.0, 1e-12);
  EXPECT_EQ(locate_on_grid(g, 0.14 + 1e-13).upper, 4);  // above the last point: clamped, f > 1
  EXPECT_GT(locate_on_grid(g, 0.14 + 1e-13).f, 1.0);
}

TEST(FineAnalysis, SkippedEnergiesAndSmallGrids) {
  const RunSetup& s = hardwall_small();
  FineAnalysisGrid g = make_fine_analysis_grid(s.config, 1400);
  const BkSeries series = make_series(s, g, random_records(s, 1400, 4));
  const ScatteringResults spectrum =
      fine_spectrum(s.config, g.z_energies(), s.basis(), s.gas, half_fourier_transform(series, g, s.fft));
  const AnalysisEnergies energies = make_analysis_energies(s.config, s.gas);
  const FineAnalysisResults r = interpolate_to_analysis_energies(s.config, g, s.basis(), s.gas, energies, spectrum);
  // hardwall_small: energies 0.2 .. 0.6 eV; the grid covers 0.2035 .. 0.5983 eV.
  EXPECT_EQ(r.interpolated.energies.size() + r.skipped.size(), static_cast<std::size_t>(energies.size()));
  EXPECT_EQ(r.skipped.size(), 2U);  // as the Fortran (two messages in fine_full/stdout.txt)

  // One grid point: nothing can be interpolated (Fortran reads out of bounds).
  g.last_index = g.first_index;
  ScatteringResults one;
  one.energies.push_back(spectrum.energies[0]);
  const FineAnalysisResults r1 = interpolate_to_analysis_energies(s.config, g, s.basis(), s.gas, energies, one);
  EXPECT_TRUE(r1.interpolated.energies.empty());
  EXPECT_EQ(r1.skipped.size(), static_cast<std::size_t>(energies.size()));
}

// --- the Fortran interpolation bugs ----------------------------------------------------------

namespace {

// Fortran 5.19 smatrix_calc (nread = 2, non-NAG) literally, bugs included: the documented 5.19
// model (docs/fortran_fine_analysis_bugs.md). The gate 5 test emulates the reference 5.24 r1713,
// which fixed the s1_real typo.
cplx fortran519_interpolation_with_bugs(cplx s1, cplx s2, double f) {
  const double prob1 = s1.real() * s1.real() + s1.imag() * s1.imag();
  const double prob2 = s2.real() * s2.real() + s2.imag() * s2.imag();
  const double prob = prob1 + (prob2 - prob1) * f;
  const double phase1 = fortran_compare::rne(s1.real(), 0.0) ? std::atan(s1.imag() / s1.real()) : 0.0;
  const double phase2 = fortran_compare::rne(s2.real(), 0.0) ? std::atan(s2.imag() / s2.real()) : 0.0;
  const double phase = phase1 + (phase2 - phase1) * f;
  const double re_lin = s1.real() + (s2.real() - s1.real()) * f;
  const double im_lin = s1.imag() + (s2.imag() - s1.real()) * f;  // typo: s1_real
  auto sign = [](double a, double b) { return b >= 0.0 ? std::abs(a) : -std::abs(a); };
  return {sign(std::sqrt(prob) * std::cos(phase), re_lin), sign(std::sqrt(prob) * std::sin(phase), im_lin)};
}

// What the Fortran code intends: |S|^2 and the (unwrapped) phase interpolated linearly.
cplx intended_interpolation(cplx s1, cplx s2, double f) {
  const double prob = std::norm(s1) + (std::norm(s2) - std::norm(s1)) * f;
  double d = std::arg(s2) - std::arg(s1);
  d = d - 2.0 * Constants::pi * std::round(d / (2.0 * Constants::pi));
  return std::polar(std::sqrt(prob), std::arg(s1) + d * f);
}

}  // namespace

// On the BK coefficients of the full Fortran run (hardwall_small): (1) the Fortran result differs
// from the intended interpolation only where predicted (Re S or Im S changes sign between the
// grid points, or the typo changes the sign of the interpolated Im S); (2) the C++
// interpolation is close to the exact S (the transform evaluated at the analysis energy),
// the Fortran one is not. The numbers are quoted in docs/fortran_fine_analysis_bugs.md.
TEST(FortranInterpolationBugs, OnlyPredictedPointsDifferAndFixIsAccurate) {
  const RunSetup& s = hardwall_small();
  const double chaev = units::legacy.energy_factor();
  const FineAnalysisResults r = run_fine_analysis(s.config, s.analysis, s.gas, BkInput{fortran_bk, true}, s.fft);
  const FineAnalysisGrid& g = r.grid;
  const BkSeries series = fortran_series(s, g);

  std::vector<double> z_energies;
  for (const EnergyResult& e : r.interpolated.energies) z_energies.push_back(e.z_energy);
  const ScatteringResults exact =
      fine_spectrum(s.config, z_energies, s.basis(), s.gas, direct_amplitudes(series, g, z_energies));

  double s_max = 0.0;
  for (const EnergyResult& e : r.spectrum.energies)
    for (const cplx x : e.s_matrix) s_max = nan_aware_max(s_max, std::abs(x));
  const double populated = 1e-3 * s_max;

  int n_pairs = 0, n_differ = 0, n_atan = 0, n_sign = 0, n_typo = 0;
  double err_fixed = 0.0, err_fortran = 0.0, err_fixed_p = 0.0, err_fortran_phase = 0.0;
  for (std::size_t ie = 0; ie < r.interpolated.energies.size(); ++ie) {
    const auto [upper, f] = locate_on_grid(g, z_energies[ie]);
    for (std::size_t k = 0; k < r.interpolated.energies[ie].s_matrix.size(); ++k) {
      const cplx s1 = r.spectrum.energies[upper - 1].s_matrix[k];
      const cplx s2 = r.spectrum.energies[upper].s_matrix[k];
      const cplx s_exact = exact.energies[ie].s_matrix[k];
      if (std::abs(s1) < populated || std::abs(s2) < populated) continue;
      n_pairs = n_pairs + 1;

      const cplx fixed = r.interpolated.energies[ie].s_matrix[k];
      const cplx fortran = fortran519_interpolation_with_bugs(s1, s2, f);
      const cplx intended = intended_interpolation(s1, s2, f);
      err_fixed = nan_aware_max(err_fixed, std::abs(fixed - s_exact) / std::abs(s_exact));
      err_fortran = nan_aware_max(err_fortran, std::abs(fortran - s_exact) / std::abs(s_exact));
      err_fixed_p = nan_aware_max(err_fixed_p, std::abs(std::norm(fixed) - std::norm(s_exact)) / std::norm(s_exact));
      err_fortran_phase = nan_aware_max(err_fortran_phase, std::abs(std::arg(fortran / s_exact)));

      const bool atan_case = s1.real() * s2.real() < 0.0;
      const bool sign_case = s1.imag() * s2.imag() < 0.0;
      const bool typo_case =
          (s1.imag() + (s2.imag() - s1.real()) * f >= 0.0) != (s1.imag() + (s2.imag() - s1.imag()) * f >= 0.0);
      if (std::abs(fortran - intended) > 1e-12 * s_max) {
        n_differ = n_differ + 1;
        EXPECT_TRUE(atan_case || sign_case || typo_case)
            << "unpredicted difference at E_z = " << z_energies[ie] << ", state " << k;
      }
      if (atan_case) n_atan = n_atan + 1;
      if (sign_case) n_sign = n_sign + 1;
      if (typo_case) n_typo = n_typo + 1;
    }
  }
  std::cout << "  populated (energy, state) pairs: " << n_pairs << "; Fortran differs from the intended result at "
            << n_differ << " (Re S sign change: " << n_atan << ", Im S sign change: " << n_sign
            << ", typo flips the sign: " << n_typo << ")\n"
            << "  largest relative error against the exact S: C++ " << err_fixed << " (|S|^2: " << err_fixed_p
            << "), Fortran " << err_fortran << " (phase " << err_fortran_phase << " rad)\n"
            << "  spectrum: " << r.large_phase_steps << " populated phase steps > pi/2; "
            << g.spectral_resolution * 1000.0 << " meV per grid step; E_z range " << z_energies.front() << " .. "
            << z_energies.back() << " eV (" << chaev << ")\n";
  EXPECT_GT(n_pairs, 0);
  EXPECT_LT(err_fixed, 0.02);
  EXPECT_LT(err_fixed, err_fortran);
}

// A synthetic S(E) = a(E) exp(i phi(E)) with a fast phase (0.9 rad per grid step), so that all
// three Fortran error classes occur. The Fortran result has O(1) errors exactly in these
// classes; the C++ interpolation stays close to S(E).
TEST(FortranInterpolationBugs, SyntheticFastPhase) {
  auto exact = [](double x) {  // x: energy in grid steps
    return std::polar(0.5 + 0.3 * std::sin(0.2 * x), 0.9 * x + 0.3);
  };
  int n = 0, n_atan = 0, n_sign = 0, n_typo = 0, n_wrong = 0;
  double err_fixed = 0.0;
  for (int i = 0; i < 40; ++i)
    for (const double f : {0.1, 0.3, 0.5, 0.7, 0.9}) {
      const cplx s1 = exact(i);
      const cplx s2 = exact(i + 1);
      const cplx s = exact(i + f);
      const cplx fixed = interpolate_s_matrix(s1, s2, f);
      const cplx fortran = fortran519_interpolation_with_bugs(s1, s2, f);
      const cplx intended = intended_interpolation(s1, s2, f);
      err_fixed = nan_aware_max(err_fixed, std::abs(fixed - s) / std::abs(s));
      const bool atan_case = s1.real() * s2.real() < 0.0;
      const bool sign_case = s1.imag() * s2.imag() < 0.0;
      const bool typo_case =
          (s1.imag() + (s2.imag() - s1.real()) * f >= 0.0) != (s1.imag() + (s2.imag() - s1.imag()) * f >= 0.0);
      const bool wrong = std::abs(fortran - intended) > 1e-12;
      if (wrong) {
        EXPECT_TRUE(atan_case || sign_case || typo_case) << "i = " << i << ", f = " << f;
        n_wrong = n_wrong + 1;
        n_atan = n_atan + (atan_case ? 1 : 0);
        n_sign = n_sign + (sign_case ? 1 : 0);
        n_typo = n_typo + (typo_case ? 1 : 0);
      }
      n = n + 1;
    }
  std::cout << "  synthetic: " << n << " interpolations, Fortran wrong at " << n_wrong << " (with Re S sign change "
            << n_atan << ", Im S sign change " << n_sign << ", typo sign flip " << n_typo
            << "); C++ largest relative error " << err_fixed << "\n";
  EXPECT_GT(n_atan, 0);
  EXPECT_GT(n_sign, 0);
  EXPECT_GT(n_typo, 0);
  EXPECT_LT(err_fixed, 0.12);  // second order in the phase step (0.9 rad)
}
