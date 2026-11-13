// Tests: the flux analysis in the time loop (WavePacketPropagation): its coefficients are the sums of the
// flux surfaces of every analysis step, whatever the slot and thread count, a restart continues them
// bitwise, and the rest of the run is the same with the flux analysis on; its result files and table
// (the writer and report of run.cpp) do not depend on the thread count.
#include <gtest/gtest.h>
#include <omp.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "backend_setup.hpp"
#include "core/array_index.hpp"
#include "flux_propagation_run.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "propagation_reference.hpp"
#include "report/flux_files.hpp"
#include "report/flux_report.hpp"
#include "symmetry/grid_symmetry.hpp"
#include "symmetry/point_group.hpp"
#include "symmetry/surface_lattice.hpp"

using namespace spodvr;
using namespace spodvr::testing;
using cplx = std::complex<double>;

namespace {

// 12 analysis steps (0 .. 120 a.u.); averages at 0 and 100 a.u.
constexpr int n_analysis_steps = 12;
constexpr double analysis_step = 10.0;  // a.u.

// A run with `threads` OpenMP threads (the caller's count restored afterwards).
FluxPropagationRun run_with_threads(const BackendSetup& s, int slots, int threads) {
  const int caller_threads = omp_get_max_threads();
  omp_set_num_threads(threads);
  FluxPropagationRun run = run_flux_propagation(s, slots);
  omp_set_num_threads(caller_threads);
  return run;
}

const BackendSetup& flux_setup() {
  static const BackendSetup s(flux_run_case, flux_run_input(0.0, 120.0, 60));
  return s;
}

// The reference run of these tests: the budget's slots, 4 threads.
const FluxPropagationRun& budget_run() {
  static const FluxPropagationRun run = run_with_threads(flux_setup(), flux_budget_slots(flux_setup()), 4);
  return run;
}

}  // namespace

// The coefficients after run() are sum_t c_n(t) psi(t) over all 12 analysis steps, recomputed here from the
// flux surfaces the run handed to its observer (each at the time at the end of its step, after the BK
// analysis of that step), bitwise. 8 slots for 12 steps: one flush when full, 4 steps flushed at the end
// of run(), so nothing is pending after it.
TEST(FluxPropagation, CoefficientsAreTheSumsOfEveryAnalysisStep) {
  const FluxPropagationRun& run = budget_run();
  ASSERT_EQ(run.flux_times.size(), std::size_t(n_analysis_steps));
  EXPECT_EQ(run.flux_times, run.bk_times);
  for (int k = 0; k < n_analysis_steps; ++k) EXPECT_EQ(run.flux_times[k], analysis_step * (k + 1)) << "step " << k;
  EXPECT_EQ(run.slot_count, max_flux_slots);
  ASSERT_NE(n_analysis_steps % run.slot_count, 0);  // steps are left for the flush at the end of run()
  EXPECT_EQ(run.pending_after_run, 0);

  ASSERT_TRUE(run.coefficients && run.recomputed);
  EXPECT_EQ(bit_differences(*run.coefficients, *run.recomputed), 0);
  // Values of order 1 (from the test wave packet), not rounding noise.
  const auto large = std::ranges::count_if(run.coefficients->phi, [](cplx c) { return std::abs(c) > 1e-3; });
  EXPECT_GT(large, std::ssize(run.coefficients->phi) / 2);
}

// The slots only set how many steps one flush adds: 1 slot (a flush per step) and 5 (two flushes when full,
// two steps at the end of run()) give bitwise the coefficients of the budget's 8.
TEST(FluxPropagation, CoefficientsDoNotDependOnTheSlotCount) {
  for (const int slots : {1, 5}) {
    SCOPED_TRACE(::testing::Message() << slots << " slots");
    const FluxPropagationRun run = run_with_threads(flux_setup(), slots, 4);
    EXPECT_EQ(run.slot_count, slots);
    ASSERT_TRUE(run.coefficients);
    EXPECT_EQ(bit_differences(*run.coefficients, *budget_run().coefficients), 0);
  }
}

// Propagation, flux surface and flush: 1 and 4 threads give bitwise the same coefficients.
TEST(FluxPropagation, CoefficientsDoNotDependOnTheThreadCount) {
  const FluxPropagationRun one_thread = run_with_threads(flux_setup(), flux_budget_slots(flux_setup()), 1);
  ASSERT_TRUE(one_thread.coefficients);
  EXPECT_EQ(bit_differences(*one_thread.coefficients, *budget_run().coefficients), 0);
}

// The flux analysis only reads the wave packet (as in the Fortran, where fort.21 is bitwise the same with
// it): the BK coefficients and the averages of the run without it are bitwise those of the run with it.
TEST(FluxPropagation, TheFluxAnalysisIsPassive) {
  const BackendSetup without(flux_run_case, flux_run_input(0.0, 120.0, 60, false));
  const FluxPropagationRun run = run_with_threads(without, 1, 4);
  EXPECT_FALSE(run.coefficients);
  EXPECT_TRUE(run.flux_times.empty());
  const FluxPropagationRun& with = budget_run();
  ASSERT_EQ(run.bk.size(), with.bk.size());
  for (std::size_t k = 0; k < run.bk.size(); ++k) EXPECT_EQ(bit_differences(run.bk[k], with.bk[k]), 0) << "step " << k;
  ASSERT_EQ(run.averages.size(), with.averages.size());
  for (std::size_t k = 0; k < run.averages.size(); ++k) {
    EXPECT_EQ(run.averages[k].norm, with.averages[k].norm);
    EXPECT_EQ(run.averages[k].z_mean, with.averages[k].z_mean);
    EXPECT_EQ(run.averages[k].r_mean, with.averages[k].r_mean);
    EXPECT_EQ(run.averages[k].rest_norm, with.averages[k].rest_norm);
  }
}

// flux Q4: c_n weighs a step with its length nanspo dt, not with the input dtAna. dtAna = 10.5 a.u. at
// dtStep = 2 gives nanspo = INT(5.25) = 5, steps of 10 a.u. at the same times as with dtAna = 10 (tEnd
// 126, a multiple of dtAna; the last 3 time steps are not propagated, as in the Fortran): bitwise the
// coefficients of the run with dtAna = 10. With the Fortran's weight dtAna they would be 1.05 times those.
// This is the non-multiple case of flux Q4; the truncating and rounding cases are listed there.
TEST(FluxPropagation, FluxQ4WeightIsTheEffectiveAnalysisStep) {
  nlohmann::json changes = flux_run_input(0.0, 126.0, 63);
  changes["dtAna"] = 10.5;
  changes["dtAve"] = 105.0;
  const BackendSetup s(flux_run_case, changes);
  ASSERT_EQ(s.config.derived.steps_per_analysis, 5);
  ASSERT_EQ(s.config.derived.analysis_time_step, analysis_step);
  const FluxPropagationRun run = run_with_threads(s, flux_budget_slots(s), 4);
  ASSERT_EQ(run.flux_times, budget_run().flux_times);
  ASSERT_TRUE(run.coefficients && run.recomputed);
  EXPECT_EQ(bit_differences(*run.coefficients, *run.recomputed), 0);
  EXPECT_EQ(bit_differences(*run.coefficients, *budget_run().coefficients), 0);
}

// Restart: a run split at 60 a.u. continues from the first part's wave packet, analysis surface and flux
// coefficients (read_flux_coefficients, before run()) to bitwise the coefficients of the unsplit run; the
// reaction probability from them is the one of the unsplit coefficients.
TEST(FluxPropagation, SplitRunContinuesTheCoefficientsBitwise) {
  const BackendSetup first_setup(flux_run_case, flux_run_input(0.0, 60.0, 30));
  WavePacketPropagation first(first_setup.config, first_setup.grids, first_setup.angular, first_setup.potential,
                              first_setup.gas, first_setup.fft, make_cpu_backend, 3);
  start_from_test_wave_packet(first);
  first.run(PropagationObserver{});

  nlohmann::json restart = flux_run_input(60.0, 120.0, 30);
  restart["NewCalc"] = false;
  restart["RestartDir"] = "unused";  // the restart data come from `first`, in memory
  const BackendSetup second_setup(flux_run_case, restart);
  WavePacketPropagation second(second_setup.config, second_setup.grids, second_setup.angular, second_setup.potential,
                               second_setup.gas, second_setup.fft, make_cpu_backend, 3);
  second.read_wave_packet([&](WavePacket& psi) {
    std::ranges::copy(first.wave_packet_on_host().scattering.values, psi.scattering.values.begin());
    std::ranges::copy(first.wave_packet_on_host().specular.values, psi.specular.values.begin());
  });
  second.restore_analysis_surface(first.analysis_surface().state());
  second.read_flux_coefficients([&](FluxCoefficients& c) {
    c.phi = first.flux()->coefficients().phi;
    c.phi_derivative = first.flux()->coefficients().phi_derivative;
  });
  second.run(PropagationObserver{});

  const FluxCoefficients& unsplit = *budget_run().coefficients;
  EXPECT_EQ(bit_differences(second.flux()->coefficients(), unsplit), 0);
  const FluxResults results = second.flux_results();
  const FluxResults expected =
      flux_results(unsplit, second.flux()->energies(),
                   make_flux_grid(second_setup.config, second_setup.grids, second_setup.basis.size()));
  EXPECT_EQ(results.probability, expected.probability);
  EXPECT_EQ(results.density, expected.density);
}

// The comparison of these tests fails on non-finite values: the same NaN or infinity on both sides
// (as a backend writing NaN into a slot would leave in the run and in its recomputation alike) counts
// as a difference, and so does -0 against +0; equal finite values do not.
TEST(FluxTestHelpers, BitDifferencesCountNonFiniteValues) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<cplx> finite = {{1.5, -2.0}, {0.0, 3.25}};
  EXPECT_EQ(bit_differences(finite, finite), 0);
  const std::vector<cplx> with_nan = {{1.5, nan}, {0.0, 3.25}};
  EXPECT_EQ(bit_differences(with_nan, with_nan), 1);
  const std::vector<cplx> with_inf = {{inf, -2.0}, {-inf, 3.25}};
  EXPECT_EQ(bit_differences(with_inf, with_inf), 2);
  const std::vector<cplx> negative_zero = {{1.5, -2.0}, {-0.0, 3.25}};
  EXPECT_EQ(bit_differences(finite, negative_zero), 1);
  EXPECT_EQ(bit_differences(finite, std::vector<cplx>(1, cplx(1.5, -2.0))), 4);  // other sizes
}

// Without doFluxAnalysis there is no flux accumulator, and asking for its results or reading
// coefficients into it is a programming error.
TEST(FluxPropagation, WithoutTheFluxAnalysisThereIsNone) {
  const BackendSetup s(flux_run_case, flux_run_input(0.0, 120.0, 60, false));
  WavePacketPropagation propagation(s.config, s.grids, s.angular, s.potential, s.gas, s.fft, make_cpu_backend, 1);
  EXPECT_FALSE(propagation.flux().has_value());
  EXPECT_THROW((void)propagation.flux_results(), std::logic_error);
  EXPECT_THROW(propagation.read_flux_coefficients([](FluxCoefficients&) {}), std::logic_error);
}

namespace {

std::string file_contents(const std::filesystem::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream text;
  text << file.rdbuf();
  return text.str();
}

// What run.cpp does with the flux results of a propagation, for the input of these tests from the test
// wave packet with `threads` threads: run_cpp_propagation, then the result table and the files of
// run.cpp (print_flux_results, write_flux_files) into dir. Returns the results and the printed table.
std::pair<FluxResults, std::string> flux_output(int threads, const std::filesystem::path& dir) {
  nlohmann::json input = reference_case_input(flux_run_case);
  input.update(flux_run_input(0.0, 120.0, 60));
  const int caller_threads = omp_get_max_threads();
  omp_set_num_threads(threads);
  const PropagationRun run = run_cpp_propagation(input, make_cpu_backend, start_from_test_wave_packet);
  omp_set_num_threads(caller_threads);
  EXPECT_TRUE(run.flux.has_value());
  if (!run.flux) return {};
  const Config config = build_config(parse_input(input));
  std::filesystem::create_directories(dir);
  (void)write_flux_files(dir, config, make_spatial_grids(config), *run.flux);
  std::ostringstream table;
  print_flux_results(table, Style{false}, *run.flux);
  return {*run.flux, table.str()};
}

}  // namespace

// End to end: the propagation's flux results through the writer and report of run.cpp give the three
// files and the table, byte-identical with 1 and 4 threads; the results are those of the coefficients of
// the run above (budget_run) and are finite and nonzero (from the test wave packet, not rounding noise).
TEST(FluxOutput, FilesAndTableOfOneAndFourThreadsAreByteIdentical) {
  const std::filesystem::path root =
      std::filesystem::temp_directory_path() / ("spodvr_flux_output_" + std::to_string(::getpid()));
  std::filesystem::remove_all(root);
  const auto [one, one_table] = flux_output(1, root / "threads_1");
  const auto [four, four_table] = flux_output(4, root / "threads_4");
  for (const char* name : {"flux_reaction.tsv", "flux_density.tsv", "flux_analysis.json"}) {
    ASSERT_TRUE(std::filesystem::exists(root / "threads_1" / name)) << name;
    const std::string text = file_contents(root / "threads_1" / name);
    EXPECT_FALSE(text.empty()) << name;
    EXPECT_EQ(text, file_contents(root / "threads_4" / name)) << name;
  }
  EXPECT_EQ(one_table, four_table);
  EXPECT_EQ(one.probability, four.probability);
  std::filesystem::remove_all(root);

  const BackendSetup& s = flux_setup();
  const FluxResults expected = flux_results(
      *budget_run().coefficients, make_flux_energies(s.config, s.tables.parallel_energy, s.gas.initial.energy),
      make_flux_grid(s.config, s.grids, s.basis.size()));
  EXPECT_EQ(one.probability, expected.probability);
  EXPECT_EQ(one.density, expected.density);
  ASSERT_EQ(one.probability.size(), 3U);
  for (const double p : one.probability) EXPECT_TRUE(std::isfinite(p) && p != 0.0) << p;
}

namespace {

// The (x, y) point of the grid index ixy = iy + ny ix (GridSymmetry::apply_to_xy).
FluxGridXY grid_point(int ixy, int ny) { return {ixy / ny, ixy % ny}; }

// The largest |P(x', y'; E) - P(x, y; E)| / max |P| over every (x, y) and energy, where (x', y') is the
// image of (x, y) under operation `op` of `symmetry`, with ix shifted by `ix_shift` (0: none).
double density_asymmetry(const FluxResults& r, const GridSymmetry& symmetry, int op, int ix_shift) {
  const int nx = symmetry.nx();
  const int ny = symmetry.ny();
  double worst = 0.0;
  double largest = 0.0;
  for (std::size_t n = 0; n < r.probability.size(); ++n)
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy) {
        const FluxGridXY image = grid_point(symmetry.apply_to_xy(op, ix, iy), ny);
        const double p = r.density[r.density_index(ix, iy, n)];
        const double q = r.density[r.density_index((image.ix + ix_shift) % nx, image.iy, n)];
        worst = nan_aware_max(worst, std::abs(q - p));
        largest = nan_aware_max(largest, std::abs(p));
      }
  return worst / largest;
}

}  // namespace

// C6v invariance of P(x, y; E) of the gate-6 case (cu111_dailight_short: j = 0, normal incidence, a C6v
// PES). With nphi = 10 the usable grid operations are C2v: the identity, C2 and two mirrors. The
// discrete kinetic operator of the 60-degree grid is not invariant under them: the cross term
// -2 cos(skew) kx ky / 2M differs in sign between a Nyquist momentum's grid image and its true image
// (C2), and a mirror maps momenta out of the FFT box (aliasing). So P(x, y; E) of the run, as the
// Fortran's, is symmetric only to 1.7e-4 (C2) and 2.1e-4 (mirrors) of its largest value (both codes,
// measured 2026-11-13). Dropping that cross term on the Nyquist row and column makes the propagation
// C2-symmetric; then P(x, y; E) is too, to rounding: the flux surface, c_n and flux_results keep the
// wave function's symmetry. Control: the same check with ix shifted by one (a density layout
// off by one) fails. A full run of the case: Long suite.
TEST(FluxPropagation, DensityIsC2InvariantWithAC2SymmetricKineticOperatorLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  BackendSetup s("cu111_dailight_short");
  const InputParameters& p = s.config.input;
  const DerivedParameters& d = s.config.derived;
  ASSERT_EQ(p.e_parallel, 0.0);  // k0 = 0: momenta (i or i - n) dk
  auto momentum = [](int i, int n, double dk) { return (i <= n / 2 ? i : i - n) * dk; };
  const Index2D xy{std::size_t(p.ny), std::size_t(p.nx)};  // ixy = iy + ny ix
  for (int ix = 0; ix < p.nx; ++ix)
    for (int iy = 0; iy < p.ny; ++iy) {
      if (ix != p.nx / 2 && iy != p.ny / 2) continue;
      const double kx = momentum(ix, p.nx, d.delta_kx);
      const double ky = momentum(iy, p.ny, d.delta_ky);
      const double energy = (0.5 / d.total_mass) * (kx * kx) + (0.5 / d.total_mass) * (ky * ky);
      s.tables.xy_kinetic.full[xy(iy, ix)] = exp_minus_i(energy * p.time_step);
      s.tables.xy_kinetic.half[xy(iy, ix)] = exp_minus_i(energy * p.time_step / 2.0);
    }

  const std::unique_ptr<PropagationBackend> backend = s.cpu();
  const FluxEnergies energies = make_flux_energies(s.config, s.tables.parallel_energy, s.gas.initial.energy);
  const int n_blocks = s.basis.size() * p.nx * p.ny;
  FluxCoefficients sums(n_blocks, s.grids.nz_potential, energies.size());
  FluxSurfaceWave surface(n_blocks, s.grids.nz_potential);
  const int n_steps = static_cast<int>(std::lround((p.t_end - p.t_start) / d.analysis_time_step));
  for (int step = 1; step <= n_steps; ++step) {
    backend->propagate_analysis_step();
    backend->flux_surface(surface);
    add_flux_surface_by_hand(sums, energies, p.t_start + step * d.analysis_time_step, d.analysis_time_step, surface);
  }
  const FluxResults r = flux_results(sums, energies, make_flux_grid(s.config, s.grids, s.basis.size()));

  const GridSymmetry symmetry(SurfaceLattice(d.system), PointGroup(d.system.point_group), p.nx, p.ny,
                              s.angular.n_phi());
  ASSERT_EQ(symmetry.operations().size(), 4U);
  const auto c2 = std::ranges::find(symmetry.operations(), SurfaceOperation{false, 6});
  ASSERT_NE(c2, symmetry.operations().end());
  const int op = static_cast<int>(c2 - symmetry.operations().begin());
  const double asymmetry = density_asymmetry(r, symmetry, op, 0);
  const double shifted = density_asymmetry(r, symmetry, op, 1);
  EXPECT_LE(asymmetry, 1e-12);
  EXPECT_GT(shifted, 1e-2);
  double mirrors = 0.0;  // printed only: not invariant on this grid
  for (std::size_t k = 0; k < symmetry.operations().size(); ++k)
    if (symmetry.operations()[k].mirror)
      mirrors = nan_aware_max(mirrors, density_asymmetry(r, symmetry, static_cast<int>(k), 0));
  std::cout << "  C2: max |P(-x, -y) - P(x, y)| / max |P| = " << asymmetry << "; with ix shifted by one " << shifted
            << "; the mirrors " << mirrors << '\n';
}
