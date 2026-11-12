// Tests: the flux analysis in the time loop (WavePacketPropagation): its coefficients are the sums of the
// flux surfaces of every analysis step, whatever the slot and thread count, a restart continues them
// bitwise, and the rest of the run is the same with the flux analysis on.
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "backend_setup.hpp"
#include "flux_propagation_run.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"

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
