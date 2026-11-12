// Tests: the flux analysis in the time loop with the GPU backend: the coefficients are the sums of the
// flux surfaces the GPU hands out, whatever the slot count, the flux analysis leaves the rest of a GPU run
// unchanged, and the coefficients agree with the CPU run's to rounding (CPU twin: test_flux_propagation.cpp).
#include <gtest/gtest.h>

#include <complex>
#include <cstddef>
#include <vector>

#include "backend_setup.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda_test_setup.hpp"
#include "flux_propagation_run.hpp"

using namespace spodvr;
using namespace spodvr::testing;

namespace {

const BackendSetup& flux_setup() {
  static const BackendSetup s(flux_run_case, flux_run_input(0.0, 120.0, 60));
  return s;
}

// The GPU run of these tests, with the budget's slots.
const FluxPropagationRun& gpu_budget_run() {
  static const FluxPropagationRun run =
      run_flux_propagation(flux_setup(), flux_budget_slots(flux_setup()), cuda::make_cuda_backend);
  return run;
}

}  // namespace

// As on the CPU: the coefficients after run() are the sums over all 12 analysis steps of c_n(t) times the
// flux surfaces the GPU backend handed to the observer, recomputed here, bitwise; 8 slots for 12 steps
// (a flush when full, 4 steps at the end of run()); nothing pending after it.
TEST(CudaFlux, CoefficientsAreTheSumsOfTheGpuFluxSurfaces) {
  SPODVR_REQUIRE_GPU();
  const FluxPropagationRun& run = gpu_budget_run();
  ASSERT_EQ(run.flux_times.size(), 12U);
  EXPECT_EQ(run.flux_times, run.bk_times);
  EXPECT_EQ(run.slot_count, 8);
  EXPECT_EQ(run.pending_after_run, 0);
  ASSERT_TRUE(run.coefficients && run.recomputed);
  EXPECT_EQ(bit_differences(*run.coefficients, *run.recomputed), 0);
}

// One slot (each GPU flux surface flushed at once) gives bitwise the coefficients of 8.
TEST(CudaFlux, CoefficientsDoNotDependOnTheSlotCount) {
  SPODVR_REQUIRE_GPU();
  const FluxPropagationRun run = run_flux_propagation(flux_setup(), 1, cuda::make_cuda_backend);
  ASSERT_TRUE(run.coefficients);
  EXPECT_EQ(bit_differences(*run.coefficients, *gpu_budget_run().coefficients), 0);
}

// The flux surface reuses the GPU read-back buffer of the BK plane; it must not change the BK
// coefficients or the averages of a GPU run, which are bitwise those of the GPU run without flux analysis.
TEST(CudaFlux, TheFluxAnalysisIsPassive) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup without(flux_run_case, flux_run_input(0.0, 120.0, 60, false));
  const FluxPropagationRun run = run_flux_propagation(without, 1, cuda::make_cuda_backend);
  const FluxPropagationRun& with = gpu_budget_run();
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

// The GPU and the CPU run of the same input: the flux surfaces agree bitwise on the same psi
// (CudaBackend.FluxSurfaceMatchesTheCpu) and the sums are the same host code, so the coefficients
// differ only by the propagation's rounding (cuFFT against the CPU FFT library). Tolerance 1e-12 of
// max |phi|: the per-step bound of the backend tests (1e-13, AnalysisStepMatchesTheCpu*) over the 12
// steps would allow 1.2e-12; measured 1.1e-14 (phi) and 1.5e-14 (phi'), 2026-11-12.
TEST(CudaFlux, CoefficientsMatchTheCpu) {
  SPODVR_REQUIRE_GPU();
  const FluxPropagationRun cpu = run_flux_propagation(flux_setup(), flux_budget_slots(flux_setup()));
  const FluxPropagationRun& gpu = gpu_budget_run();
  ASSERT_TRUE(cpu.coefficients && gpu.coefficients);
  const double phi = relative_difference(cpu.coefficients->phi, gpu.coefficients->phi);
  const double phi_derivative = relative_difference(cpu.coefficients->phi_derivative, gpu.coefficients->phi_derivative);
  EXPECT_LE(phi, 1e-12);
  EXPECT_LE(phi_derivative, 1e-12);
}
