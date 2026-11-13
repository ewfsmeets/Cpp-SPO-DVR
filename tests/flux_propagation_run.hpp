#pragma once
// Shared code of the tests of the flux analysis in the time loop (tests/test_flux_propagation.cpp, and
// tests/cuda/test_cuda_flux.cpp with the GPU backend): a short run of a reference case through
// WavePacketPropagation from a test wave packet, and the time-to-energy sums recomputed here from the
// flux surfaces the run hands to its observer. Include it only in test programs linked with
// spodvr_strict_fp: the sums here must round as the flux analysis's own (no FMA), also in build-fast.

#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "backend_setup.hpp"
#include "propagation/averages.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"

namespace spodvr::testing {

// The reference case of these tests and its flux surface: hardwall_short at r_c = 2.0 bohr (ir 8),
// energies EMinWP .. EMaxWP = 0.2 .. 0.6 eV by 0.15 eV (3 energies, flux Q3).
inline const char* const flux_run_case = "hardwall_short";

// hardwall_short from tStart to tEnd in n_steps time steps of 2 a.u. (one analysis step: 10 a.u.),
// with the flux analysis unless `flux` is false.
inline nlohmann::json flux_run_input(double t_start, double t_end, int n_steps, bool flux = true) {
  nlohmann::json changes = flux ? flux_analysis_keys(flux_run_case, 2.0, 0.15) : nlohmann::json::object();
  changes["tStart"] = t_start;
  changes["tEnd"] = t_end;
  changes["nStep"] = n_steps;
  return changes;
}

// The slot count of the run's memory estimate (the flux slot budget) for the set-up.
inline int flux_budget_slots(const BackendSetup& s) {
  return estimate_memory(s.config, s.grids, s.angular, 1).flux_slot_count();
}

// Replaces the initial wave packet of a propagation by distinct values of varying magnitude and phase
// everywhere (as a restart reads one): not a physical state, on purpose. From the reference case's own
// state the scattering grid holds only rounding noise in a short run (hardwall_short: |psi| <= 1.2e-29
// after 120 a.u., Dai-Light <= 2e-13; measured 2026-11-12), and so would the flux surface; here it holds
// values of order 1, which a comparison of two runs can tell apart.
inline void start_from_test_wave_packet(WavePacketPropagation& propagation) {
  auto fill = [](std::span<std::complex<double>> values, double phase) {
    for (std::size_t i = 0; i < values.size(); ++i) {
      const double x = static_cast<double>(i);
      values[i] = std::complex<double>(std::sin(0.37 * x + phase), std::cos(1.31 * x + phase)) *
                  (1.0 + 0.5 * std::sin(0.011 * x));
    }
  };
  propagation.read_wave_packet([&](WavePacket& psi) {
    fill(psi.scattering.values, 0.1);
    fill(psi.specular.values, 0.7);
  });
}

// phi(block, iz, n) += c_n(t) psi(block, iz) and phi' += c_n(t) d psi / dr, written out here in the
// Fortran's order (one addition per step), with c_n(t) = flux_step_coefficient at the weight time_step.
// c_n comes from the same function as in the run (its own tests: test_flux_analysis.cpp): a comparison
// with these sums checks the time loop's bookkeeping (steps, times, weights, slots, flushes), not c_n.
// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF
inline void add_flux_surface_by_hand(FluxCoefficients& sums, const FluxEnergies& energies, double time,
                                     double time_step, const FluxSurfaceWave& surface) {
  for (int n = 0; n < energies.size(); ++n) {
    const std::complex<double> c = flux_step_coefficient(energies, n, time, time_step);
    for (int iz = 0; iz < surface.n_planes(); ++iz)
      for (int block = 0; block < surface.n_blocks(); ++block) {
        const std::size_t k = sums.index(block, iz, n);
        sums.phi[k] = sums.phi[k] + c * surface.psi[surface.index(block, iz)];
        sums.phi_derivative[k] = sums.phi_derivative[k] + c * surface.derivative[surface.index(block, iz)];
      }
  }
}

// What a run hands out: per analysis step its time, BK coefficients and (flux analysis) flux
// surface; the averages; and after run() the flux coefficients, their recomputation from the
// observed surfaces, the steps still pending and the slot count.
struct FluxPropagationRun {
  std::vector<double> bk_times;
  std::vector<std::vector<std::complex<double>>> bk;
  std::vector<Averages> averages;
  std::vector<double> flux_times;
  std::optional<FluxCoefficients> coefficients;  // WavePacketPropagation::flux() after run()
  std::optional<FluxCoefficients> recomputed;    // add_flux_surface_by_hand over the observed surfaces
  int pending_after_run = -1;
  int slot_count = 0;
};

// Runs the propagation of the set-up from the test wave packet with flux_slot_count slots on `backend`.
// The recomputation uses energies built here from the set-up and the weight derived.analysis_time_step
// (flux Q4).
inline FluxPropagationRun run_flux_propagation(const BackendSetup& s, int flux_slot_count,
                                               const BackendFactory& backend = make_cpu_backend) {
  WavePacketPropagation propagation(s.config, s.grids, s.angular, s.potential, s.gas, s.fft, backend, flux_slot_count);
  start_from_test_wave_packet(propagation);
  const FluxEnergies energies = make_flux_energies(s.config, s.tables.parallel_energy, s.gas.initial.energy);
  FluxPropagationRun run;
  PropagationObserver observer;
  observer.on_averages = [&](const Averages& a) { run.averages.push_back(a); };
  observer.on_analysis_step = [&](double time, const BkCoefficients& bk) {
    run.bk_times.push_back(time);
    run.bk.push_back(bk.values);
  };
  observer.on_flux_surface = [&](double time, const FluxSurfaceWave& surface) {
    if (!run.recomputed) run.recomputed.emplace(surface.n_blocks(), surface.n_planes(), energies.size());
    run.flux_times.push_back(time);
    add_flux_surface_by_hand(*run.recomputed, energies, time, s.config.derived.analysis_time_step, surface);
  };
  propagation.run(observer);
  if (propagation.flux()) {
    run.pending_after_run = propagation.flux()->pending_steps();
    run.slot_count = propagation.flux()->slot_count();
    run.coefficients = propagation.flux()->coefficients();
  }
  return run;
}

// Whether x is a NaN or an infinity, read from its bits (all exponent bits set), in every build.
inline bool non_finite_bits(double x) {
  constexpr std::uint64_t exponent = 0x7ff0000000000000ULL;
  return (std::bit_cast<std::uint64_t>(x) & exponent) == exponent;
}

// Two doubles differ unless their bit patterns are equal and finite: a NaN or infinity on either side
// is always a difference (also the same NaN on both sides), and so is -0 against +0.
inline bool doubles_differ(double x, double y) {
  return non_finite_bits(x) || non_finite_bits(y) || std::bit_cast<std::uint64_t>(x) != std::bit_cast<std::uint64_t>(y);
}

// Doubles of a and b that differ (doubles_differ; all of them if the sizes differ).
inline int bit_differences(std::span<const std::complex<double>> a, std::span<const std::complex<double>> b) {
  if (a.size() != b.size()) return static_cast<int>(2 * std::max(a.size(), b.size()));
  int differences = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (doubles_differ(a[i].real(), b[i].real())) ++differences;
    if (doubles_differ(a[i].imag(), b[i].imag())) ++differences;
  }
  return differences;
}

// Bit differences of phi and phi' together.
inline int bit_differences(const FluxCoefficients& a, const FluxCoefficients& b) {
  return bit_differences(a.phi, b.phi) + bit_differences(a.phi_derivative, b.phi_derivative);
}

}  // namespace spodvr::testing
