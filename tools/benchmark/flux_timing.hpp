#pragma once
// The cost of the flux analysis per analysis step, shared by spodvr_operator_timing (CPU) and
// spodvr_cuda_operator_timing (GPU): one flush of the accumulator on the host per slot count, and the
// flux surface plus the flush as a share of an analysis step.

#include <omp.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <optional>
#include <vector>

#include "analysis/flux_analysis.hpp"
#include "propagation/memory_estimate.hpp"

namespace spodvr::benchmark {

// The slot counts whose flush is timed: one flush per analysis step up to twice max_flux_slots.
inline constexpr std::array<int, 5> flush_slot_counts{1, 2, 4, 8, 16};

// The flush with slot_count slots, per analysis step (one flush adds slot_count steps).
struct FlushTiming {
  int slot_count;
  double ms_per_step;
};

// For each of flush_slot_counts: fills every slot once with fill (the last fill flushes), then times
// the flush of slot_count steps (one warm-up, then `repetitions`), on the host's OpenMP threads.
// Each flush adds the same slot contents again; the cost does not depend on the values.
inline std::vector<FlushTiming> time_flux_flushes(int n_blocks, int n_planes, const FluxEnergies& energies,
                                                  double time_step, int repetitions,
                                                  const std::function<void(FluxSurfaceWave&)>& fill) {
  std::vector<FlushTiming> timings;
  for (const int slots : flush_slot_counts) {
    FluxAccumulator flux(n_blocks, n_planes, energies, time_step, slots);
    for (int k = 0; k < slots; ++k) flux.add_step(time_step * k, fill);
    const auto flush = [&] {
      for (int k = 0; k < slots; ++k) flux.add_step(time_step * k, [](FluxSurfaceWave&) {});
    };
    flush();
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repetitions; ++i) flush();
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    timings.push_back({slots, elapsed.count() / repetitions / slots * 1000.0});
  }
  return timings;
}

// Prints the flux cost per analysis step for each slot count: the flux surface (surface_ms) plus the
// flush, as a share of an analysis step of analysis_step_ms, and the slots the program uses for
// slots of slot_bytes when the physical memory does not limit them (plan_flux_memory).
inline void print_flux_cost(int n_energies, double surface_ms, const std::vector<FlushTiming>& flushes,
                            double analysis_step_ms, std::size_t slot_bytes) {
  std::printf(
      "Flux analysis per analysis step (%d energies, flush on %d host threads; analysis step %.1f ms; one slot "
      "%.1f MB, the program uses %d):\n",
      n_energies, omp_get_max_threads(), analysis_step_ms, double(slot_bytes) / 1e6,
      plan_flux_memory(0, 0, slot_bytes, std::nullopt).slot_count);
  std::printf("  %6s %16s %24s %12s\n", "slots", "flush (ms/step)", "surface + flush (ms)", "of the step");
  for (const FlushTiming& f : flushes) {
    const double sum = surface_ms + f.ms_per_step;
    std::printf("  %6d %16.1f %24.1f %11.2f %%\n", f.slot_count, f.ms_per_step, sum, 100.0 * sum / analysis_step_ms);
  }
}

}  // namespace spodvr::benchmark
