#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "symmetry/grid_symmetry.hpp"

namespace spodvr {

// One large array of a run.
struct MemoryItem {
  std::string name;
  std::size_t values;  // number of values (per thread for the work arrays)
  std::size_t bytes;   // total size (all threads for the work arrays)
  std::string note;
  bool during_setup;        // allocated while the propagation is set up
  bool during_propagation;  // allocated during the time loop
};

// Expected memory of the large arrays, computed from the input before anything is allocated.
// Small tables (kinetic factors, rotational energies, analysis basis, S-matrix) and the FFT
// plans are not counted.
struct MemoryEstimate {
  std::vector<MemoryItem> items;
  std::size_t setup_bytes = 0;        // while the propagation is set up (potential and exp table together)
  std::size_t propagation_bytes = 0;  // during the time loop
  std::size_t peak_bytes = 0;         // the larger of the two
  int flux_slot_count = 0;            // slots of the flux accumulator (plan_flux_memory); 0 without flux analysis
};

// Host memory for the slots of the flux accumulator, which collect the flux surfaces of that many
// analysis steps before they are added to the coefficients. More slots, fewer passes over the
// coefficients: on the DZ grid (one slot 250.6 MB) one flush of 22 energies took 387 ms per analysis
// step with 1 slot and 79.6 ms with 8 (docs/cuda.md, "Flux surface and flush"). 2 GiB holds those 8.
inline constexpr std::size_t flux_slot_budget_bytes = std::size_t(2) << 30;

// At most this many slots, the most measured above: the budget alone would allow thousands on small
// grids (hardwall_small: 0.87 MB per slot), where one flush per step costs little anyway.
inline constexpr int max_flux_slots = 8;

// The memory plan of the flux analysis: how many slots, and whether the run fits.
struct FluxMemoryPlan {
  int slot_count;  // >= 1
  bool fits;       // false only if the flux analysis is the cause: base fits, base + coefficients + 1 slot does not
};

// The slots for coefficients of coefficient_bytes (phi and phi') and slots of slot_bytes (> 0) each,
// next to base bytes of everything else that is allocated during the time loop: as many as
// flux_slot_budget_bytes holds, 1 to max_flux_slots, fewer if they would not fit in the physical memory
// (if known). Device-agnostic: the caller computes base for its device.
[[nodiscard]] FluxMemoryPlan plan_flux_memory(std::size_t base, std::size_t coefficient_bytes, std::size_t slot_bytes,
                                              std::optional<std::size_t> physical_bytes);

// The arrays of the run mode in config: the potential grid only for a potential-only run (6),
// otherwise also the exp(-i V dt) table, the wave functions and the per-thread work arrays of
// the propagation operators (n_threads OpenMP threads), and with the flux analysis its coefficients
// and slots (allocated last; flux_slot_count from plan_flux_memory without a physical-memory limit).
// The potential and the exp table hold the (x, y) representatives of the surface symmetry only; the
// symmetry maps are listed when they reduce anything.
[[nodiscard]] MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, const GridSymmetry& symmetry, int n_threads);

// The same with the symmetry of the run (make_grid_symmetry).
[[nodiscard]] MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, int n_threads);

}  // namespace spodvr
