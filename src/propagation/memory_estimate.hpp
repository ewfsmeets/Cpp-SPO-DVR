#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"

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
};

// The arrays of the run mode in config: the potential grid only for a potential-only run (6),
// otherwise also the exp(-i V dt) table, the wave functions and the per-thread work arrays of
// the propagation operators (n_threads OpenMP threads).
[[nodiscard]] MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, int n_threads);

}  // namespace spodvr
