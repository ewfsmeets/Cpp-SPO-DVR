#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/propagation_backend.hpp"
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

// Host memory of a run without the flux analysis: while the propagation is set up, and during the
// time loop, next to which the flux arrays are allocated (they do not exist during the set-up).
struct HostMemoryBase {
  std::size_t setup = 0;
  std::size_t time_loop = 0;
};

// The arrays of the flux analysis (host memory, allocated after the set-up).
struct FluxMemoryItems {
  std::size_t bytes_per_energy = 0;  // phi and phi' of one energy
  int n_energies = 0;
  std::size_t slot_bytes = 0;  // psi and d psi / dr of one analysis step

  [[nodiscard]] std::size_t coefficient_bytes() const { return bytes_per_energy * n_energies; }
};

// The memory plan of the flux analysis: how many slots, and whether the run fits.
struct FluxMemoryPlan {
  int slot_count;  // >= 1
  bool fits;       // false only if the flux analysis is the cause: base fits, base + coefficients + 1 slot does not
};

// The memory rule of the flux analysis (check_flux_memory) with every number it used, so that the
// printed numbers and the error (require_flux_memory) agree.
struct FluxMemoryCheck {
  Device device;
  HostMemoryBase base;  // the run without the flux analysis on this device
  FluxMemoryItems items;
  std::optional<std::size_t> physical_bytes;
  bool base_fits;       // the run fits without the flux analysis (false if the physical memory is unknown)
  FluxMemoryPlan plan;  // the slot count; fits is false only when the flux analysis is the cause
  // The most flux energies that fit with one slot next to base.time_loop; empty unless base_fits.
  std::optional<std::size_t> max_energies;
};

// Expected memory of the large arrays, computed from the input before anything is allocated.
// Small tables (kinetic factors, rotational energies, analysis basis, S-matrix) and the FFT
// plans are not counted.
struct MemoryEstimate {
  std::vector<MemoryItem> items;
  std::size_t setup_bytes = 0;        // while the propagation is set up (potential and exp table together)
  std::size_t propagation_bytes = 0;  // during the time loop
  std::size_t peak_bytes = 0;         // the larger of the two
  // The host memory of a GPU run (--device cuda) without the flux analysis, which differs from the
  // items above (they are a CPU run's): during the set-up V and exp(i A) (CudaBackend's constructor;
  // exp_potential_device_order in DevicePotentialOperator's constructor, at any table size), the
  // symmetry maps, the specular waves and the read-back plane (CudaBackend::z_plane, allocated in the
  // constructor); in the time loop the host copy of the wave packet (CudaBackend::wave_packet_on_host,
  // CudaBackend::read_wave_packet), the specular waves and the read-back plane. Zero for a
  // potential-only run.
  HostMemoryBase gpu_host;
  // For a GPU run (estimate_memory with Device::Cuda): its peak host memory, the larger of
  // gpu_host.setup and gpu_host.time_loop plus the flux arrays. The memory warning uses it instead of
  // peak_bytes, which is a CPU run's. Empty otherwise.
  std::optional<std::size_t> gpu_host_peak_bytes;
  std::optional<FluxMemoryCheck> flux_memory;  // with the flux analysis: its memory rule

  // The slot count of the flux accumulator (flux_memory's plan); 0 without the flux analysis.
  [[nodiscard]] int flux_slot_count() const { return flux_memory ? flux_memory->plan.slot_count : 0; }
};

// Host memory for the slots of the flux accumulator, which collect the flux surfaces of that many
// analysis steps before they are added to the coefficients. More slots, fewer passes over the
// coefficients: on the DZ grid (one slot 250.6 MB, 22 energies, 8 threads) the flush took 390-408 ms
// per analysis step with 1 slot, 60-75 ms with 8 and 43-56 ms with 16 (docs/performance.md, "Flux
// analysis cost on the DZ grid"). 2 GiB holds those 8; 16 would take 2 GB more and still not bring a
// GPU run's flux analysis below 2 % of its time.
inline constexpr std::size_t flux_slot_budget_bytes = std::size_t(2) << 30;

// At most this many slots, what the budget holds on DZ: the budget alone would allow thousands on
// small grids (hardwall_small: 0.87 MB per slot), where one flush per step costs little anyway.
inline constexpr int max_flux_slots = 8;

// The slots for coefficients of coefficient_bytes (phi and phi') and slots of slot_bytes (> 0) each,
// next to base bytes of everything else that is allocated during the time loop: as many as
// flux_slot_budget_bytes holds, 1 to max_flux_slots, fewer if they would not fit in the physical memory
// (if known). Device-agnostic: the caller computes base for its device.
[[nodiscard]] FluxMemoryPlan plan_flux_memory(std::size_t base, std::size_t coefficient_bytes, std::size_t slot_bytes,
                                              std::optional<std::size_t> physical_bytes);

// The memory rule of the flux analysis for a run on `device` whose estimate without the flux analysis
// is `without_flux` (CPU: its set-up and time-loop totals; CUDA: its gpu_host). The flux arrays exist
// only in the time loop, so the flux analysis is the cause of a shortage only if the run fits without
// it (max(setup, time_loop) <= physical) and not with its coefficients and one slot next to time_loop.
// Then plan.fits is false; otherwise the slots are plan_flux_memory's, within the physical memory.
// Unknown physical memory, or a run too large without the flux analysis (the estimate warns): the
// budget's slots, and fits.
[[nodiscard]] FluxMemoryCheck check_flux_memory(const MemoryEstimate& without_flux, const FluxMemoryItems& items,
                                                std::optional<std::size_t> physical_bytes, Device device);

// Throws InputError (DeltaEFluxAna) with the numbers of the check if the flux analysis does not fit
// (plan.fits false): its coefficients and one slot exceed what the run leaves of the physical memory.
void require_flux_memory(const FluxMemoryCheck& check);

// The arrays of the run mode in config: the potential grid only for a potential-only run (6),
// otherwise also the exp(-i V dt) table, the wave functions and the per-thread work arrays of
// the propagation operators (n_threads OpenMP threads), and with the flux analysis its coefficients
// and slots (allocated last; the slot count from check_flux_memory for `device` and physical_bytes).
// The potential and the exp table hold the (x, y) representatives of the surface symmetry only; the
// symmetry maps are listed when they reduce anything. The items are a CPU run's; gpu_host holds the
// host memory of a GPU run.
[[nodiscard]] MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, const GridSymmetry& symmetry, int n_threads,
                                             Device device, std::optional<std::size_t> physical_bytes);

// A CPU run with the symmetry of the run (make_grid_symmetry) and an unknown physical memory (the
// flux slots of the budget).
[[nodiscard]] MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, int n_threads);

}  // namespace spodvr
