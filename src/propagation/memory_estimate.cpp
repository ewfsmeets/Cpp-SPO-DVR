// Expected memory of the large arrays of a run, from the grid sizes and the allocation rules
// of the propagation operators (WavePacketPropagation and its stages).
#include "propagation/memory_estimate.hpp"

#include <algorithm>
#include <complex>
#include <format>
#include <stdexcept>

#include "analysis/flux_analysis.hpp"
#include "core/error.hpp"
#include "grids/rotational_basis.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"
#include "report/format.hpp"

namespace spodvr {

namespace {

constexpr std::size_t real_bytes = sizeof(double);
constexpr std::size_t complex_bytes = sizeof(std::complex<double>);

// The totals of the items while the propagation is set up, during the time loop, and their maximum.
void add_up(MemoryEstimate& estimate) {
  estimate.setup_bytes = 0;
  estimate.propagation_bytes = 0;
  for (const MemoryItem& item : estimate.items) {
    if (item.during_setup) estimate.setup_bytes = estimate.setup_bytes + item.bytes;
    if (item.during_propagation) estimate.propagation_bytes = estimate.propagation_bytes + item.bytes;
  }
  estimate.peak_bytes = std::max(estimate.setup_bytes, estimate.propagation_bytes);
}

// The coefficients phi, phi' of the flux analysis and its slots, each with psi and d psi / dr per block
// (n_blocks = n_rot nx ny) and plane below Z_asymptotic. Allocated after the set-up, next to the other
// arrays of the time loop; the slot count from the memory rule. Fortran: FluxAnalyzor.f90 phi, phiDeriv
void add_flux_items(MemoryEstimate& estimate, const Config& config, const SpatialGrids& grids, std::size_t n_blocks,
                    Device device, std::optional<std::size_t> physical_bytes) {
  const std::size_t surface_values = 2 * n_blocks * grids.nz_potential;
  const FluxMemoryItems flux{surface_values * complex_bytes, flux_energy_count(config.input),
                             surface_values * complex_bytes};
  const FluxMemoryCheck check = check_flux_memory(estimate, flux, physical_bytes, device);
  estimate.flux_memory = check;
  const std::size_t coefficient_values = surface_values * flux.n_energies;
  estimate.items.push_back({"Flux coefficients", coefficient_values, coefficient_values * complex_bytes,
                            std::to_string(flux.n_energies) + " energies", false, true});
  const std::size_t slot_values = surface_values * check.plan.slot_count;
  estimate.items.push_back({"Flux surface slots", slot_values, slot_values * complex_bytes,
                            std::to_string(check.plan.slot_count) + " analysis steps per flush", false, true});
}

// "1.5 GiB (1610612736 bytes)": the size as printed, and exact.
std::string bytes_text(std::size_t bytes) { return std::format("{} ({} bytes)", format_bytes(bytes), bytes); }

}  // namespace

FluxMemoryPlan plan_flux_memory(std::size_t base, std::size_t coefficient_bytes, std::size_t slot_bytes,
                                std::optional<std::size_t> physical_bytes) {
  if (slot_bytes == 0) throw std::invalid_argument("plan_flux_memory: slots of 0 bytes");
  const int budget_slots =
      static_cast<int>(std::clamp(flux_slot_budget_bytes / slot_bytes, std::size_t(1), std::size_t(max_flux_slots)));
  // Unknown memory, or a run too large without the flux analysis too (the estimate warns about that).
  if (!physical_bytes || base > *physical_bytes) return {budget_slots, true};
  const std::size_t free_bytes = *physical_bytes - base;
  if (coefficient_bytes > free_bytes || free_bytes - coefficient_bytes < slot_bytes) return {1, false};
  const std::size_t fitting_slots = (free_bytes - coefficient_bytes) / slot_bytes;
  return {static_cast<int>(std::min(fitting_slots, std::size_t(budget_slots))), true};
}

FluxMemoryCheck check_flux_memory(const MemoryEstimate& without_flux, const FluxMemoryItems& items,
                                  std::optional<std::size_t> physical_bytes, Device device) {
  if (items.bytes_per_energy == 0) throw std::invalid_argument("check_flux_memory: flux energies of 0 bytes");
  const HostMemoryBase base = device == Device::Cuda
                                  ? without_flux.gpu_host
                                  : HostMemoryBase{without_flux.setup_bytes, without_flux.propagation_bytes};
  const bool base_fits = physical_bytes && std::max(base.setup, base.time_loop) <= *physical_bytes;
  // With a run too large anyway the physical memory limits nothing (plan_flux_memory without it).
  const FluxMemoryPlan plan = plan_flux_memory(base.time_loop, items.coefficient_bytes(), items.slot_bytes,
                                               base_fits ? physical_bytes : std::nullopt);
  std::optional<std::size_t> max_energies;
  if (base_fits) {
    const std::size_t free_bytes = *physical_bytes - base.time_loop;
    max_energies = free_bytes < items.slot_bytes ? 0 : (free_bytes - items.slot_bytes) / items.bytes_per_energy;
  }
  return {device, base, items, physical_bytes, base_fits, plan, max_energies};
}

void require_flux_memory(const FluxMemoryCheck& check) {
  if (check.plan.fits) return;
  const FluxMemoryItems& f = check.items;
  const std::size_t physical = check.physical_bytes.value_or(0);  // known: plan.fits is false only then
  const std::size_t free_bytes = physical - check.base.time_loop;
  throw InputError(
      "DeltaEFluxAna",
      std::format("the flux analysis needs {} during the time loop ({} energies of {} and one slot of {}), but "
                  "only {} of the physical memory are left ({} minus {} for the rest of the run{}); at most {} "
                  "energies fit: choose a larger DeltaEFluxAna or a narrower EMinFluxAna .. EMaxFluxAna",
                  bytes_text(f.coefficient_bytes() + f.slot_bytes), f.n_energies, format_bytes(f.bytes_per_energy),
                  format_bytes(f.slot_bytes), bytes_text(free_bytes), format_bytes(physical),
                  format_bytes(check.base.time_loop), check.device == Device::Cuda ? ", host memory of a GPU run" : "",
                  check.max_energies.value_or(0)));
}

MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                               int n_threads) {
  return estimate_memory(config, grids, angular, make_grid_symmetry(config, angular), n_threads, Device::Cpu,
                         std::nullopt);
}

MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                               const GridSymmetry& symmetry, int n_threads, Device device,
                               std::optional<std::size_t> physical_bytes) {
  const InputParameters& p = config.input;
  const std::size_t n_xy = std::size_t(p.nx) * p.ny;
  const std::size_t n_tp = std::size_t(angular.n_theta()) * angular.n_phi();
  const bool reduced = symmetry.n_representatives() < symmetry.n_xy();
  MemoryEstimate estimate;

  // Potential on the interaction region (Z < zasymp), at the (x, y) representatives. Fortran: vdvr
  const std::size_t potential_values = n_tp * symmetry.n_representatives() * grids.nz_potential * p.nr;
  const bool precomputed = potential_values * complex_bytes <= PotentialOperator::default_table_limit;
  const bool propagation = config.derived.run_mode != 6;
  std::string potential_note;
  if (propagation) potential_note = precomputed ? "freed after setup" : "kept";
  estimate.items.push_back({"Potential V", potential_values, potential_values * real_bytes, potential_note, true,
                            propagation && !precomputed});
  std::size_t map_bytes = 0;
  if (reduced) {
    // GridSymmetry: representative and operation per (x, y) point, a phi map per operation.
    const std::size_t map_values = 2 * n_xy + symmetry.operations().size() * std::size_t(angular.n_phi());
    map_bytes = map_values * sizeof(int);
    estimate.items.push_back(
        {"Surface symmetry maps", map_values, map_bytes,
         std::to_string(symmetry.n_representatives()) + " of " + std::to_string(n_xy) + " (x, y) points", true,
         propagation});
  }

  if (propagation) {
    if (precomputed)
      estimate.items.push_back({"exp(-i V dt) table", potential_values, potential_values * complex_bytes,
                                "precomputed (at most 4 GiB)", true, true});
    else
      estimate.items.push_back({"exp(-i V dt) table", 0, 0, "computed in every step (more than 4 GiB)", true, true});

    const RotationalBasis basis(angular, p.j_in, p.mj_in);
    const std::size_t scattering_values = std::size_t(basis.size()) * n_xy * p.nz * p.nr;
    estimate.items.push_back(
        {"Wave function, scattering grid", scattering_values, scattering_values * complex_bytes, "", true, true});
    // The propagated specular wave plus the initial wave packet.
    const std::size_t specular_values = std::size_t(p.nr) * grids.z_specular.size();
    estimate.items.push_back(
        {"Wave function, specular grid", specular_values, 2 * specular_values * complex_bytes, "2 copies", true, true});

    // Per thread: the (r, Z) unit of T_Z + T_r, or the (x, y) group of T_xy together with the
    // DVR blocks of the potential step (the stages run one after the other).
    const ScatteringShape shape{basis.size(), p.nx, p.ny, p.nr, p.nz};
    const std::size_t zr_unit = std::size_t(zr_chunk_size(shape.n_rot * shape.ny)) * p.nr * p.nz;
    const std::size_t xy_group = std::size_t(shape.n_rot) * n_xy * xy_group_points(shape);
    const std::size_t dvr_blocks =
        std::size_t(AngularTransform::dvr_block_size(angular.n_theta(), angular.n_phi())) * n_xy;
    // On the fly with symmetry: exp(i A) of the representatives of one (Z, r) point.
    const std::size_t exp_scratch = !precomputed && reduced ? n_tp * symmetry.n_representatives() : 0;
    const std::size_t work_values = std::max(zr_unit, xy_group + dvr_blocks + exp_scratch);
    estimate.items.push_back({"Work arrays", work_values, work_values * complex_bytes * n_threads,
                              "per thread, x " + std::to_string(n_threads) + " threads", false, true});

    // A GPU run's host memory (MemoryEstimate::gpu_host): V with exp(i A) during the set-up, the host
    // copy of psi and a third specular wave in the time loop; the read-back plane (n_rot nx ny nr) in both.
    const std::size_t specular_bytes = 2 * specular_values * complex_bytes;  // the item above
    const std::size_t plane_bytes = std::size_t(basis.size()) * n_xy * p.nr * complex_bytes;
    estimate.gpu_host.setup =
        potential_values * (real_bytes + complex_bytes) + map_bytes + specular_bytes + plane_bytes;
    estimate.gpu_host.time_loop =
        scattering_values * complex_bytes + specular_bytes + specular_values * complex_bytes + plane_bytes;
    add_up(estimate);
    if (p.do_flux_analysis)
      add_flux_items(estimate, config, grids, std::size_t(basis.size()) * n_xy, device, physical_bytes);
  }

  add_up(estimate);
  if (propagation && device == Device::Cuda) {
    const std::size_t flux_bytes =
        estimate.flux_memory ? estimate.flux_memory->items.coefficient_bytes() +
                                   estimate.flux_memory->plan.slot_count * estimate.flux_memory->items.slot_bytes
                             : 0;
    estimate.gpu_host_peak_bytes = std::max(estimate.gpu_host.setup, estimate.gpu_host.time_loop + flux_bytes);
  }
  return estimate;
}

}  // namespace spodvr
