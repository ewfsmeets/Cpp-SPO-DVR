// Expected memory of the large arrays of a run, from the grid sizes and the allocation rules
// of the propagation operators (WavePacketPropagation and its stages).
#include "propagation/memory_estimate.hpp"

#include <algorithm>
#include <complex>

#include "grids/rotational_basis.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"

namespace spodvr {

namespace {

constexpr std::size_t real_bytes = sizeof(double);
constexpr std::size_t complex_bytes = sizeof(std::complex<double>);

}  // namespace

MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                               int n_threads) {
  return estimate_memory(config, grids, angular, make_grid_symmetry(config, angular), n_threads);
}

MemoryEstimate estimate_memory(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                               const GridSymmetry& symmetry, int n_threads) {
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
  if (reduced) {
    // GridSymmetry: representative and operation per (x, y) point, a phi map per operation.
    const std::size_t map_values = 2 * n_xy + symmetry.operations().size() * std::size_t(angular.n_phi());
    estimate.items.push_back(
        {"Surface symmetry maps", map_values, map_values * sizeof(int),
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
  }

  for (const MemoryItem& item : estimate.items) {
    if (item.during_setup) estimate.setup_bytes = estimate.setup_bytes + item.bytes;
    if (item.during_propagation) estimate.propagation_bytes = estimate.propagation_bytes + item.bytes;
  }
  estimate.peak_bytes = std::max(estimate.setup_bytes, estimate.propagation_bytes);
  return estimate;
}

}  // namespace spodvr
