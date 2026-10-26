// Per-operator timing of the propagator for an input file: the (r, Z) kinetic step, the (x, y)
// kinetic step, the potential step, the specular grid, and one analysis step (all operators).
// Usage: spodvr_operator_timing <input.json> [repetitions]   (threads: OMP_NUM_THREADS)
#include <omp.h>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>

#include <nlohmann/json.hpp>

#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "propagation/propagation.hpp"

using namespace spodvr;

namespace {

// Runs `step` once (warm-up), then `repetitions` times, and prints the mean time per call.
void time_call(const char* name, int repetitions, const std::function<void()>& step) {
  step();
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < repetitions; ++i) step();
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  std::printf("  %-28s %10.1f ms per call\n", name, elapsed.count() / repetitions * 1000.0);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <input.json> [repetitions]\n", argv[0]);
    return 1;
  }
  const int repetitions = argc > 2 ? std::stoi(argv[2]) : 3;
  std::ifstream file(argv[1]);
  const Config config = build_config(parse_input(nlohmann::json::parse(file)));
  const InputParameters& p = config.input;

  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular(p.j_max, config.derived.mj_max_basis, p.homonuclear, p.j_in, p.mj_in);
  PotentialGrid potential = compute_potential_grid(*make_potential(p.system_id), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());

  const RotationalBasis basis(angular, p.j_in, p.mj_in);
  const AngularTransform transform(angular, basis, fft, p.nx * p.ny);
  const PropagationTables tables = make_propagation_tables(config, grids, basis, potential.asymptotic);
  premultiply_with_time_step(potential.values, p.time_step);
  const ScatteringShape shape{basis.size(), p.nx, p.ny, p.nr, p.nz};
  const ScatteringKineticOperator kinetic(shape, tables, fft);
  const PotentialOperator potential_op(shape, grids.nz_potential, basis.initial_index(), tables, transform,
                                       potential.values);
  const SpecularOperator specular(p.nr, grids.z_specular.size(), grids.nz_potential, basis.initial_index(),
                                  *grids.z_specular_optical, grids.r_optical, tables, fft);
  const SplitOperatorPropagator propagator(config.derived.steps_per_analysis, kinetic, potential_op, specular);

  WavePacket psi{ScatteringWave(shape.n_rot, p.nx * p.ny, p.nr * p.nz),
                 make_initial_wave_packet(config, grids, gas, fft).specular};
  for (auto& c : psi.scattering.values) c = 1e-3;  // arbitrary non-zero data

  std::printf("%s\n", numerical_libraries_summary().c_str());
  std::printf("%d OpenMP threads, %d rotational states, nx = ny = %d, nr = %d, nz = %d (%.2f GB wave function)\n",
              omp_get_max_threads(), shape.n_rot, p.nx, p.nr, p.nz, double(psi.scattering.values.size()) * 16.0 / 1e9);
  time_call("T_Z + T_r (2D FFT)", repetitions, [&] { kinetic.apply_zr(psi.scattering, TimeStep::Full); });
  time_call("T_xy (2D FFT)", repetitions, [&] { kinetic.apply_xy(psi.scattering, TimeStep::Full); });
  time_call("potential incl. transforms", repetitions, [&] { potential_op.apply(psi); });
  time_call("specular grid", repetitions, [&] { specular.apply(psi.specular, TimeStep::Full, false); });
  const std::string label = "analysis step (" + std::to_string(config.derived.steps_per_analysis) + " steps)";
  time_call(label.c_str(), repetitions, [&] { propagator.propagate_analysis_step(psi); });
  return 0;
}
