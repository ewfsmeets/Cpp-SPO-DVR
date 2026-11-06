// Per-operator and per-kernel timing of the GPU backend for an input file: the specular grid,
// T_Z + T_r, T_xy, a whole analysis step (as the program runs it), the potential stage kernel
// by kernel, and the per-step copies to the host.
// Every call is followed by cudaDeviceSynchronize, so the times are those of the GPU work.
// Usage: spodvr_cuda_operator_timing <input.json> [repetitions]
#include <cuda_runtime_api.h>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "cuda/cuda_backend.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/kernels.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"

using namespace spodvr;
using namespace spodvr::cuda;

namespace {

// Runs `step` once (warm-up), then `repetitions` times; prints the mean time per call.
double time_call(const std::string& name, int repetitions, const std::function<void()>& step) {
  step();
  cudaDeviceSynchronize();
  const auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < repetitions; ++i) step();
  cudaDeviceSynchronize();
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  const double ms = elapsed.count() / repetitions * 1000.0;
  std::printf("  %-44s %10.2f ms per call\n", name.c_str(), ms);
  return ms;
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
  initialise_cuda();
  std::printf("%s\n", cuda_device_summary().c_str());

  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  PotentialGrid potential =
      compute_potential_grid(*make_potential(p.system_id), config, grids, angular, make_grid_symmetry(config, angular));
  const FftwFft fft;
  const RotationalBasis basis(angular, p.j_in, p.mj_in);
  const PropagationTables tables = make_propagation_tables(config, grids, basis, potential.asymptotic);
  premultiply_with_time_step(potential.values, p.time_step);
  const ScatteringShape shape{basis.size(), p.nx, p.ny, p.nr, p.nz};
  const std::size_t fft_work_limit = std::size_t(512) << 20;

  const DeviceLayout layout{std::int64_t(p.nr) * p.nz, shape.n_rot, std::int64_t(p.nx) * p.ny};
  DeviceWavePacket psi(layout, p.nr, grids.z_specular.size());
  {
    std::vector<std::complex<double>> host(layout.size(), {1e-3, 0.0});  // arbitrary non-zero data
    psi.scattering.upload(as_device_complex(std::span<const std::complex<double>>(host)));
    psi.specular.set_zero();
  }
  DeviceKineticOperator kinetic(shape, tables, fft_work_limit);
  DeviceSpecularOperator specular(p.nr, grids.z_specular.size(), grids.nz_potential, basis.initial_index(),
                                  *grids.z_specular_optical, grids.r_optical, tables, fft_work_limit);
  const AngularTransform transform(angular, basis, fft);
  DevicePotentialOperator potential_op(shape, grids.nz_potential, basis.initial_index(), tables, transform,
                                       potential.values, potential.symmetry, fft_work_limit);
  const std::vector<DeviceFftPlan*> kinetic_plans = kinetic.fft_plans();  // (Z, r) pair, then (x, y) pair if any
  std::vector<DeviceFftPlan*> plans = kinetic_plans;
  for (DeviceFftPlan* q : specular.fft_plans()) plans.push_back(q);
  for (DeviceFftPlan* q : potential_op.fft_plans()) plans.push_back(q);
  std::size_t work_bytes = 1;
  for (const DeviceFftPlan* q : plans) work_bytes = std::max(work_bytes, q->work_size());
  DeviceBuffer<std::byte> work(work_bytes);
  for (DeviceFftPlan* q : plans) q->set_work_area(work.data());

  std::printf(
      "%d rotational states, nx = %d, ny = %d, nr = %d, nz = %d (%.2f GB wave function); %lld (Z, r) "
      "points per potential slab; cuFFT work area %.0f MiB\n",
      shape.n_rot, p.nx, p.ny, p.nr, p.nz, double(layout.size()) * 16.0 / 1e9,
      static_cast<long long>(potential_op.slab_points()), double(work_bytes) / double(1 << 20));

  std::printf("Stages (one full time step):\n");
  double step_ms = 0.0;
  step_ms += time_call("specular grid", repetitions, [&] { specular.apply(psi, TimeStep::Full, false); });
  step_ms += time_call("T_Z + T_r", repetitions, [&] { kinetic.apply_zr(psi, TimeStep::Full); });
  time_call("  of which the two (Z, r) cuFFTs", repetitions, [&] {
    kinetic_plans[0]->execute(psi.scattering.data());
    kinetic_plans[1]->execute(psi.scattering.data());
  });
  step_ms += time_call("T_xy", repetitions, [&] { kinetic.apply_xy(psi, TimeStep::Full); });
  if (kinetic_plans.size() == 4)
    time_call("  of which the two (x, y) cuFFTs", repetitions, [&] {
      kinetic_plans[2]->execute(psi.scattering.data());
      kinetic_plans[3]->execute(psi.scattering.data());
    });
  step_ms += time_call("potential (all kernels)", repetitions, [&] { potential_op.apply(psi); });
  std::printf("  %-44s %10.2f ms\n", "sum = one time step", step_ms);

  // One analysis step as the program runs it: T/2 V (T V)^(n-1) T/2, i.e. one kinetic stage more
  // than n time steps. Same order as CudaBackend::propagate_analysis_step (cuda_backend.cpp).
  const int n_per_analysis = config.derived.steps_per_analysis;
  auto kinetic_and_potential = [&](TimeStep step, bool kinetic_only) {
    specular.apply(psi, step, kinetic_only);
    kinetic.apply_zr(psi, step);
    kinetic.apply_xy(psi, step);
    if (!kinetic_only) potential_op.apply(psi);
  };
  const double analysis_ms =
      time_call("one analysis step (" + std::to_string(n_per_analysis) + " time steps)", repetitions, [&] {
        kinetic_and_potential(TimeStep::Half, false);
        for (int s = 1; s < n_per_analysis; ++s) kinetic_and_potential(TimeStep::Full, false);
        kinetic_and_potential(TimeStep::Half, true);
      });
  std::printf("  %-44s %10.2f ms\n", "  per time step", analysis_ms / n_per_analysis);

  // The potential stage kernel by kernel, summed over the slabs.
  std::printf("Potential stage, kernel by kernel (summed over the slabs):\n");
  const PotentialArgs base = potential_op.kernel_args();
  Complex* dvr = potential_op.scratch();
  const std::vector<DeviceFftPlan*> phi = potential_op.fft_plans();
  auto over_slabs = [&](const std::function<void(const PotentialArgs&)>& kernel) {
    return [&, kernel] {
      PotentialArgs slab = base;
      for (slab.izr_first = 0; slab.izr_first < base.n_zr_potential; slab.izr_first += base.n_slab) kernel(slab);
    };
  };
  time_call("apply_gas_phase_step", repetitions, [&] { apply_gas_phase_step(base, psi.scattering.data()); });
  time_call("absorb_and_transform_to_dvr", repetitions, over_slabs([&](const PotentialArgs& a) {
              absorb_and_transform_to_dvr(a, psi.scattering.data(), psi.specular.data(), dvr);
            }));
  if (!phi.empty())
    time_call("phi cuFFT (to DVR)", repetitions, over_slabs([&](const PotentialArgs&) { phi[0]->execute(dvr); }));
  time_call("multiply_by_exp_potential", repetitions,
            over_slabs([&](const PotentialArgs& a) { multiply_by_exp_potential(a, dvr); }));
  if (!phi.empty())
    time_call("phi cuFFT (to FBR)", repetitions, over_slabs([&](const PotentialArgs&) { phi[1]->execute(dvr); }));
  time_call("transform_to_fbr_and_rotate", repetitions,
            over_slabs([&](const PotentialArgs& a) { transform_to_fbr_and_rotate(a, dvr, psi.scattering.data()); }));
  time_call("project_onto_specular", repetitions, over_slabs([&](const PotentialArgs& a) {
              project_onto_specular(a, psi.scattering.data(), psi.specular.data());
            }));

  std::printf("Per analysis step (host side):\n");
  DeviceBuffer<Complex> plane(std::size_t(layout.n_rot * layout.n_xy * p.nr));
  std::vector<std::complex<double>> host_plane(plane.size());
  time_call("BK plane: pack + copy to host", repetitions, [&] {
    pack_z_plane(psi.scattering.data(), layout, p.nr, grids.iz_analysis, plane.data());
    plane.download(as_device_complex(std::span<std::complex<double>>(host_plane)));
  });
  DeviceBuffer<double> densities(std::size_t(layout.n_zr));
  time_call("densities", repetitions, [&] { block_densities(psi.scattering.data(), layout, densities.data()); });
  return 0;
}
