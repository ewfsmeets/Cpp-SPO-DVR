// Per-operator and per-kernel timing of the GPU backend for an input file: the specular grid,
// the kinetic stage T, a whole analysis step (as the program runs it), the potential stage kernel
// by kernel, the per-step copies to the host, and the flux analysis: the flux surface (kernel and
// chunked copies) and, with doFluxAnalysis, one flush of the accumulator on the host per slot count
// and the flux cost as a share of the analysis step (with its BK plane).
// Every call is followed by cudaDeviceSynchronize, so the times are those of the GPU work.
// Usage: spodvr_cuda_operator_timing <input.json> [repetitions]   (host threads: OMP_NUM_THREADS)
#include <cuda_runtime_api.h>

#include <chrono>
#include <complex>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda/device_flux_surface.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/kernels.hpp"
#include "fft/fftw_fft.hpp"
#include "flux_timing.hpp"
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
      compute_potential_grid(*make_potential(p), config, grids, angular, make_grid_symmetry(config, angular));
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
                                  grids.r_optical, tables, fft_work_limit);
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
  step_ms += time_call("T (Z, r, x, y; one stage)", repetitions, [&] { kinetic.apply(psi, TimeStep::Full); });
  time_call("  of which the two (Z, r) cuFFTs", repetitions, [&] {
    kinetic_plans[0]->execute(psi.scattering.data());
    kinetic_plans[1]->execute(psi.scattering.data());
  });
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
    kinetic.apply(psi, step);
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
  // With the FFT callback, exp(i A) is applied inside the inverse phi FFT; the kernel is timed for
  // comparison only (the callback reads the last slab's izr_first here: the timing is the same).
  const bool in_callback = potential_op.exp_potential_in_fft_callback();
  time_call(in_callback ? "multiply_by_exp_potential (not used)" : "multiply_by_exp_potential", repetitions,
            over_slabs([&](const PotentialArgs& a) { multiply_by_exp_potential(a, dvr); }));
  if (!phi.empty())
    time_call(in_callback ? "phi cuFFT (to FBR) with exp(iA) callback" : "phi cuFFT (to FBR)", repetitions,
              over_slabs([&](const PotentialArgs&) { phi[1]->execute(dvr); }));
  time_call("transform_to_fbr_and_rotate", repetitions,
            over_slabs([&](const PotentialArgs& a) { transform_to_fbr_and_rotate(a, dvr, psi.scattering.data()); }));
  time_call("project_onto_specular", repetitions, over_slabs([&](const PotentialArgs& a) {
              project_onto_specular(a, psi.scattering.data(), psi.specular.data());
            }));

  std::printf("Per analysis step (host side):\n");
  // The backend's read-back buffer: the BK plane, or a chunk of nr / 2 planes of the flux surface.
  DeviceBuffer<Complex> readback(std::size_t(layout.n_rot * layout.n_xy * p.nr));
  std::vector<std::complex<double>> host_plane(readback.size());
  const double bk_plane_ms = time_call("BK plane: pack + copy to host", repetitions, [&] {
    pack_z_plane(psi.scattering.data(), layout, p.nr, grids.iz_analysis, readback.data());
    readback.download(as_device_complex(std::span<std::complex<double>>(host_plane)));
  });
  DeviceBuffer<double> densities(std::size_t(layout.n_zr));
  // Only with the averages (every dtAve), so not part of the analysis step below.
  time_call("densities (every dtAve)", repetitions,
            [&] { block_densities(psi.scattering.data(), layout, densities.data()); });

  // The flux surface at rFluxAnalysis, or in the middle of the r grid without flux analysis (the
  // cost does not depend on r_c): the whole read, then its kernels and its copies to the host alone.
  const int ir_flux = grids.ir_flux_analysis.value_or(p.nr / 2);
  const DeviceBuffer<double> flux_weights(
      std::span<const double>(flux_derivative_weights(p.nr, grids.r.spacing(), ir_flux)));
  const int n_blocks = static_cast<int>(layout.n_rot * layout.n_xy);
  const int n_planes = grids.nz_potential;
  const int chunk_planes = flux_surface_chunk_planes(readback.size(), n_blocks);
  FluxSurfaceWave slot(n_blocks, n_planes);
  std::printf(
      "Flux surface (per analysis step): %d planes of %.2f MB each for psi and d psi / dr, %d chunks of <= %d "
      "planes\n",
      n_planes, 16.0 * n_blocks / 1e6, (n_planes + chunk_planes - 1) / chunk_planes, chunk_planes);
  const auto fill = [&](FluxSurfaceWave& s) { read_flux_surface(psi, ir_flux, flux_weights, readback, s); };
  const double surface_ms = time_call("flux surface: kernels + copies to host", repetitions, [&] { fill(slot); });
  auto over_chunks = [&](const std::function<void(int first, int count)>& part) {
    return [&, part] {
      for (int first = 0; first < n_planes; first += chunk_planes)
        part(first, std::min(chunk_planes, n_planes - first));
    };
  };
  time_call("  of which the kernels", repetitions, over_chunks([&](int first, int count) {
              flux_surface_chunk(psi.scattering.data(), layout, p.nr, ir_flux, flux_weights.data(), first, count,
                                 readback.data());
            }));
  time_call("  of which the copies to host (cudaMemcpy)", repetitions, over_chunks([&](int first, int count) {
              const std::size_t values = std::size_t(n_blocks) * count;
              const auto into = [&](ComplexVector& v) {
                return as_device_complex(std::span<std::complex<double>>(v).subspan(slot.index(0, first), values));
              };
              readback.download(into(slot.psi));
              readback.download(into(slot.derivative), values);
            }));

  // The flush of the accumulator (host, OpenMP) per slot count. The slots hold the surface read
  // above (filling them is the row above); the energies and c_n are the input's.
  if (!p.do_flux_analysis) {
    std::printf("Flux flush: skipped (doFluxAnalysis is off)\n");
    return 0;
  }
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());
  const FluxEnergies energies = make_flux_energies(config, tables.parallel_energy, gas.initial.energy);
  const std::vector<benchmark::FlushTiming> flushes =
      benchmark::time_flux_flushes(n_blocks, n_planes, energies, config.derived.analysis_time_step, repetitions, fill);
  benchmark::print_flux_cost(energies.size(), surface_ms, flushes, analysis_ms + bk_plane_ms,
                             2 * slot.psi.size() * sizeof(std::complex<double>));
  return 0;
}
