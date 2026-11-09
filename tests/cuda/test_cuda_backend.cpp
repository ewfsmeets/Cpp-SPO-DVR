// Tests: the GPU backend against the CPU backend. Both start from the same random wave packet and
// propagate one analysis step; the wave packets must agree to rounding (cuFFT and the
// CPU FFT library differ in the last bits). Also: copies between host and GPU are exact, the BK
// plane and the densities equal the CPU's, and two GPU runs are bitwise identical.
#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

#include "cuda/cuda_backend.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/kernels.hpp"
#include "cuda_test_setup.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/averages.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/memory_estimate.hpp"

using namespace spodvr;
using namespace spodvr::testing;

namespace {

using cplx = std::complex<double>;

// Everything a backend is built from, for a reference case (optionally with changed keys).
struct BackendSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  FftwFft fft;
  RotationalBasis basis;
  PotentialGrid potential;
  GasPhase gas;
  PropagationTables tables;
  SpecularWave initial;

  static Config read(const std::string& case_name, const nlohmann::json& changes) {
    std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + case_name + "/input.json");
    nlohmann::json json = nlohmann::json::parse(file);
    json.update(changes);
    return build_config(parse_input(json));
  }

  explicit BackendSetup(const std::string& case_name, const nlohmann::json& changes = nlohmann::json::object())
      : config(read(case_name, changes)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        basis(angular, config.input.j_in, config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input), config, grids, angular,
                                         make_grid_symmetry(config, angular))),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())),
        tables(make_propagation_tables(config, grids, basis, potential.asymptotic)),
        initial(make_initial_wave_packet(config, grids, gas, fft).specular) {}

  [[nodiscard]] BackendInput input() const { return {config, grids, angular, basis, tables, fft}; }
  [[nodiscard]] std::unique_ptr<PropagationBackend> cpu() const {
    return make_cpu_backend(input(), potential, initial);
  }
  [[nodiscard]] std::unique_ptr<PropagationBackend> gpu() const {
    return cuda::make_cuda_backend(input(), potential, initial);
  }
};

// Loads the same pseudo-random wave packet into a backend.
void load_random(PropagationBackend& backend, std::uint64_t seed = 1) {
  backend.read_wave_packet([&](WavePacket& psi) {
    fill_pseudo_random(psi.scattering.values, seed);
    fill_pseudo_random(psi.specular.values, seed + 1);
  });
}

// One analysis step on the CPU and on the GPU from the same random wave packet.
void expect_same_step(const BackendSetup& s, double tolerance) {
  const auto cpu = s.cpu();
  const auto gpu = s.gpu();
  load_random(*cpu);
  load_random(*gpu);
  cpu->propagate_analysis_step();
  gpu->propagate_analysis_step();
  const WavePacket& a = cpu->wave_packet_on_host();
  const WavePacket& b = gpu->wave_packet_on_host();
  EXPECT_LT(relative_difference(a.scattering.values, b.scattering.values), tolerance);
  EXPECT_LT(relative_difference(a.specular.values, b.specular.values), tolerance);
  EXPECT_LT(relative_difference(cpu->specular().values, gpu->specular().values), tolerance);
}

// The GPU kinetic stage (T_Z + T_r and T_xy merged into one stage: methodology exception in
// docs/cuda.md) against the CPU's two stages, apply_zr then apply_xy, on a random scattering wave.
void expect_same_kinetic_stage(const BackendSetup& s, double tolerance) {
  const InputParameters& p = s.config.input;
  const ScatteringShape shape{s.basis.size(), p.nx, p.ny, p.nr, p.nz};
  const int n_xy = p.nx * p.ny;
  const int n_zr = p.nr * p.nz;
  const ScatteringKineticOperator cpu(shape, s.tables, s.fft);
  cuda::DeviceKineticOperator gpu(shape, s.tables, std::size_t(512) << 20);
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan* plan : gpu.fft_plans()) work_bytes = std::max(work_bytes, plan->work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan* plan : gpu.fft_plans()) plan->set_work_area(work.data());
  const cuda::DeviceLayout layout{n_zr, shape.n_rot, n_xy};
  cuda::DeviceWavePacket device(layout, p.nr, static_cast<int>(s.grids.z_specular.size()));

  for (const TimeStep step : {TimeStep::Full, TimeStep::Half}) {
    ScatteringWave psi(shape.n_rot, n_xy, n_zr);
    fill_pseudo_random(psi.values, 7);
    // Host layout -> device layout through a staging copy on the GPU.
    cuda::DeviceBuffer<cuda::Complex> staging(cuda::as_device_complex(std::span<const cplx>(psi.values)));
    cuda::host_chunk_to_device(staging.data(), layout, 0, n_zr, device.scattering.data());

    cpu.apply_zr(psi, step);
    cpu.apply_xy(psi, step);
    gpu.apply(device, step);

    cuda::device_to_host_chunk(device.scattering.data(), layout, 0, n_zr, staging.data());
    std::vector<cplx> result(psi.values.size());
    staging.download(cuda::as_device_complex(std::span<cplx>(result)));
    EXPECT_LT(relative_difference(psi.values, result), tolerance) << (step == TimeStep::Full ? "full" : "half");
  }
}

}  // namespace

// exp(i A) applied by the load callback of the inverse phi FFT gives the result of the kernel
// multiply_by_exp_potential: one potential stage on a random wave packet, with surface symmetry
// (representatives, phi maps) and several slabs (the callback's izr_first changes per slab). The
// multiplication is the same (explicitly rounded), but for a plan with a callback cuFFT may use a
// differently rounding FFT kernel: here (nphi = 12) 62% of the values differ by <= 4e-16 relative;
// for the DZ phi layout (nphi = 60) the results were bitwise equal.
TEST(CudaBackend, ExpPotentialCallbackMatchesTheKernel) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short_sym");
  const InputParameters& p = s.config.input;
  const ScatteringShape shape{s.basis.size(), p.nx, p.ny, p.nr, p.nz};
  const int n_xy = p.nx * p.ny;
  const int n_zr = p.nr * p.nz;
  auto potential = s.potential.values;  // copy, premultiplied below
  premultiply_with_time_step(potential, p.time_step);
  const AngularTransform transform(s.angular, s.basis, s.fft);
  const std::size_t scratch_limit = std::size_t(64) << 20;  // several slabs
  const std::size_t fft_limit = std::size_t(512) << 20;
  auto make = [&](cuda::ExpPotentialMode mode) {
    return cuda::DevicePotentialOperator(shape, s.grids.nz_potential, s.basis.initial_index(), s.tables, transform,
                                         potential, s.potential.symmetry, fft_limit, scratch_limit, mode);
  };
  cuda::DevicePotentialOperator with_callback = make(cuda::ExpPotentialMode::FftCallback);
  cuda::DevicePotentialOperator with_kernel = make(cuda::ExpPotentialMode::Kernel);
  ASSERT_TRUE(with_callback.exp_potential_in_fft_callback());
  ASSERT_FALSE(with_kernel.exp_potential_in_fft_callback());
  ASSERT_LT(with_callback.slab_points(), std::int64_t(s.grids.nz_potential) * p.nr);  // more than one slab

  std::vector<cuda::DeviceFftPlan*> plans = with_callback.fft_plans();
  for (cuda::DeviceFftPlan* plan : with_kernel.fft_plans()) plans.push_back(plan);
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan* plan : plans) work_bytes = std::max(work_bytes, plan->work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan* plan : plans) plan->set_work_area(work.data());

  const cuda::DeviceLayout layout{n_zr, shape.n_rot, n_xy};
  const int nz_specular = static_cast<int>(s.grids.z_specular.size());
  std::vector<cplx> scattering(layout.size());
  std::vector<cplx> specular(std::size_t(p.nr) * nz_specular);
  fill_pseudo_random(scattering, 3);
  fill_pseudo_random(specular, 4);
  auto run = [&](const cuda::DevicePotentialOperator& op) {
    cuda::DeviceWavePacket psi(layout, p.nr, nz_specular);
    psi.scattering.upload(cuda::as_device_complex(std::span<const cplx>(scattering)));
    psi.specular.upload(cuda::as_device_complex(std::span<const cplx>(specular)));
    op.apply(psi);
    std::vector<cplx> out(scattering.size() + specular.size());
    psi.scattering.download(cuda::as_device_complex(std::span<cplx>(out).first(scattering.size())));
    psi.specular.download(cuda::as_device_complex(std::span<cplx>(out).subspan(scattering.size())));
    return out;
  };
  const std::vector<cplx> a = run(with_callback);
  const std::vector<cplx> b = run(with_kernel);
  std::size_t differing = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::memcmp(&a[i], &b[i], sizeof(cplx)) != 0) ++differing;
  EXPECT_LT(relative_difference(a, b), 1e-15) << differing << " of " << a.size() << " values differ";
  EXPECT_GT(relative_difference(scattering, a), 1e-3);  // the stage did change the wave function
}

// The merged GPU kinetic stage equals the CPU's T_Z + T_r stage followed by its T_xy stage, also
// with the momentum shift (parallel momentum) and on a corrugated case's grid.
TEST(CudaBackend, MergedKineticStageMatchesTheCpuStages) {
  SPODVR_REQUIRE_GPU();
  expect_same_kinetic_stage(BackendSetup("hardwall_short"), 1e-13);
  expect_same_kinetic_stage(BackendSetup("hardwall_short", {{"EParallel", 0.05}, {"PhiParallel", 30.0}}), 1e-13);
  expect_same_kinetic_stage(BackendSetup("cu111_dailight_short"), 1e-13);
}

// Upload, download and the BK plane are exact copies; the densities equal the CPU's up to an ulp.
TEST(CudaBackend, CopiesAreExactAndDensitiesMatchTheCpu) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("hardwall_short");
  const auto cpu = s.cpu();
  const auto gpu = s.gpu();
  load_random(*cpu);
  load_random(*gpu);
  const WavePacket& a = cpu->wave_packet_on_host();
  const WavePacket& b = gpu->wave_packet_on_host();
  EXPECT_EQ(relative_difference(a.scattering.values, b.scattering.values), 0.0);
  EXPECT_EQ(relative_difference(a.specular.values, b.specular.values), 0.0);
  for (const int iz : {0, 3, s.config.input.nz - 1}) {
    const auto pa = cpu->z_plane(iz);
    const auto pb = gpu->z_plane(iz);
    ASSERT_EQ(pa.size(), pb.size());
    EXPECT_EQ(relative_difference(pa, pb), 0.0) << "iz = " << iz;
  }
  const std::vector<double> da = cpu->block_densities();
  const std::vector<double> db = gpu->block_densities();
  ASSERT_EQ(da.size(), db.size());
  for (std::size_t i = 0; i < da.size(); ++i) EXPECT_NEAR(db[i], da[i], 1e-14 * da[i]) << "izr = " << i;
}

// The hard wall (no coupling between theta, phi, x, y): every stage on a random wave packet.
TEST(CudaBackend, AnalysisStepMatchesTheCpuHardWall) {
  SPODVR_REQUIRE_GPU();
  expect_same_step(BackendSetup("hardwall_short"), 1e-13);
}

// Parallel momentum: the momentum shift around the (x, y) transform.
TEST(CudaBackend, AnalysisStepMatchesTheCpuWithParallelMomentum) {
  SPODVR_REQUIRE_GPU();
  expect_same_step(BackendSetup("hardwall_short", {{"EParallel", 0.05}, {"PhiParallel", 30.0}}), 1e-13);
}

// Corrugated PES (Dai-Light, 60 deg cell), without and with C6v surface symmetry in the potential.
TEST(CudaBackend, AnalysisStepMatchesTheCpuDaiLight) {
  SPODVR_REQUIRE_GPU();
  expect_same_step(BackendSetup("cu111_dailight_short"), 1e-13);
}

TEST(CudaBackend, AnalysisStepMatchesTheCpuDaiLightWithSymmetry) {
  SPODVR_REQUIRE_GPU();
  expect_same_step(BackendSetup("cu111_dailight_short_sym"), 1e-13);
}

// The GPU results do not change from run to run (no atomics, fixed launch shapes).
TEST(CudaBackend, TwoRunsAreBitwiseIdentical) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short_sym");
  const auto a = s.gpu();
  const auto b = s.gpu();
  load_random(*a);
  load_random(*b);
  for (int i = 0; i < 2; ++i) {
    a->propagate_analysis_step();
    b->propagate_analysis_step();
  }
  EXPECT_EQ(relative_difference(a->wave_packet_on_host().scattering.values, b->wave_packet_on_host().scattering.values),
            0.0);
}

// The GPU memory estimate (printed before the potential is computed, and checked against the free
// GPU memory) matches what the backend allocates: at least the exactly known arrays, at most the
// estimate (which counts the cuFFT work area at its limit and the copy buffer).
TEST(CudaBackend, MemoryEstimateMatchesTheAllocations) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short_sym");
  const MemoryEstimate e =
      cuda::estimate_gpu_memory(s.config, s.grids, s.angular, s.potential.symmetry.n_representatives());
  std::size_t exact = 0;
  for (const MemoryItem& item : e.items)
    if (item.name != "cuFFT work area" && item.name != "Copy buffer") exact += item.bytes;
  const std::size_t before = cuda::free_gpu_memory();
  const auto gpu = s.gpu();
  const std::size_t used = before - cuda::free_gpu_memory();
  const std::size_t slack = std::size_t(64) << 20;  // allocation granularity, cuFFT plans
  EXPECT_GE(used + slack, exact);
  EXPECT_LE(used, e.peak_bytes + slack);
}
