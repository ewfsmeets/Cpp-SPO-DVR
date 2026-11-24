// Tests: the GPU backend against the CPU backend. Both start from the same random wave packet and
// propagate one analysis step; the wave packets must agree to rounding (cuFFT and the
// CPU FFT library differ in the last bits). Also: copies between host and GPU are exact, the BK
// plane and the densities equal the CPU's, and two GPU runs are bitwise identical.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

#include "cuda/cuda_backend.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/kernels.hpp"
#include "cuda_test_setup.hpp"
#include "designer_cases.hpp"
#include "fft/fftw_fft.hpp"
#include "frames/frame_source.hpp"
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
#include "propagation/potential_operator.hpp"

using namespace spodvr;
using namespace spodvr::testing;

namespace {

using cplx = std::complex<double>;

// Everything a backend is built from, for a reference case (optionally with changed keys).
struct BackendSetup {
  std::unique_ptr<Potential> pes;  // built once, as the program does
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  FftwFft fft;
  RotationalBasis basis;
  PotentialGrid potential;
  GasPhase gas;
  PropagationTables tables;
  SpecularWave initial;

  static InputParameters read(const std::string& case_name, const nlohmann::json& changes) {
    std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + case_name + "/input.json");
    nlohmann::json json = nlohmann::json::parse(file);
    json.update(changes);
    return parse_input(json);
  }

  explicit BackendSetup(const std::string& case_name, const nlohmann::json& changes = nlohmann::json::object())
      : pes(make_potential(read(case_name, changes))),
        config(build_config(read(case_name, changes), *pes)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        basis(angular, config.input.j_in, config.input.mj_in),
        potential(compute_potential_grid(*pes, config, grids, angular, make_grid_symmetry(config, angular))),
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
  for (const cuda::DeviceFftPlan& plan : gpu.fft_plans()) work_bytes = std::max(work_bytes, plan.work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan& plan : gpu.fft_plans()) plan.set_work_area(work.data());
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

// The CPU-GPU difference of these tests is NaN when a value on either side is NaN, so the
// EXPECT_LT and EXPECT_EQ on it fail (std::max would drop the NaN difference). Host only.
TEST(CudaTestSetup, RelativeDifferenceIsNaNWithANaNValue) {
  const std::vector<cplx> a = {{1.0, 0.5}, {-0.25, 2.0}};
  std::vector<cplx> b = a;
  EXPECT_EQ(relative_difference(a, b), 0.0);
  b[1] = {std::numeric_limits<double>::quiet_NaN(), 2.0};
  EXPECT_TRUE(std::isnan(relative_difference(a, b)));
  EXPECT_TRUE(std::isnan(relative_difference(b, a)));
}

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

  cuda::DeviceFftPlans plans = with_callback.fft_plans();
  for (cuda::DeviceFftPlan& plan : with_kernel.fft_plans()) plans.emplace_back(plan);
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan& plan : plans) work_bytes = std::max(work_bytes, plan.work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan& plan : plans) plan.set_work_area(work.data());

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

// The GPU potential stage executed in several slabs of (Z, r) points equals the CPU's
// PotentialOperator (all points at once): one stage on a random wave packet, with surface symmetry,
// the coupling to the specular grid and the exp(i A) callback (the default mode).
TEST(CudaPotential, SlabsEqualTheCpu) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short_sym");
  const InputParameters& p = s.config.input;
  const ScatteringShape shape{s.basis.size(), p.nx, p.ny, p.nr, p.nz};
  const int n_xy = p.nx * p.ny;
  const int n_zr = p.nr * p.nz;
  const int nz_specular = static_cast<int>(s.grids.z_specular.size());
  auto potential = s.potential.values;  // copy, premultiplied below
  premultiply_with_time_step(potential, p.time_step);
  const AngularTransform transform(s.angular, s.basis, s.fft);

  const std::size_t scratch_limit = std::size_t(64) << 20;  // several slabs
  const std::size_t fft_limit = std::size_t(512) << 20;
  cuda::DevicePotentialOperator gpu(shape, s.grids.nz_potential, s.basis.initial_index(), s.tables, transform,
                                    potential, s.potential.symmetry, fft_limit, scratch_limit);
  const std::int64_t interaction_points = std::int64_t(s.grids.nz_potential) * p.nr;
  ASSERT_LT(gpu.slab_points(), interaction_points);  // more than one slab
  ASSERT_TRUE(gpu.exp_potential_in_fft_callback());
  const PotentialOperator cpu(shape, s.grids.nz_potential, s.basis.initial_index(), s.tables, transform, potential,
                              s.potential.symmetry);

  WavePacket host{ScatteringWave(shape.n_rot, n_xy, n_zr), SpecularWave(p.nr, nz_specular)};
  fill_pseudo_random(host.scattering.values, 5);
  fill_pseudo_random(host.specular.values, 6);
  const WavePacket before = host;

  // GPU: host layout -> device layout through a staging copy, one stage, and back.
  cuda::DeviceFftPlans plans = gpu.fft_plans();
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan& plan : plans) work_bytes = std::max(work_bytes, plan.work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan& plan : plans) plan.set_work_area(work.data());
  const cuda::DeviceLayout layout{n_zr, shape.n_rot, n_xy};
  cuda::DeviceWavePacket device(layout, p.nr, nz_specular);
  cuda::DeviceBuffer<cuda::Complex> staging(cuda::as_device_complex(std::span<const cplx>(host.scattering.values)));
  cuda::host_chunk_to_device(staging.data(), layout, 0, n_zr, device.scattering.data());
  device.specular.upload(cuda::as_device_complex(std::span<const cplx>(host.specular.values)));
  gpu.apply(device);
  cuda::device_to_host_chunk(device.scattering.data(), layout, 0, n_zr, staging.data());
  std::vector<cplx> gpu_scattering(host.scattering.values.size());
  std::vector<cplx> gpu_specular(host.specular.values.size());
  staging.download(cuda::as_device_complex(std::span<cplx>(gpu_scattering)));
  device.specular.download(cuda::as_device_complex(std::span<cplx>(gpu_specular)));

  cpu.apply(host);
  const double scattering_difference = relative_difference(host.scattering.values, gpu_scattering);
  const double specular_difference = relative_difference(host.specular.values, gpu_specular);
  EXPECT_LT(scattering_difference, 1e-13) << gpu.slab_points() << " points per slab of " << interaction_points;
  EXPECT_LT(specular_difference, 1e-13);
  // The stage changed both parts (the specular one through its coupling to the scattering grid).
  EXPECT_GT(relative_difference(before.scattering.values, host.scattering.values), 1e-3);
  EXPECT_GT(relative_difference(before.specular.values, host.specular.values), 1e-3);
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
  // site_column (the energy-resolved sums): the same bytes, at the first, an inner and the last grid point
  const std::size_t column = a.scattering.values.size() / std::size_t(s.config.input.nx * s.config.input.ny);
  for (const int ixy : {0, 1, s.config.input.nx * s.config.input.ny - 1}) {
    std::vector<std::complex<double>> ca(column), cb(column);
    cpu->site_column(ixy, ca);
    gpu->site_column(ixy, cb);
    EXPECT_EQ(relative_difference(ca, cb), 0.0) << "ixy = " << ixy;
  }
}

// densities/ (writeDensities) of the same wave packet: the specular cells are bitwise equal (an exact copy),
// the scattering cells come from block_densities() (CUDA hypot) and agree to 1e-14 each.
TEST(CudaBackend, DensityRecordsOfTheSameWaveMatchTheCpu) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short");
  const auto cpu = s.cpu();
  const auto gpu = s.gpu();
  load_random(*cpu);
  load_random(*gpu);
  const CellProbabilities a = cell_probabilities(cpu->block_densities(), cpu->specular(), s.config);
  const CellProbabilities b = cell_probabilities(gpu->block_densities(), gpu->specular(), s.config);
  EXPECT_EQ(a.specular, b.specular);
  ASSERT_EQ(a.scattering.size(), b.scattering.size());
  for (std::size_t i = 0; i < a.scattering.size(); ++i)
    EXPECT_NEAR(b.scattering[i], a.scattering[i], 1e-14 * a.scattering[i]) << "izr = " << i;
}

// frames/ (dtFrames) of the same wave packet, with every FrameExtras extra: the reductions run on the host on exact
// copies of the planes and of the specular wave, so the CPU and the GPU frame are bitwise equal.
TEST(CudaBackend, FrameOfTheSameWaveIsTheCpus) {
  SPODVR_REQUIRE_GPU();
  const nlohmann::json extras = {
      {"dtFrames", 100.0},
      {"FrameExtras", nlohmann::json::array({"angular", "diffraction", "vibration", "rotation_moments", "momentum"})}};
  const BackendSetup s("cu111_dailight_short", extras);
  const auto cpu = s.cpu();
  const auto gpu = s.gpu();
  load_random(*cpu);
  load_random(*gpu);
  FrameSource cpu_frames(s.config, s.grids, s.angular, s.basis, s.gas, s.tables, s.fft);
  FrameSource gpu_frames(s.config, s.grids, s.angular, s.basis, s.gas, s.tables, s.fft);
  const FrameRecord& a = cpu_frames.frame(*cpu);
  const FrameRecord& b = gpu_frames.frame(*gpu);
  EXPECT_EQ(a.prob_zrxy, b.prob_zrxy);
  EXPECT_EQ(a.prob_tail_zr, b.prob_tail_zr);
  EXPECT_EQ(a.pop_jmj_z, b.pop_jmj_z);
  EXPECT_EQ(a.prob_angle_z, b.prob_angle_z);
  EXPECT_EQ(a.pop_nm_z, b.pop_nm_z);
  EXPECT_EQ(a.pop_vj_z, b.pop_vj_z);
  EXPECT_EQ(a.p2_numerator_z, b.p2_numerator_z);
  EXPECT_EQ(a.rot_numerator_z, b.rot_numerator_z);
  EXPECT_EQ(a.prob_kr, b.prob_kr);
  EXPECT_EQ(a.prob_kz_spec, b.prob_kz_spec);
  EXPECT_GT(*std::max_element(a.prob_zrxy.begin(), a.prob_zrxy.end()), 0.0);
  EXPECT_GT(*std::max_element(a.prob_angle_z.begin(), a.prob_angle_z.end()), 0.0);
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

namespace {

// Both low-Z absorbers of upstream 5.24 on, from ZStart to Zana (scattering) and to Z = 6
// (specular): on random data every point inside them is absorbed.
const nlohmann::json low_z_absorbers = {{"ZOptLowStart", -2.0},   {"ZOptLowEnd", 8.0},   {"ZOptLowEnergy", 0.3},
                                        {"ZSpOptLowStart", -2.0}, {"ZSpOptLowEnd", 6.0}, {"ZSpOptLowFactor", 0.01}};

// The specular step (kinetic part, gas phase, absorbers) of the CPU and the GPU twin on one
// random specular wave.
std::vector<cplx> cpu_specular_step(const BackendSetup& s, const std::vector<cplx>& initial) {
  const SpecularOperator op(s.config.input.nr, s.grids.z_specular.size(), s.grids.nz_potential, s.basis.initial_index(),
                            s.grids.r_optical, s.tables, s.fft);
  SpecularWave psi(s.config.input.nr, s.grids.z_specular.size());
  std::copy(initial.begin(), initial.end(), psi.values.begin());
  op.apply(psi, TimeStep::Full, false);
  return {psi.values.begin(), psi.values.end()};
}

// The same step on the GPU twin (DeviceSpecularOperator).
std::vector<cplx> gpu_specular_step(const BackendSetup& s, const std::vector<cplx>& initial) {
  const int nr = s.config.input.nr;
  const int nz = s.grids.z_specular.size();
  cuda::DeviceSpecularOperator op(nr, nz, s.grids.nz_potential, s.basis.initial_index(), s.grids.r_optical, s.tables,
                                  std::size_t(512) << 20);
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan& plan : op.fft_plans()) work_bytes = std::max(work_bytes, plan.work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan& plan : op.fft_plans()) plan.set_work_area(work.data());
  cuda::DeviceWavePacket psi(cuda::DeviceLayout{1, 1, 1}, nr, nz);
  psi.specular.upload(cuda::as_device_complex(std::span<const cplx>(initial)));
  op.apply(psi, TimeStep::Full, false);
  std::vector<cplx> result(initial.size());
  psi.specular.download(cuda::as_device_complex(std::span<cplx>(result)));
  return result;
}

}  // namespace

// The specular step applies the Z factors at every Z point on the GPU as on the CPU, so a low-Z
// absorber on the specular grid acts on both (the two differ only by the FFT libraries' rounding).
// Without the absorber the result differs by far more than that, so the comparison sees it.
TEST(CudaBackend, SpecularStepMatchesTheCpuWithALowZAbsorber) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup with("hardwall_short", low_z_absorbers);
  const BackendSetup without("hardwall_short");
  std::vector<cplx> initial(std::size_t(with.config.input.nr) * with.grids.z_specular.size());
  fill_pseudo_random(initial, 11);
  const std::vector<cplx> cpu = cpu_specular_step(with, initial);
  EXPECT_LT(relative_difference(cpu, gpu_specular_step(with, initial)), 1e-13);
  EXPECT_GT(relative_difference(cpu, cpu_specular_step(without, initial)), 1e-3);
}

// One analysis step of the whole backend with both low-Z absorbers on.
TEST(CudaBackend, AnalysisStepMatchesTheCpuWithLowZAbsorbers) {
  SPODVR_REQUIRE_GPU();
  expect_same_step(BackendSetup("hardwall_short", low_z_absorbers), 1e-13);
}

// The wave-packet designer copies the GPU memory estimate too (designer_cases.hpp). Only arithmetic on the
// input's sizes: no GPU needed. The program passes the representative count of the run's symmetry.
TEST(CudaBackend, MemoryEstimateMatchesTheWavePacketDesignerCases) {
  for (const nlohmann::json& designer_case : designer_cases()) {
    SCOPED_TRACE("designer case " + designer_case.at("name").get<std::string>());
    const Config config = build_config(parse_input(designer_case_input(designer_case)));
    if (config.derived.run_mode == 6) {  // a potential-only run prints no GPU estimate
      EXPECT_TRUE(designer_case.at("gpu").is_null());
      continue;
    }
    const SpatialGrids grids = make_spatial_grids(config);
    const AngularGrid angular = make_angular_grid(config);
    const int n_representatives = make_grid_symmetry(config, angular).n_representatives();
    expect_designer_estimate(cuda::estimate_gpu_memory(config, grids, angular, n_representatives),
                             designer_case.at("gpu"));
  }
}
