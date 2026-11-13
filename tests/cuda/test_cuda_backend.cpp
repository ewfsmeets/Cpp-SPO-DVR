// Tests: the GPU backend against the CPU backend. Both start from the same random wave packet and
// propagate one analysis step; the wave packets must agree to rounding (cuFFT and the
// CPU FFT library differ in the last bits). Also: copies between host and GPU are exact, the BK
// plane, the densities and the flux surface equal the CPU's, and two GPU runs are bitwise identical.
#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "backend_setup.hpp"
#include "cuda/cuda_backend.hpp"
#include "cuda/device_flux_surface.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/kernels.hpp"
#include "cuda_test_setup.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/averages.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/memory_estimate.hpp"

using namespace spodvr;
using namespace spodvr::testing;

namespace {

using cplx = std::complex<double>;

// The GPU backend of a set-up (the CPU one is BackendSetup::cpu).
std::unique_ptr<PropagationBackend> gpu_backend(const BackendSetup& s) {
  return cuda::make_cuda_backend(s.input(), s.potential, s.initial);
}

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
  const auto gpu = gpu_backend(s);
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
  const auto gpu = gpu_backend(s);
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
  const auto a = gpu_backend(s);
  const auto b = gpu_backend(s);
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
  const auto gpu = gpu_backend(s);
  const std::size_t used = before - cuda::free_gpu_memory();
  const std::size_t slack = std::size_t(64) << 20;  // allocation granularity, cuFFT plans
  EXPECT_GE(used + slack, exact);
  EXPECT_LE(used, e.peak_bytes + slack);
}

// The host memory of a GPU run in the estimate (MemoryEstimate::gpu_host) is that of the GPU backend's
// host arrays: during the set-up V (one double) and exp(i A) (one complex value) per value of the potential
// grid, its symmetry maps, two specular waves and the read-back plane of z_plane; in the time loop the
// host copy of psi (wave_packet_on_host), three specular waves and the read-back plane.
TEST(CudaBackend, HostMemoryEstimateMatchesTheAllocations) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s("cu111_dailight_short_sym");
  const MemoryEstimate e = estimate_memory(s.config, s.grids, s.angular, 1);
  const GridSymmetry& symmetry = s.potential.symmetry;
  ASSERT_LT(symmetry.n_representatives(), symmetry.n_xy());  // the maps are counted
  std::size_t maps = 2 * std::size_t(symmetry.n_xy());       // representative and operation per (x, y)
  for (int iop = 0; iop < static_cast<int>(symmetry.operations().size()); ++iop) maps += symmetry.phi_index(iop).size();

  const auto gpu = gpu_backend(s);
  const std::size_t complex_bytes = sizeof(cplx);
  const std::size_t plane = gpu->z_plane(0).size() * complex_bytes;
  const WavePacket& host = gpu->wave_packet_on_host();
  const std::size_t specular = host.specular.values.size() * complex_bytes;
  EXPECT_EQ(e.gpu_host.setup,
            s.potential.values.size() * (sizeof(double) + complex_bytes) + maps * sizeof(int) + 2 * specular + plane);
  EXPECT_EQ(e.gpu_host.time_loop, host.scattering.values.size() * complex_bytes + 3 * specular + plane);
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

std::vector<cplx> gpu_specular_step(const BackendSetup& s, const std::vector<cplx>& initial) {
  const int nr = s.config.input.nr;
  const int nz = s.grids.z_specular.size();
  cuda::DeviceSpecularOperator op(nr, nz, s.grids.nz_potential, s.basis.initial_index(), s.grids.r_optical, s.tables,
                                  std::size_t(512) << 20);
  std::size_t work_bytes = 1;
  for (const cuda::DeviceFftPlan* plan : op.fft_plans()) work_bytes = std::max(work_bytes, plan->work_size());
  cuda::DeviceBuffer<std::byte> work(work_bytes);
  for (cuda::DeviceFftPlan* plan : op.fft_plans()) plan->set_work_area(work.data());
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

// ------------------------------------------------------------------------------------------------
// The flux surface: psi and d psi / dr at r_c on the planes below Z_asymptotic (flux analysis)
// ------------------------------------------------------------------------------------------------

namespace {

// hardwall_short with the flux analysis at r_c = 2.0 bohr (ir 8 of nr = 56).
const std::string flux_case = "hardwall_short";
nlohmann::json flux_keys() { return flux_analysis_keys(flux_case, 2.0, 0.15); }

// NaN, read from the bits (also in fast-math builds, where std::isnan may be folded to false).
bool is_nan(double x) {
  const std::uint64_t bits = std::bit_cast<std::uint64_t>(x);
  return (bits & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL && (bits & 0x000fffffffffffffULL) != 0;
}

// Doubles of a and b that differ: equal bit patterns, or NaN on both sides (the CPU's and the GPU's
// NaN differ in their bits), count as equal. No floating-point comparison, so NaN on one side only
// is a difference in every build.
int value_differences(std::span<const cplx> a, std::span<const cplx> b) {
  auto differ = [](double x, double y) {
    return std::bit_cast<std::uint64_t>(x) != std::bit_cast<std::uint64_t>(y) && !(is_nan(x) && is_nan(y));
  };
  int differences = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (differ(a[i].real(), b[i].real())) ++differences;
    if (differ(a[i].imag(), b[i].imag())) ++differences;
  }
  return differences;
}

int nan_count(std::span<const cplx> values) {
  int count = 0;
  for (const cplx& v : values) count += int(is_nan(v.real())) + int(is_nan(v.imag()));
  return count;
}

// A surface whose every element is NaN: an element the backend does not write stays NaN.
FluxSurfaceWave nan_surface(int n_blocks, int n_planes) {
  FluxSurfaceWave out(n_blocks, n_planes);
  const cplx nan(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
  std::ranges::fill(out.psi, nan);
  std::ranges::fill(out.derivative, nan);
  return out;
}

#ifndef __FAST_MATH__  // the fast build compares with a tolerance, which such mistakes may stay within
// The porting mistakes of the r-derivative sum that the CPU-GPU comparison must catch.
enum class Mistake { DescendingSum, SkippedZeroTerm };

// d psi / dr at r_c of every column of the host wave psi on the planes iz < n_planes, summed with
// `mistake` (flux_r_derivative's order otherwise: from 0, ir ascending, the zero term included).
std::vector<cplx> mistaken_derivative(const ScatteringWave& psi, int ir_flux, std::span<const double> weights,
                                      int n_planes, Mistake mistake) {
  const int nr = static_cast<int>(weights.size());
  std::vector<cplx> derivative(std::size_t(psi.block_size()) * n_planes);
  for (int iz = 0; iz < n_planes; ++iz)
    for (int block = 0; block < psi.block_size(); ++block) {
      cplx sum(0.0, 0.0);
      for (int k = 0; k < nr; ++k) {
        const int ir = mistake == Mistake::DescendingSum ? nr - 1 - k : k;
        if (mistake == Mistake::SkippedZeroTerm && ir == ir_flux) continue;
        sum = sum + weights[ir] * psi.block(ir + nr * iz)[block];
      }
      derivative[block + std::size_t(psi.block_size()) * iz] = sum;
    }
  return derivative;
}
#endif

}  // namespace

// The GPU backend's flux surface equals the CPU backend's on the same wave packet, value by value,
// on every plane below Z_asymptotic; the surface spans several chunks of the read-back buffer.
// Bitwise with -fmad=false. The wave holds psi = inf at r_c in one column, so the zero term
// D(r_c) psi(r_c) makes that d psi / dr NaN (0 * inf), as in the Fortran. The same comparison finds
// the porting mistakes of the sum (strict builds): ir descending changes most derivatives, and
// leaving out the zero term loses that one NaN.
TEST(CudaBackend, FluxSurfaceMatchesTheCpu) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s(flux_case, flux_keys());
  const int nr = s.config.input.nr;
  const int n_planes = s.grids.nz_potential;
  const auto cpu = s.cpu();
  const auto gpu = gpu_backend(s);
  const int n_blocks = cpu->wave_packet_on_host().scattering.block_size();
  ASSERT_LT(cuda::flux_surface_chunk_planes(std::size_t(n_blocks) * nr, n_blocks), n_planes);  // several chunks
#ifndef __FAST_MATH__  // -ffinite-math-only: no defined arithmetic with inf
  const int ir_flux = s.grids.ir_flux_analysis.value();
  const int inf_plane = n_planes - 1;  // in the last chunk
  const int inf_block = n_blocks - 1;  // the last thread of a partly used thread block
#endif
  auto load = [&](WavePacket& psi) {
    fill_pseudo_random(psi.scattering.values, 1);
    fill_pseudo_random(psi.specular.values, 2);
#ifndef __FAST_MATH__
    psi.scattering.block(ir_flux + nr * inf_plane)[inf_block] = cplx(std::numeric_limits<double>::infinity(), 0.0);
#endif
  };
  cpu->read_wave_packet(load);
  gpu->read_wave_packet(load);

  FluxSurfaceWave on_cpu = nan_surface(n_blocks, n_planes);
  FluxSurfaceWave on_gpu = nan_surface(n_blocks, n_planes);
  std::as_const(*cpu).flux_surface(on_cpu);
  std::as_const(*gpu).flux_surface(on_gpu);
  EXPECT_EQ(value_differences(on_cpu.psi, on_gpu.psi), 0);
#ifdef __FAST_MATH__
  // FMA in the kernels and in the host code: each side's sum moves by <= nr eps sum |D psi| per
  // component, so CPU and GPU differ by <= 2 nr eps sum |D psi| per component. The factor 4 instead
  // of 2 covers the complex modulus (up to sqrt 2 times a component) and the rounded products in
  // sum |D psi|. No NaN may appear (is_nan reads the bits, which fast-math cannot fold away).
  EXPECT_EQ(nan_count(on_gpu.derivative), 0);
  const ScatteringWave& psi = cpu->wave_packet_on_host().scattering;
  const std::span<const double> weights = s.tables.flux_derivative;
  double worst = 0.0;
  for (int iz = 0; iz < n_planes; ++iz)
    for (int block = 0; block < n_blocks; ++block) {
      double magnitude = 0.0;
      for (int ir = 0; ir < nr; ++ir) magnitude += std::abs(weights[ir] * psi.block(ir + nr * iz)[block]);
      const std::size_t k = on_cpu.index(block, iz);
      worst = nan_aware_max(worst, std::abs(on_cpu.derivative[k] - on_gpu.derivative[k]) /
                                       (4.0 * nr * std::numeric_limits<double>::epsilon() * magnitude));
    }
  EXPECT_LE(worst, 1.0);
#else
  EXPECT_EQ(value_differences(on_cpu.derivative, on_gpu.derivative), 0);
  EXPECT_EQ(nan_count(on_gpu.derivative), 1);  // the real part at (inf_block, inf_plane)
  EXPECT_TRUE(is_nan(on_gpu.derivative[on_gpu.index(inf_block, inf_plane)].real()));

  const ScatteringWave& psi = cpu->wave_packet_on_host().scattering;
  const std::vector<cplx> descending =
      mistaken_derivative(psi, ir_flux, s.tables.flux_derivative, n_planes, Mistake::DescendingSum);
  const std::vector<cplx> no_zero_term =
      mistaken_derivative(psi, ir_flux, s.tables.flux_derivative, n_planes, Mistake::SkippedZeroTerm);
  const int n_values = 2 * n_blocks * n_planes;
  EXPECT_GT(value_differences(descending, on_gpu.derivative), n_values / 2) << "of " << n_values;
  EXPECT_EQ(value_differences(no_zero_term, on_gpu.derivative), 1);
#endif
}

// The chunk loop: the size of the read-back buffer sets how many planes one chunk holds, never the
// values. One chunk, three chunks (the last one shorter) and one plane per chunk give the same bits,
// every element written (in every build: each element is computed on its own). Shapes that do not
// fit, and a buffer too small for one plane, are refused.
TEST(CudaBackend, FluxSurfaceDoesNotDependOnTheChunkSize) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup s(flux_case, flux_keys());
  const InputParameters& p = s.config.input;
  const int ir_flux = s.grids.ir_flux_analysis.value();
  const int n_planes = s.grids.nz_potential;
  const cuda::DeviceLayout layout{std::int64_t(p.nr) * p.nz, s.basis.size(), std::int64_t(p.nx) * p.ny};
  const int n_blocks = static_cast<int>(layout.n_rot * layout.n_xy);
  cuda::DeviceWavePacket psi(layout, p.nr, static_cast<int>(s.grids.z_specular.size()));
  std::vector<cplx> values(layout.size());
  fill_pseudo_random(values, 5);
  psi.scattering.upload(cuda::as_device_complex(std::span<const cplx>(values)));
  const cuda::DeviceBuffer<double> weights(std::span<const double>(s.tables.flux_derivative));

  auto read = [&](int chunk_planes) {
    cuda::DeviceBuffer<cuda::Complex> readback(2 * std::size_t(n_blocks) * chunk_planes);
    FluxSurfaceWave out = nan_surface(n_blocks, n_planes);
    cuda::read_flux_surface(psi, ir_flux, weights, readback, out);
    return out;
  };
  const int third = n_planes / 3 + 1;
  ASSERT_EQ((n_planes + third - 1) / third, 3);
  ASSERT_NE(n_planes % third, 0);  // the last of the three chunks is shorter
  const FluxSurfaceWave one_chunk = read(n_planes);
  EXPECT_EQ(nan_count(one_chunk.psi) + nan_count(one_chunk.derivative), 0);
  for (const int chunk_planes : {third, 1}) {
    SCOPED_TRACE(::testing::Message() << chunk_planes << " planes per chunk");
    const FluxSurfaceWave chunked = read(chunk_planes);
    EXPECT_EQ(value_differences(one_chunk.psi, chunked.psi), 0);
    EXPECT_EQ(value_differences(one_chunk.derivative, chunked.derivative), 0);
  }

  cuda::DeviceBuffer<cuda::Complex> readback(2 * std::size_t(n_blocks));  // one plane
  FluxSurfaceWave fits(n_blocks, n_planes);
  EXPECT_NO_THROW(cuda::read_flux_surface(psi, ir_flux, weights, readback, fits));
  cuda::DeviceBuffer<cuda::Complex> too_small(2 * std::size_t(n_blocks) - 1);
  EXPECT_THROW(cuda::read_flux_surface(psi, ir_flux, weights, too_small, fits), std::invalid_argument);
  EXPECT_THROW(cuda::read_flux_surface(psi, p.nr, weights, readback, fits), std::invalid_argument);
  EXPECT_THROW(cuda::read_flux_surface(psi, -1, weights, readback, fits), std::invalid_argument);
  const cuda::DeviceBuffer<double> fewer_weights(std::span<const double>(s.tables.flux_derivative).first(p.nr - 1));
  EXPECT_THROW(cuda::read_flux_surface(psi, ir_flux, fewer_weights, readback, fits), std::invalid_argument);
  FluxSurfaceWave other_blocks(n_blocks + 1, n_planes);
  EXPECT_THROW(cuda::read_flux_surface(psi, ir_flux, weights, readback, other_blocks), std::invalid_argument);
  FluxSurfaceWave too_many_planes(n_blocks, p.nz + 1);
  EXPECT_THROW(cuda::read_flux_surface(psi, ir_flux, weights, readback, too_many_planes), std::invalid_argument);
}

// The seam's contract, as on the CPU (test_flux_surface.cpp): a run without flux analysis, and a
// slot of another shape than nz_potential planes of n_rot nx ny blocks, are refused.
TEST(CudaBackend, FluxSurfaceRefusesLikeTheCpu) {
  SPODVR_REQUIRE_GPU();
  const BackendSetup without(flux_case);
  FluxSurfaceWave slot(without.basis.size() * without.config.input.nx * without.config.input.ny,
                       without.grids.nz_potential);
  EXPECT_THROW(gpu_backend(without)->flux_surface(slot), std::logic_error);

  const BackendSetup s(flux_case, flux_keys());
  const auto gpu = gpu_backend(s);
  const int blocks = s.basis.size() * s.config.input.nx * s.config.input.ny;
  const int planes = s.grids.nz_potential;
  FluxSurfaceWave fits(blocks, planes);
  EXPECT_NO_THROW(gpu->flux_surface(fits));
  FluxSurfaceWave fewer_planes(blocks, planes - 1);
  EXPECT_THROW(gpu->flux_surface(fewer_planes), std::logic_error);
  FluxSurfaceWave other_blocks(blocks + 1, planes);
  EXPECT_THROW(gpu->flux_surface(other_blocks), std::logic_error);
}
