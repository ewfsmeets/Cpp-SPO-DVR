// The GPU propagation backend: the wave packet stays in GPU memory; per analysis step only the
// specular wave, the BK plane and (for the averages) the per-point densities are copied back.
// CPU twin: propagation/cpu_backend.cpp. Fortran: SplitOperator.f90 Propagate
#include "cuda/cuda_backend.hpp"

#include <cuda_runtime_api.h>

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <optional>
#include <utility>

#include "cuda/cuda_error.hpp"
#include "cuda/cuda_fft.hpp"
#include "cuda/device_kinetic_operator.hpp"
#include "cuda/device_potential_operator.hpp"
#include "cuda/device_wave_packet.hpp"
#include "cuda/kernels.hpp"
#include "propagation/angular_transform.hpp"

namespace spodvr::cuda {

namespace {

using cplx = std::complex<double>;

// Limit of the cuFFT work area (one area shared by all plans) and of the staging buffer used to
// copy the wave function between the host and the GPU layout.
constexpr std::size_t fft_work_limit = std::size_t(512) << 20;  // 512 MiB
constexpr std::size_t staging_bytes = std::size_t(64) << 20;    // 64 MiB

ScatteringShape scattering_shape(const Config& config, int n_rot) {
  const InputParameters& p = config.input;
  return {n_rot, p.nx, p.ny, p.nr, p.nz};
}

DeviceLayout device_layout(const ScatteringShape& s) {
  return {std::int64_t(s.nr) * s.nz, s.n_rot, std::int64_t(s.nx) * s.ny};
}

// (Z, r) points per staging chunk (at least one).
std::int64_t staging_points(const DeviceLayout& layout) {
  const std::size_t point_bytes = sizeof(Complex) * layout.n_rot * layout.n_xy;
  return std::max<std::int64_t>(1, staging_bytes / point_bytes);
}

class CudaBackend final : public PropagationBackend {
 public:
  CudaBackend(const BackendInput& in, PotentialGrid potential, SpecularWave initial);

  void read_wave_packet(const std::function<void(WavePacket&)>& read) override;
  [[nodiscard]] const WavePacket& wave_packet_on_host() override;
  void propagate_analysis_step() override;
  [[nodiscard]] const SpecularWave& specular() const override { return host_specular_; }
  [[nodiscard]] std::span<const std::complex<double>> z_plane(int iz) const override;
  [[nodiscard]] std::vector<double> block_densities() const override;

 private:
  void upload(const WavePacket& host);
  void download(WavePacket& host) const;
  // One T(dt') block (specular grid, T_Z + T_r, T_xy), followed by V unless kinetic_only.
  void kinetic_and_potential_step(TimeStep step, bool kinetic_only);
  void share_fft_work_area();

  ScatteringShape shape_;
  int steps_per_analysis_;  // Fortran: nanspo
  DeviceWavePacket psi_;
  DeviceKineticOperator kinetic_;
  DeviceSpecularOperator specular_;
  DevicePotentialOperator potential_;
  DeviceBuffer<std::byte> fft_work_;
  SpecularWave host_specular_;  // copy of the specular wave after the last analysis step
  // Buffers of the per-step copies (filled on request, hence mutable).
  mutable DeviceBuffer<Complex> plane_;
  mutable std::vector<cplx> host_plane_;
  mutable DeviceBuffer<double> densities_;
  std::optional<WavePacket> host_psi_;  // only for wave_packet_on_host (checkpoint, dump)
};

DevicePotentialOperator make_potential_operator(const BackendInput& in, PotentialGrid& potential) {
  premultiply_with_time_step(potential.values, in.config.input.time_step);
  const AngularTransform transform(in.angular, in.basis, in.fft);  // its tables only
  return DevicePotentialOperator(scattering_shape(in.config, in.basis.size()), in.grids.nz_potential,
                                 in.basis.initial_index(), in.tables, transform, potential.values, potential.symmetry,
                                 fft_work_limit);
}

CudaBackend::CudaBackend(const BackendInput& in, PotentialGrid potential, SpecularWave initial)
    : shape_(scattering_shape(in.config, in.basis.size())),
      steps_per_analysis_(in.config.derived.steps_per_analysis),
      psi_(device_layout(shape_), in.config.input.nr, in.grids.z_specular.size()),
      kinetic_(shape_, in.tables, fft_work_limit),
      specular_(in.config.input.nr, in.grids.z_specular.size(), in.grids.nz_potential, in.basis.initial_index(),
                *in.grids.z_specular_optical, in.grids.r_optical, in.tables, fft_work_limit),
      potential_(make_potential_operator(in, potential)),
      host_specular_(std::move(initial)),
      plane_(std::size_t(shape_.n_rot) * shape_.nx * shape_.ny * shape_.nr),
      host_plane_(plane_.size()),
      densities_(std::size_t(psi_.layout.n_zr)) {
  share_fft_work_area();
  // The initial wave packet lives on the specular grid; the scattering grid starts empty.
  psi_.scattering.set_zero();
  psi_.specular.upload(as_device_complex(std::span<const cplx>(host_specular_.values)));
}

void CudaBackend::share_fft_work_area() {
  std::vector<DeviceFftPlan*> plans = kinetic_.fft_plans();
  for (DeviceFftPlan* p : specular_.fft_plans()) plans.push_back(p);
  for (DeviceFftPlan* p : potential_.fft_plans()) plans.push_back(p);
  std::size_t bytes = 1;
  for (const DeviceFftPlan* p : plans) bytes = std::max(bytes, p->work_size());
  fft_work_ = DeviceBuffer<std::byte>(bytes);
  for (DeviceFftPlan* p : plans) p->set_work_area(fft_work_.data());
}

// Host layout (irot, ixy, izr) -> device layout, chunk by chunk of (Z, r) points through a
// staging buffer on the GPU (the reordering is done there).
void CudaBackend::upload(const WavePacket& host) {
  const DeviceLayout& layout = psi_.layout;
  const std::int64_t chunk = staging_points(layout);
  const std::int64_t point_values = layout.n_rot * layout.n_xy;
  DeviceBuffer<Complex> staging(std::size_t(chunk * point_values));
  const auto values = as_device_complex(std::span<const cplx>(host.scattering.values));
  for (std::int64_t first = 0; first < layout.n_zr; first += chunk) {
    const std::int64_t n = std::min(chunk, layout.n_zr - first);
    staging.upload(values.subspan(first * point_values, n * point_values));
    host_chunk_to_device(staging.data(), layout, first, n, psi_.scattering.data());
  }
  psi_.specular.upload(as_device_complex(std::span<const cplx>(host.specular.values)));
  host_specular_ = host.specular;
}

void CudaBackend::download(WavePacket& host) const {
  const DeviceLayout& layout = psi_.layout;
  const std::int64_t chunk = staging_points(layout);
  const std::int64_t point_values = layout.n_rot * layout.n_xy;
  DeviceBuffer<Complex> staging(std::size_t(chunk * point_values));
  const auto values = as_device_complex(std::span<cplx>(host.scattering.values));
  for (std::int64_t first = 0; first < layout.n_zr; first += chunk) {
    const std::int64_t n = std::min(chunk, layout.n_zr - first);
    device_to_host_chunk(psi_.scattering.data(), layout, first, n, staging.data());
    staging.download(values.subspan(first * point_values, n * point_values));
  }
  psi_.specular.download(as_device_complex(std::span<cplx>(host.specular.values)));
}

void CudaBackend::read_wave_packet(const std::function<void(WavePacket&)>& read) {
  WavePacket host{ScatteringWave(shape_.n_rot, shape_.nx * shape_.ny, shape_.nr * shape_.nz),
                  SpecularWave(psi_.nr, psi_.nz_specular)};
  read(host);
  upload(host);
}

const WavePacket& CudaBackend::wave_packet_on_host() {
  if (!host_psi_)
    host_psi_.emplace(WavePacket{ScatteringWave(shape_.n_rot, shape_.nx * shape_.ny, shape_.nr * shape_.nz),
                                 SpecularWave(psi_.nr, psi_.nz_specular)});
  download(*host_psi_);
  return *host_psi_;
}

// Fortran: SplitOperator.f90 KineticAndPotentialOp (lastStep: kinetic only)
void CudaBackend::kinetic_and_potential_step(TimeStep step, bool kinetic_only) {
  specular_.apply(psi_, step, kinetic_only);
  kinetic_.apply_zr(psi_, step);
  kinetic_.apply_xy(psi_, step);
  if (!kinetic_only) potential_.apply(psi_);
}

// T(dt/2) V (T V)^(n-1) T(dt/2), as SplitOperatorPropagator::propagate_analysis_step.
void CudaBackend::propagate_analysis_step() {
  kinetic_and_potential_step(TimeStep::Half, false);
  for (int step = 1; step < steps_per_analysis_; ++step) kinetic_and_potential_step(TimeStep::Full, false);
  kinetic_and_potential_step(TimeStep::Half, true);
  // The analysis surface and the BK coefficients need the specular wave on the host.
  psi_.specular.download(as_device_complex(std::span<cplx>(host_specular_.values)));
}

std::span<const std::complex<double>> CudaBackend::z_plane(int iz) const {
  pack_z_plane(psi_.scattering.data(), psi_.layout, shape_.nr, iz, plane_.data());
  plane_.download(as_device_complex(std::span<cplx>(host_plane_)));
  return host_plane_;
}

std::vector<double> CudaBackend::block_densities() const {
  cuda::block_densities(psi_.scattering.data(), psi_.layout, densities_.data());
  std::vector<double> out(densities_.size());
  densities_.download(out);
  return out;
}

}  // namespace

void initialise_cuda() {
  setenv("CUDA_MODULE_LOADING", "LAZY", 0);
  int count = 0;
  if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) throw CudaError("no CUDA-capable GPU found");
  check_cuda_status(cudaSetDevice(0), "cudaSetDevice");
}

std::string cuda_device_summary() {
  cudaDeviceProp p{};
  check_cuda_status(cudaGetDeviceProperties(&p, 0), "cudaGetDeviceProperties");
  int runtime = 0;
  int driver = 0;
  cudaRuntimeGetVersion(&runtime);
  cudaDriverGetVersion(&driver);
  auto version = [](int v) { return std::to_string(v / 1000) + "." + std::to_string(v % 1000 / 10); };
  char memory[32];
  std::snprintf(memory, sizeof(memory), "%.1f GiB", double(p.totalGlobalMem) / double(1 << 30));
  return std::string("GPU ") + p.name + " (" + memory + ", compute capability " + std::to_string(p.major) + "." +
         std::to_string(p.minor) + "), CUDA runtime " + version(runtime) + ", driver " + version(driver) + ", " +
         cufft_version();
}

std::size_t free_gpu_memory() {
  std::size_t free = 0;
  std::size_t total = 0;
  check_cuda_status(cudaMemGetInfo(&free, &total), "cudaMemGetInfo");
  return free;
}

std::size_t estimate_gpu_memory(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                                int n_representatives) {
  const InputParameters& p = config.input;
  const RotationalBasis basis(angular, p.j_in, p.mj_in);
  const ScatteringShape shape = scattering_shape(config, basis.size());
  const DeviceLayout layout = device_layout(shape);
  const std::size_t psi = sizeof(Complex) * layout.size();
  const std::size_t specular = sizeof(Complex) * std::size_t(p.nr) * grids.z_specular.size();
  const std::size_t plane = sizeof(Complex) * std::size_t(layout.n_rot) * layout.n_xy * p.nr;
  const std::size_t potential = DevicePotentialOperator::memory_bytes(shape, grids.nz_potential, angular.n_theta(),
                                                                      angular.n_phi(), n_representatives);
  return psi + specular + plane + potential + fft_work_limit + staging_bytes;
}

std::unique_ptr<PropagationBackend> make_cuda_backend(const BackendInput& input, PotentialGrid potential,
                                                      SpecularWave initial) {
  return std::make_unique<CudaBackend>(input, std::move(potential), std::move(initial));
}

}  // namespace spodvr::cuda
