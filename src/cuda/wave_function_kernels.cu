// Kernels on the whole scattering wave function: densities for the averages, the BK plane, and
// the conversion between the host layout (irot, ixy, izr) and the device layout (izr, irot, ixy).
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// One thread per (Z, r) point. The sum runs over (irot, ixy) in the host memory order (irot
// fastest), like the CPU block_densities; neighbouring threads read neighbouring izr (coalesced).
// |c|^2 = hypot(re, im)^2 as the CPU's std::abs(c)^2 (CUDA's hypot may differ by an ulp).
__global__ void block_densities_kernel(const Complex* psi, DeviceLayout layout, double* densities) {
  const std::int64_t izr = thread_index();
  if (izr >= layout.n_zr) return;
  double sum = 0.0;
  for (std::int64_t ixy = 0; ixy < layout.n_xy; ++ixy)
    for (std::int64_t irot = 0; irot < layout.n_rot; ++irot) {
      const Complex c = psi[layout(izr, irot, ixy)];
      const double a = hypot(c.re, c.im);
      sum = sum + a * a;
    }
  densities[izr] = sum;
}

// One thread per output element (irot, ixy, ir) of the plane.
__global__ void pack_z_plane_kernel(const Complex* psi, DeviceLayout layout, int nr, int iz, Complex* plane) {
  const std::int64_t i = thread_index();
  const std::int64_t k_count = layout.n_rot * layout.n_xy;  // values per (Z, r) point
  if (i >= k_count * nr) return;
  const std::int64_t k = i % k_count;  // irot + n_rot * ixy
  const std::int64_t ir = i / k_count;
  plane[i] = psi[(ir + std::int64_t(nr) * iz) + layout.n_zr * k];
}

// One thread per element; i runs izr fastest, so the device-layout side is coalesced (the host
// chunk side is strided; these copies run once per run).
__global__ void host_chunk_to_device_kernel(const Complex* chunk, DeviceLayout layout, std::int64_t izr_first,
                                            std::int64_t n_points, Complex* psi) {
  const std::int64_t i = thread_index();
  const std::int64_t k_count = layout.n_rot * layout.n_xy;
  if (i >= k_count * n_points) return;
  const std::int64_t point = i % n_points;
  const std::int64_t k = i / n_points;
  psi[(izr_first + point) + layout.n_zr * k] = chunk[k + k_count * point];
}

__global__ void device_to_host_chunk_kernel(const Complex* psi, DeviceLayout layout, std::int64_t izr_first,
                                            std::int64_t n_points, Complex* chunk) {
  const std::int64_t i = thread_index();
  const std::int64_t k_count = layout.n_rot * layout.n_xy;
  if (i >= k_count * n_points) return;
  const std::int64_t point = i % n_points;
  const std::int64_t k = i / n_points;
  chunk[k + k_count * point] = psi[(izr_first + point) + layout.n_zr * k];
}

}  // namespace

// Fortran: Averages.f90 complexInnerProd(1, izr)
void block_densities(const Complex* psi, const DeviceLayout& layout, double* densities) {
  block_densities_kernel<<<blocks_for(layout.n_zr), threads_per_block>>>(psi, layout, densities);
  check_kernel_launch("block_densities");
}

void pack_z_plane(const Complex* psi, const DeviceLayout& layout, int nr, int iz, Complex* plane) {
  const std::int64_t n = layout.n_rot * layout.n_xy * nr;
  pack_z_plane_kernel<<<blocks_for(n), threads_per_block>>>(psi, layout, nr, iz, plane);
  check_kernel_launch("pack_z_plane");
}

void host_chunk_to_device(const Complex* chunk, const DeviceLayout& layout, std::int64_t izr_first,
                          std::int64_t n_points, Complex* psi) {
  const std::int64_t n = layout.n_rot * layout.n_xy * n_points;
  host_chunk_to_device_kernel<<<blocks_for(n), threads_per_block>>>(chunk, layout, izr_first, n_points, psi);
  check_kernel_launch("host_chunk_to_device");
}

void device_to_host_chunk(const Complex* psi, const DeviceLayout& layout, std::int64_t izr_first,
                          std::int64_t n_points, Complex* chunk) {
  const std::int64_t n = layout.n_rot * layout.n_xy * n_points;
  device_to_host_chunk_kernel<<<blocks_for(n), threads_per_block>>>(psi, layout, izr_first, n_points, chunk);
  check_kernel_launch("device_to_host_chunk");
}

}  // namespace spodvr::cuda
