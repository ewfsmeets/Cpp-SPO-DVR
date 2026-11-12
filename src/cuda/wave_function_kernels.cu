// Kernels on the whole scattering wave function: densities for the averages, the BK plane, the flux
// surface, and the conversion between the host layout (irot, ixy, izr) and the device layout (izr, irot, ixy).
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
      const double abs_c = hypot(c.re, c.im);
      sum = sum + abs_c * abs_c;
    }
  densities[izr] = sum;
}

// One thread per output element (irot, ixy, ir) of the plane.
__global__ void pack_z_plane_kernel(const Complex* psi, DeviceLayout layout, int nr, int iz, Complex* plane) {
  const std::int64_t i = thread_index();
  const std::int64_t n_rot_xy = layout.n_rot * layout.n_xy;  // values per (Z, r) point: (irot, ixy)
  if (i >= n_rot_xy * nr) return;
  const std::int64_t rot_xy = i % n_rot_xy;  // irot + n_rot * ixy
  const std::int64_t ir = i / n_rot_xy;
  plane[i] = psi[(ir + std::int64_t(nr) * iz) + layout.n_zr * rot_xy];
}

// x over the blocks (irot, ixy) of a plane, y over the planes of the chunk (grid-stride). The block
// index b = irot + n_rot ixy is the host order, and layout(izr, irot, ixy) = izr + n_zr b, so each
// thread walks its own contiguous r column. The reads are not coalesced (neighbouring threads read
// n_zr * 16 B apart), the writes are (neighbouring threads, neighbouring blocks). About 110 GB/s on
// the DZ grid (docs/cuda.md, "Flux surface and flush"), probably because of the uncoalesced reads and
// the one serial chain of nr terms per thread (inferred, not profiled).
// The sum from 0, ir ascending, incl. the zero term D(ir_flux) psi(ir_flux): flux_r_derivative's order.
__global__ void flux_surface_kernel(const Complex* psi, DeviceLayout layout, int nr, int ir_flux, const double* weights,
                                    int first_plane, int n_planes, Complex* chunk) {
  const std::int64_t block = thread_index();
  const std::int64_t n_blocks = layout.n_rot * layout.n_xy;
  if (block >= n_blocks) return;
  const std::int64_t derivative_offset = n_blocks * n_planes;  // chunk = [psi_c | d psi / dr]
  for (std::int64_t plane = blockIdx.y; plane < n_planes; plane += gridDim.y) {
    const std::int64_t column = std::int64_t(nr) * (first_plane + plane) + layout.n_zr * block;  // izr of ir = 0
    Complex derivative{0.0, 0.0};
    for (int ir = 0; ir < nr; ++ir) derivative = derivative + weights[ir] * psi[column + ir];
    const std::int64_t out = block + n_blocks * plane;  // (block, iz - first_plane) in host layout
    chunk[out] = psi[column + ir_flux];
    chunk[derivative_offset + out] = derivative;
  }
}

// One thread per element; i runs izr fastest, so the device-layout side is coalesced (the host
// chunk side is strided; these copies run once per run).
__global__ void host_chunk_to_device_kernel(const Complex* chunk, DeviceLayout layout, std::int64_t izr_first,
                                            std::int64_t n_points, Complex* psi) {
  const std::int64_t i = thread_index();
  const std::int64_t n_rot_xy = layout.n_rot * layout.n_xy;
  if (i >= n_rot_xy * n_points) return;
  const std::int64_t point = i % n_points;
  const std::int64_t rot_xy = i / n_points;
  psi[(izr_first + point) + layout.n_zr * rot_xy] = chunk[rot_xy + n_rot_xy * point];
}

__global__ void device_to_host_chunk_kernel(const Complex* psi, DeviceLayout layout, std::int64_t izr_first,
                                            std::int64_t n_points, Complex* chunk) {
  const std::int64_t i = thread_index();
  const std::int64_t n_rot_xy = layout.n_rot * layout.n_xy;
  if (i >= n_rot_xy * n_points) return;
  const std::int64_t point = i % n_points;
  const std::int64_t rot_xy = i / n_points;
  chunk[rot_xy + n_rot_xy * point] = psi[(izr_first + point) + layout.n_zr * rot_xy];
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

// Fortran: WaveFunction.f90 getWaveFunctionAtZr, getWaveFuncDerivAtZr (one iz each)
void flux_surface_chunk(const Complex* psi, const DeviceLayout& layout, int nr, int ir_flux, const double* weights,
                        int first_plane, int n_planes, Complex* chunk) {
  flux_surface_kernel<<<grid_2d(layout.n_rot * layout.n_xy, n_planes), threads_per_block>>>(
      psi, layout, nr, ir_flux, weights, first_plane, n_planes, chunk);
  check_kernel_launch("flux_surface_chunk");
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
