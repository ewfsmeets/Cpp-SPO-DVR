// psi and its r-derivative at the flux surface r = r_c on the GPU, copied to the host chunk by chunk
// of planes. CPU twin: propagation/flux_surface.cpp. Fortran: WaveFunction.f90 getWaveFunctionAtZr,
// getWaveFuncDerivAtZr
#include "cuda/device_flux_surface.hpp"

#include <algorithm>
#include <complex>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>

#include "cuda/kernels.hpp"

namespace spodvr::cuda {

namespace {

using cplx = std::complex<double>;

// Throws std::invalid_argument unless out (n_blocks per plane x n_planes) lies on psi's grid with nr
// weights, ir_flux is a point of it, and readback holds at least one plane.
void require_fitting_shapes(const DeviceWavePacket& psi, int ir_flux, std::size_t n_weights, std::size_t readback_size,
                            const FluxSurfaceWave& out) {
  const DeviceLayout& layout = psi.layout;
  const int nr = psi.nr;
  const bool fits = nr > 0 && layout.n_zr % nr == 0 && ir_flux >= 0 && ir_flux < nr && n_weights == std::size_t(nr) &&
                    out.n_blocks() == layout.n_rot * layout.n_xy && out.n_planes() <= layout.n_zr / nr &&
                    flux_surface_chunk_planes(readback_size, out.n_blocks()) > 0;
  if (!fits)
    throw std::invalid_argument(std::format(
        "read_flux_surface: a surface of {} blocks x {} planes at ir = {} with {} weights and a read-back buffer of "
        "{} values does not fit a wave function of {} x {} blocks x {} (Z, r) points with nr = {}",
        out.n_blocks(), out.n_planes(), ir_flux, n_weights, readback_size, layout.n_rot, layout.n_xy, layout.n_zr, nr));
}

}  // namespace

int flux_surface_chunk_planes(std::size_t readback_size, int n_blocks) {
  if (n_blocks < 1) return 0;
  const std::size_t planes = readback_size / (2 * std::size_t(n_blocks));
  return static_cast<int>(std::min<std::size_t>(planes, std::numeric_limits<int>::max()));
}

void read_flux_surface(const DeviceWavePacket& psi, int ir_flux, const DeviceBuffer<double>& weights,
                       DeviceBuffer<Complex>& readback, FluxSurfaceWave& out) {
  require_fitting_shapes(psi, ir_flux, weights.size(), readback.size(), out);
  const int n_blocks = out.n_blocks();
  const int chunk_planes = flux_surface_chunk_planes(readback.size(), n_blocks);
  for (int first = 0; first < out.n_planes(); first += chunk_planes) {
    const int n_planes = std::min(chunk_planes, out.n_planes() - first);
    flux_surface_chunk(psi.scattering.data(), psi.layout, psi.nr, ir_flux, weights.data(), first, n_planes,
                       readback.data());
    // The chunk is [psi_c | d psi / dr], each in host layout: two contiguous copies into the slot.
    // cudaMemcpy waits for the kernel and returns when the copy is done, before the next chunk.
    const std::size_t values = std::size_t(n_blocks) * n_planes;
    readback.download(as_device_complex(std::span<cplx>(out.psi).subspan(out.index(0, first), values)));
    readback.download(as_device_complex(std::span<cplx>(out.derivative).subspan(out.index(0, first), values)), values);
  }
}

}  // namespace spodvr::cuda
