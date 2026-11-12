#pragma once

// psi and d psi / dr at the flux surface r = r_c on the GPU, read back to the host in chunks of
// planes: the CUDA twin of propagation/flux_surface.hpp. Host code; needs no CUDA headers.

#include <cstddef>

#include "analysis/flux_analysis.hpp"
#include "cuda/device_buffer.hpp"
#include "cuda/device_wave_packet.hpp"

namespace spodvr::cuda {

// Planes of the flux surface that one chunk in a read-back buffer of readback_size values holds:
// each plane takes psi and d psi / dr of its n_blocks blocks. 0 if not even one plane fits.
[[nodiscard]] int flux_surface_chunk_planes(std::size_t readback_size, int n_blocks);

// psi and d psi / dr at the flux surface r = r(ir_flux) (0-based) of the device wave packet psi, on
// the planes iz < out.n_planes(), into every element of out (host layout (block, iz)). weights: the
// nr weights D(ir) (PropagationTables::flux_derivative) on the GPU. Chunk by chunk of
// flux_surface_chunk_planes(readback.size(), n_blocks) planes: flux_surface_chunk writes a chunk into
// readback, then two copies move its psi_c and d psi / dr into out. Every element is computed on its
// own in the CPU twin's order (flux_surface), so the result does not depend on the chunk size and,
// with -fmad=false, equals the CPU's bitwise. Throws std::invalid_argument if the shapes do not fit
// or readback cannot hold one plane.
// Fortran: WaveFunction.f90 getWaveFunctionAtZr, getWaveFuncDerivAtZr
void read_flux_surface(const DeviceWavePacket& psi, int ir_flux, const DeviceBuffer<double>& weights,
                       DeviceBuffer<Complex>& readback, FluxSurfaceWave& out);

}  // namespace spodvr::cuda
