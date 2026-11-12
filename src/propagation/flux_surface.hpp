#pragma once

#include <span>
#include <string_view>

#include "analysis/flux_analysis.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Blocks per tile of flux_surface's derivative sum: d psi / dr of one tile (16 KiB) stays in the L1
// cache while all nr points in r are added to it. Chosen from cache sizes (a plane of the DZ grid,
// 435 rotational states x 12 x 12 blocks x 16 B = 1.0 MB, docs/cuda.md, exceeds a core's L2), not measured.
inline constexpr int flux_surface_tile_blocks = 1024;

// psi and d psi / dr at the flux surface r = r(ir_flux) (ir_flux 0-based) of the scattering wave
// psi, on the planes iz < out.n_planes() (nz_potential), into every element of out. weights: the
// nr weights D(ir) at r_c (PropagationTables::flux_derivative). Each d psi / dr is summed in
// flux_r_derivative's order, so it equals that function on the column of psi(irot, ixy, :, iz); the
// result does not depend on the thread count. Throws std::invalid_argument if the shapes do not fit.
// The CPU backend's flux surface (CpuBackend::flux_surface; GPU twin cuda::read_flux_surface).
// Fortran: WaveFunction.f90 getWaveFunctionAtZr, getWaveFuncDerivAtZr
void flux_surface(const ScatteringWave& psi, int ir_flux, std::span<const double> weights, FluxSurfaceWave& out);

// The contract of PropagationBackend::flux_surface, checked by both backends (`backend` names the
// caller in the message): throws std::logic_error without flux analysis, or unless out has n_planes
// (nz_potential) planes of n_blocks (n_rot nx ny) blocks.
void require_flux_surface_slot(std::string_view backend, bool has_flux_analysis, const FluxSurfaceWave& out,
                               int n_blocks, int n_planes);

}  // namespace spodvr
