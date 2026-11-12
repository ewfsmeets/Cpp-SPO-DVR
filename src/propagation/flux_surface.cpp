// psi and its r-derivative at the flux surface r = r_c of the scattering wave, on the host (OpenMP).
// CUDA twin: CudaBackend::flux_surface (not ported yet, it refuses). Fortran: WaveFunction.f90
// getWaveFunctionAtZr, getWaveFuncDerivAtZr
#include "propagation/flux_surface.hpp"

#include <algorithm>
#include <complex>
#include <format>
#include <stdexcept>

namespace spodvr {

namespace {

using cplx = std::complex<double>;

// Throws std::invalid_argument unless out (n_blocks per plane x n_planes) lies on psi's grid with
// nr = weights points in r and ir_flux is one of them.
void require_fitting_shapes(const ScatteringWave& psi, int ir_flux, int nr, const FluxSurfaceWave& out) {
  const int n_zr = psi.n_blocks();  // (Z, r) points of psi
  const bool fits = nr > 0 && n_zr % nr == 0 && ir_flux >= 0 && ir_flux < nr && out.n_blocks() == psi.block_size() &&
                    out.n_planes() <= n_zr / nr;
  if (!fits)
    throw std::invalid_argument(std::format(
        "flux_surface: a surface of {} blocks x {} planes at ir = {} does not fit a wave function of {} blocks x "
        "{} (Z, r) points with {} r-derivative weights",
        out.n_blocks(), out.n_planes(), ir_flux, psi.block_size(), n_zr, nr));
}

// The (Z, r) point izr = ir + nr iz of the scattering wave (the index of ScatteringWave::block).
int zr_point(int ir, int iz, int nr) { return ir + nr * iz; }

// d psi / dr at r_c of the blocks first_block .. first_block + derivative.size() - 1 on plane iz,
// with the nr weights D(ir).
void sum_derivative_tile(const ScatteringWave& psi, int iz, int first_block, std::span<const double> weights,
                         std::span<cplx> derivative) {
  const int nr = static_cast<int>(weights.size());
  std::ranges::fill(derivative, cplx(0.0, 0.0));
  // From ir = 0 upward, incl. the zero D(ir_c) term: flux_r_derivative's (and the Fortran's) order
  // per element, so d psi / dr is bitwise equal to it. Must not be reordered or fused into FMA (strict
  // builds; build-fast has a tolerance path in the tests).
  for (int ir = 0; ir < nr; ++ir) {
    const double weight = weights[ir];
    const std::span<const cplx> values = psi.block(zr_point(ir, iz, nr)).subspan(first_block, derivative.size());
    for (std::size_t b = 0; b < derivative.size(); ++b) derivative[b] = derivative[b] + weight * values[b];
  }
}

}  // namespace

// Fortran: WaveFunction.f90 getWaveFunctionAtZr (wfAtZr = ScatWave(:, :, irFluxAnalysis, iz)) and
// getWaveFuncDerivAtZr (wf_deriv), called per iz by FluxAnalyzor.f90 storeFluxCoeffs
void flux_surface(const ScatteringWave& psi, int ir_flux, std::span<const double> weights, FluxSurfaceWave& out) {
  const int nr = static_cast<int>(weights.size());
  require_fitting_shapes(psi, ir_flux, nr, out);
  const int n_blocks = out.n_blocks();
  // Static over iz: the partition that zero-filled the slot and that flushes it (first touch). The
  // reads of psi span all r of a plane, not the propagator's (Z, r) groups: like apply_zr they are
  // not NUMA-local. Every element is computed by one thread in a fixed order: no thread dependence.
#pragma omp parallel for schedule(static)
  for (int iz = 0; iz < out.n_planes(); ++iz) {
    const std::span<cplx> psi_c = std::span(out.psi).subspan(out.index(0, iz), n_blocks);
    std::ranges::copy(psi.block(zr_point(ir_flux, iz, nr)), psi_c.begin());
    for (int first = 0; first < n_blocks; first += flux_surface_tile_blocks) {
      const int tile_size = std::min(flux_surface_tile_blocks, n_blocks - first);
      sum_derivative_tile(psi, iz, first, weights, std::span(out.derivative).subspan(out.index(first, iz), tile_size));
    }
  }
}

}  // namespace spodvr
