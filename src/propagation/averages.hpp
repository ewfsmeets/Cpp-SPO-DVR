#pragma once

#include <span>
#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Norm and position averages of the wave packet (scattering plus specular grid).
// Fortran: Averages.f90 ComputeAverages (nproj = 2 output columns)
struct Averages {
  double time = 0.0;
  double norm = 0.0;
  double z_mean = 0.0;     // <Z>
  double r_mean = 0.0;     // <r>
  double rest_norm = 0.0;  // specular-grid norm between ZAnalysis and the analysis surface
};

// sum |psi|^2 over (irot, ixy) at every (Z, r) point izr, in memory order (the points in
// parallel). Fortran: Averages.f90 complexInnerProd(1, izr)
[[nodiscard]] std::vector<double> block_densities(const ScatteringWave& psi);

// densities: block_densities of the scattering wave; analysis_surface_index:
// AnalysisSurface::index(). The sums run in a fixed order, so the result does not depend on the
// number of threads (nor on the backend that computed the densities).
[[nodiscard]] Averages compute_averages(double time, std::span<const double> densities, const SpecularWave& specular,
                                        const Config& config, const SpatialGrids& grids, int analysis_surface_index);

// The same for a wave packet in host memory.
[[nodiscard]] Averages compute_averages(double time, const WavePacket& psi, const Config& config,
                                        const SpatialGrids& grids, int analysis_surface_index);

}  // namespace spodvr
