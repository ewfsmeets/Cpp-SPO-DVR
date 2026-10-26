#pragma once

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

// analysis_surface_index: AnalysisSurface::index(). The sums run in a fixed order, so the
// result does not depend on the number of threads.
[[nodiscard]] Averages compute_averages(double time, const WavePacket& psi, const Config& config,
                                        const SpatialGrids& grids, int analysis_surface_index);

}  // namespace spodvr
