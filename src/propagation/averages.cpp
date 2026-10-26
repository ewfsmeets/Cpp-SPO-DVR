// Norm and averages of the wave packet. Fortran: Averages.f90 ComputeAverages, ComplexInnerProduct
#include "propagation/averages.hpp"

#include <cmath>
#include <vector>

namespace spodvr {

namespace {

// |psi|^2 as ABS(psi)**2 (hypot, then squared), like the Fortran code.
double abs_squared(const std::complex<double>& c) {
  const double a = std::abs(c);
  return a * a;
}

}  // namespace

Averages compute_averages(double time, const WavePacket& psi, const Config& config, const SpatialGrids& grids,
                          int analysis_surface_index) {
  const InputParameters& p = config.input;
  const int nr = p.nr;
  const ScatteringWave& scat = psi.scattering;
  const SpecularWave& spec = psi.specular;

  // Density summed over (irot, x, y) at every (Z, r) point. Fortran: complexInnerProd(1, izr)
  std::vector<double> density(scat.n_blocks(), 0.0);
#pragma omp parallel for schedule(static)
  for (int izr = 0; izr < scat.n_blocks(); ++izr) {
    double sum = 0.0;
    for (const auto& c : scat.block(izr)) sum = sum + abs_squared(c);
    density[izr] = sum;
  }

  double norm_scat = 0.0;
  double z_scat = 0.0;
  double r_scat = 0.0;
  for (int izr = 0; izr < scat.n_blocks(); ++izr) {
    const int ir = izr % nr;
    const int iz = izr / nr;
    norm_scat = norm_scat + density[izr];
    z_scat = z_scat + grids.z_scattering[iz] * density[izr];
    r_scat = r_scat + grids.r[ir] * density[izr];
  }
  const double volume_scat = p.dz * p.dr * config.derived.xy_volume_element;  // Fortran: intVolScat
  norm_scat = norm_scat * volume_scat;
  z_scat = z_scat * volume_scat;
  r_scat = r_scat * volume_scat;

  double norm_spec = 0.0;
  double z_spec = 0.0;
  double r_spec = 0.0;
  for (int iz = 0; iz < spec.n_z(); ++iz)
    for (int ir = 0; ir < nr; ++ir) {
      const double d = abs_squared(spec(ir, iz));
      norm_spec = norm_spec + d;
      z_spec = z_spec + grids.z_specular[iz] * d;
      r_spec = r_spec + grids.r[ir] * d;
    }
  const double volume_spec = p.dz * p.dr;  // Fortran: intVolSpec
  norm_spec = norm_spec * volume_spec;
  z_spec = z_spec * volume_spec;
  r_spec = r_spec * volume_spec;

  // Norm on the specular grid between ZAnalysis and the analysis surface.
  double rest = 0.0;
  for (int iz = grids.iz_analysis + 1; iz <= analysis_surface_index; ++iz)
    for (int ir = 0; ir < nr; ++ir) rest = rest + abs_squared(spec(ir, iz));

  Averages a;
  a.time = time;
  a.norm = norm_scat + norm_spec;
  a.z_mean = (z_scat + z_spec) / a.norm;
  a.r_mean = (r_scat + r_spec) / a.norm;
  a.rest_norm = rest * volume_spec;
  return a;
}

}  // namespace spodvr
