#pragma once

#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Diagnostics of the initial wave packet, as printed by the Fortran code.
// Indices are 0-based. Fortran: WaveFunction.f90 CoherentStateFunc
struct InitialWavePacketInfo {
  double norm = 0.0;               // norm of the Gaussian on the whole specular Z grid
  double k_space_at_zero = 0.0;    // |psi(k)| at k = 0. Fortran: wfcbeg
  double k_space_at_edge = 0.0;    // |psi(k)| at the largest |k|. Fortran: wfcend
  double k_space_maximum = 0.0;    // max |psi(k)|. Fortran: wfkmax
  int i_low = 0;                   // first point with Z >= ZAnalysis. Fortran: iLow
  int i_high = 0;                  // last point before the specular CAP (or ZAnaSp). Fortran: iHigh
  double norm_inside = 0.0;        // norm between i_low and i_high (must be > 1 - tolsp)
};

struct InitialWavePacket {
  SpecularWave specular;  // the scattering grid part starts at zero
  InitialWavePacketInfo info;
};

// Gaussian wave packet in Z, centred at Z0 and moving towards the surface with momentum kz,
// times the initial rovibrational state in r, on the specular grid:
//   psi(Z, r) = (2 alpha / pi)^(1/4) exp(-alpha (Z - Z0)^2 - i kz Z) phi_{v j}(r).
// Throws InputError (TolSp) if more than tolsp of the norm lies outside [i_low, i_high].
// Fortran: WaveFunction.f90 CoherentStateFunc
[[nodiscard]] InitialWavePacket make_initial_wave_packet(const Config& config, const SpatialGrids& grids,
                                                         const GasPhase& gas, const Fft& fft);

}  // namespace spodvr
