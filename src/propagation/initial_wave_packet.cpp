// Initial wave packet (coherent state in Z times the rovibrational state) on the specular grid.
// Fortran: WaveFunction.f90 CoherentStateFunc
#include "propagation/initial_wave_packet.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <sstream>
#include <vector>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

// Norm sum of |psi|^2 over [first, last], evaluated like the Fortran conjg(psi)*psi sum.
double norm_sum(const std::vector<cplx>& psi, int first, int last) {
  double sum = 0.0;
  for (int i = first; i <= last; ++i) sum = sum + (psi[i].real() * psi[i].real() + psi[i].imag() * psi[i].imag());
  return sum;
}

// Last point allowed to carry the initial wave packet: before the specular absorber, or before
// ZAnaSp if that is fixed in the input. Fortran: CoherentStateFunc (iHigh, nproj = 2)
int last_point_inside(const Config& config, const SpatialGrids& grids) {
  const auto& z = grids.z_specular;
  if (config.derived.fixed_z_analysis_specular && config.input.z_analysis_specular &&
      *config.input.z_analysis_specular > 0.0) {
    for (int iz = 0; iz < z.size(); ++iz)
      if (fortran_compare::rge(z[iz], *config.input.z_analysis_specular)) return iz - 1;
  }
  return grids.z_specular_optical->first - 1;  // Fortran: nschmi
}

}  // namespace

InitialWavePacket make_initial_wave_packet(const Config& config, const SpatialGrids& grids, const GasPhase& gas,
                                           const Fft& fft) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const UniformGrid& z = grids.z_specular;
  const int nzsp = z.size();
  constexpr double pi = Constants::pi;

  // Gaussian in Z. Fortran: SpecWF
  const double alpha = d.wave_packet_alpha;
  const double amplitude = std::pow(2.0 * alpha / pi, 0.25);
  std::vector<cplx> gaussian(nzsp);
  for (int iz = 0; iz < nzsp; ++iz) {
    const double dz = z[iz] - p.z0;
    const cplx exponent(-(alpha * (dz * dz)), -(d.kz * z[iz]));
    gaussian[iz] = amplitude * std::exp(exponent);
  }

  InitialWavePacketInfo info;
  info.norm = norm_sum(gaussian, 0, nzsp - 1) * p.dz;

  // Momentum distribution. Fortran: k_spectrum (sign -1, scale dz / sqrt(2 pi))
  std::vector<cplx> spectrum = gaussian;
  fft.plan_complex({nzsp}, FftSign::Minus, p.dz / std::sqrt(2.0 * pi))->execute(spectrum);
  for (const cplx& c : spectrum) info.k_space_maximum = std::max(info.k_space_maximum, std::abs(c));
  info.k_space_at_zero = std::abs(spectrum[0]);
  info.k_space_at_edge = std::abs(spectrum[nzsp / 2]);

  // The packet must lie between the analysis surface and the absorber.
  for (int iz = 0; iz < nzsp; ++iz)
    if (fortran_compare::rge(z[iz], grids.z_analysis)) {
      info.i_low = iz;
      break;
    }
  info.i_high = last_point_inside(config, grids);
  info.norm_inside = norm_sum(gaussian, info.i_low, info.i_high) * p.dz;
  if ((1.0 - info.norm_inside) > p.tol_specular) {
    std::ostringstream message;
    message << "too much norm of the initial wave packet outside the allowed region (Z = " << z[info.i_low] << " .. "
            << z[info.i_high] << "): 1 - norm = " << (1.0 - info.norm_inside)
            << " > tolsp; move Z0, widen the energy range or enlarge the specular grid";
    throw InputError("TolSp", message.str());
  }

  // psi(Z, r) = gaussian(Z) * phi_{v j}(r).
  const auto phi = gas.states.wave_function(p.v_in, p.j_in);
  SpecularWave specular(p.nr, nzsp);
  for (int iz = 0; iz < nzsp; ++iz)
    for (int ir = 0; ir < p.nr; ++ir) specular(ir, iz) = phi[ir] * gaussian[iz];
  return {std::move(specular), info};
}

}  // namespace spodvr
