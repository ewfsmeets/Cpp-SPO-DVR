// The building blocks of the flux analysis: its energies with the wave packet's b(k) at them, and
// the r-derivative at the flux surface. Fortran: FluxAnalyzor.f90 setupFluxAnalyzor,
// storeTimeIndependentWF; WaveFunction.f90 getWaveFuncDerivAtZr
#include "analysis/flux_analysis.hpp"

#include <cmath>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/fortran_math.hpp"
#include "core/unit_conversion.hpp"
#include "input/wave_packet_params.hpp"

namespace spodvr {

namespace {

// FluxAnalyzor.f90 converts its energies and the wave packet's with Constants_eV2Ha (flux Q5).
constexpr const UnitConversion& flux_units = units::precise;

// The last energy "reaches" EMaxFluxAna within this fraction of a step: far above the rounding of
// E_min + n dE (~1e-16 relative), far below a step.
constexpr double end_point_tolerance = 1e-6;

// The warning of flux Q3, or an empty string if the last energy reaches EMaxFluxAna.
std::string end_point_warning(const InputParameters& p, const FluxEnergies& e) {
  const double e_max = flux_units.ev_to_hartree(p.e_max_flux_analysis);
  const double step = flux_units.ev_to_hartree(p.delta_e_flux_analysis);
  const double last = e.z_energy.back();
  if (e_max - last <= end_point_tolerance * step) return {};
  return std::format(
      "the flux energies end at {:.6g} eV, below EMaxFluxAna = {:.6g} eV: their number "
      "INT((EMaxFluxAna - EMinFluxAna) / DeltaEFluxAna) + 1 = {} is the Fortran's, which drops the end point "
      "even when DeltaEFluxAna divides the range (by rounding)",
      flux_units.hartree_to_ev(last), p.e_max_flux_analysis, e.size());
}

}  // namespace

// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (numScatEnergies)
int flux_energy_count(const InputParameters& p) {
  if (!p.do_flux_analysis) return 0;
  const double e_min = flux_units.ev_to_hartree(p.e_min_flux_analysis);
  const double e_max = flux_units.ev_to_hartree(p.e_max_flux_analysis);
  const double step = flux_units.ev_to_hartree(p.delta_e_flux_analysis);
  const double quotient = (e_max - e_min) / step;
  if (!(quotient < max_flux_energies))  // also NaN: bound it before the conversion to int
    throw InputError("DeltaEFluxAna",
                     std::format("{} eV gives {:.3g} flux energies, more than {}; choose a larger DeltaEFluxAna",
                                 p.delta_e_flux_analysis, quotient + 1.0, max_flux_energies));
  return static_cast<int>(quotient) + 1;
}

// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (ZScatEnergy, totalScatEnergy) and
// storeTimeIndependentWF (kZ_initial, bkz)
FluxEnergies make_flux_energies(const Config& config, double parallel_energy, double initial_state_energy) {
  const InputParameters& p = config.input;
  const double mass = config.derived.total_mass;  // Fortran: rmass
  const double e_min = flux_units.ev_to_hartree(p.e_min_flux_analysis);
  const double step = flux_units.ev_to_hartree(p.delta_e_flux_analysis);
  const double packet_low = flux_units.ev_to_hartree(p.e_min_wave_packet);
  const double packet_high = flux_units.ev_to_hartree(p.e_max_wave_packet);

  FluxEnergies e;
  const int count = flux_energy_count(p);
  for (int n = 0; n < count; ++n) {
    const double z_energy = e_min + n * step;
    const double momentum = std::sqrt(2.0 * mass * z_energy);
    e.z_energy.push_back(z_energy);
    e.total.push_back((z_energy + parallel_energy) + initial_state_energy);
    e.momentum.push_back(momentum);
    e.amplitude.push_back(wave_packet_momentum_amplitude(momentum, packet_low, packet_high, mass, p.z0));
  }
  if (count > 0) {
    std::string warning = end_point_warning(p, e);
    if (!warning.empty()) e.warnings.push_back(std::move(warning));
  }
  return e;
}

// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr (the ddr table, column irFluxAnalysis)
std::vector<double> flux_derivative_weights(int nr, double dr, int ir_flux) {
  const double c = Constants::pi / nr;  // Fortran: const
  std::vector<double> weights(nr);
  for (int ir = 0; ir < nr; ++ir) {
    const int d = ir_flux - ir;  // Fortran: distance = jr - ir
    // Argument c d, not pi d / nr; 2c / dr applied last (the Fortran premultiplies its table).
    // Must not change: either alternative rounds some weights differently.
    const double x = (d == 0) ? 0.0 : fortran_integer_power(-1.0, d) / (2.0 * std::tan(c * d));
    weights[ir] = ((2.0 * c) * x) / dr;
  }
  return weights;
}

// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr (the wf_deriv sum). The Fortran multiplies by
// D as the complex number (D, 0), which adds only signed zeros: the same values, up to the sign of
// a component that is exactly 0.
std::complex<double> flux_r_derivative(std::span<const double> weights, std::span<const std::complex<double>> column) {
  if (column.size() != weights.size())
    throw std::invalid_argument(
        std::format("flux_r_derivative: {} values for {} weights", column.size(), weights.size()));
  std::complex<double> derivative(0.0, 0.0);
  for (std::size_t ir = 0; ir < weights.size(); ++ir) derivative = derivative + weights[ir] * column[ir];
  return derivative;
}

}  // namespace spodvr
