// Tests: the building blocks of the flux analysis (FluxAnalyzor.f90): the wave packet's momentum
// amplitude b(k), the flux energies (flux Q3, flux Q5) and the r-derivative weights D at the flux
// surface (WaveFunction.f90 getWaveFuncDerivAtZr).
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "comparison.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "input/wave_packet_params.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/initial_wave_packet.hpp"

using namespace spodvr;
using spodvr::testing::nan_aware_max;
using cplx = std::complex<double>;

namespace {

constexpr double pi = Constants::pi;

nlohmann::json hardwall_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// hardwall_small (EMinWP 0.2, EMaxWP 0.6 eV, Z0 11.44) with the flux analysis at r = 2.0 bohr.
Config flux_config(double e_min, double e_max, double delta_e, double e_min_wp = 0.2, double e_max_wp = 0.6) {
  nlohmann::json json = hardwall_json();
  json["EMinWP"] = e_min_wp;
  json["EMaxWP"] = e_max_wp;
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 2.0;
  json["EMinFluxAna"] = e_min;
  json["EMaxFluxAna"] = e_max;
  json["DeltaEFluxAna"] = delta_e;
  return build_config(parse_input(json));
}

// b(k) of the hardwall_small wave packet, energies converted as in FluxAnalyzor.f90.
cplx hardwall_amplitude(const Config& config, double k) {
  const InputParameters& p = config.input;
  return wave_packet_momentum_amplitude(k, units::precise.ev_to_hartree(p.e_min_wave_packet),
                                        units::precise.ev_to_hartree(p.e_max_wave_packet), config.derived.total_mass,
                                        p.z0);
}

struct PacketMomenta {
  double kmin, kmax, k0;
};

PacketMomenta packet_momenta(const Config& config) {
  const double mass = config.derived.total_mass;
  const double kmin = std::sqrt(2.0 * mass * units::precise.ev_to_hartree(config.input.e_min_wave_packet));
  const double kmax = std::sqrt(2.0 * mass * units::precise.ev_to_hartree(config.input.e_max_wave_packet));
  return {kmin, kmax, 0.5 * (kmin + kmax)};
}

}  // namespace

// ------------------------------------------------------------------------------------------------
// b(k)
// ------------------------------------------------------------------------------------------------

namespace {

// getBkCoeffWavePacket written out as the Fortran source and its gfortran code compute it, with the
// real and imaginary parts of the exponent as separate products.
cplx fortran_amplitude(double momentum, double low, double high, double mass, double z0) {
  const double kmin = std::sqrt(2.0 * mass * low);
  const double kmax = std::sqrt(2.0 * mass * high);
  const double k0 = (kmin + kmax) / 2.0;
  const double alpha = ((k0 - kmin) * (k0 - kmin)) / 4.0;
  const double delk = k0 - momentum;
  const double xi2 = 1.0 / (4.0 * alpha);
  const double norm = std::pow(2.0 * xi2 / pi, 0.25);
  const cplx arg(-(xi2 * (delk * delk)), delk * z0);  // - xi2 * delk**2 + Constants_I * delk * z0
  const cplx e = std::exp(arg);
  return {norm * e.real(), norm * e.imag()};
}

}  // namespace

// Bitwise the Fortran expression at 401 momenta over kmin - w .. kmax + w (the same values as the
// reference Fortran's getBkCoeffWavePacket, 2026-11-12); the order (xi2 dk) dk, which
// wave_packet_momentum_amplitude rules out, changes some of them.
TEST(FluxMomentumAmplitude, IsBitwiseTheFortranExpression) {
  const Config config = flux_config(0.2, 0.6, 0.1);
  const InputParameters& p = config.input;
  const double low = units::precise.ev_to_hartree(p.e_min_wave_packet);
  const double high = units::precise.ev_to_hartree(p.e_max_wave_packet);
  const double mass = config.derived.total_mass;
  const auto [kmin, kmax, k0] = packet_momenta(config);
  const double w = k0 - kmin;

  int other_order = 0;
  for (int i = 0; i <= 400; ++i) {
    const double k = kmin - w + i * (kmax - kmin + 2.0 * w) / 400.0;
    const cplx b = wave_packet_momentum_amplitude(k, low, high, mass, p.z0);
    const cplx expected = fortran_amplitude(k, low, high, mass, p.z0);
    EXPECT_EQ(b.real(), expected.real()) << "k = " << k;
    EXPECT_EQ(b.imag(), expected.imag()) << "k = " << k;

    const double dk = k0 - k;
    const double xi2 = 1.0 / (4.0 * (((k0 - kmin) * (k0 - kmin)) / 4.0));
    const cplx swapped = std::pow(2.0 * xi2 / pi, 0.25) * std::exp(cplx(-((xi2 * dk) * dk), dk * p.z0));
    if (swapped != b) ++other_order;
  }
  EXPECT_GT(other_order, 0);
}

// b(k) = (2 / (pi w^2))^(1/4) exp(-((k0 - k) / w)^2) exp(i (k0 - k) z0) with w = k0 - kmin, written
// out independently of the Fortran's xi2 and alpha. The two forms are equal in exact arithmetic and
// round differently by a few ulp; tolerance 1e-14 of the peak (about 90 ulp; measured 2026-11-12:
// 1.1e-16 absolute, the peak 0.54).
TEST(FluxMomentumAmplitude, MatchesTheClosedForm) {
  const Config config = flux_config(0.2, 0.6, 0.1);
  const auto [kmin, kmax, k0] = packet_momenta(config);
  const double w = k0 - kmin;
  const double z0 = config.input.z0;
  const double peak = std::pow(2.0 / (pi * w * w), 0.25);

  double worst = 0.0;
  for (int i = 0; i <= 200; ++i) {
    const double k = kmin - w + i * (kmax - kmin + 2.0 * w) / 200.0;
    const double dk = k0 - k;
    const cplx expected = peak * std::exp(-(dk / w) * (dk / w)) * std::polar(1.0, dk * z0);
    worst = nan_aware_max(worst, std::abs(hardwall_amplitude(config, k) - expected));
  }
  EXPECT_LE(worst, 1e-14 * peak);
}

// The integral of |b|^2 over k is 1, and |b|^2 at kmin and kmax is e^-2 of its peak at k0.
TEST(FluxMomentumAmplitude, IsNormalisedWithItsEdgesAtEMinusTwo) {
  const Config config = flux_config(0.2, 0.6, 0.1);
  const auto [kmin, kmax, k0] = packet_momenta(config);
  const double w = k0 - kmin;
  const double peak = std::norm(hardwall_amplitude(config, k0));

  // Trapezoid rule over k0 +- 10 w (the tails beyond are below e^-200; for a Gaussian the rule is
  // exact to far below rounding). The tolerance covers the rounding of 4001 terms (random walk
  // sqrt(4001) * 1.1e-16 = 7e-15; measured 2026-11-12: 2.9e-15).
  const int n = 4000;
  const double dk = 20.0 * w / n;
  double integral = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double weight = (i == 0 || i == n) ? 0.5 : 1.0;
    integral += weight * std::norm(hardwall_amplitude(config, k0 - 10.0 * w + i * dk));
  }
  EXPECT_NEAR(integral * dk, 1.0, 1e-13);
  EXPECT_NEAR(std::norm(hardwall_amplitude(config, kmin)), std::exp(-2.0) * peak, 1e-14 * peak);
  EXPECT_NEAR(std::norm(hardwall_amplitude(config, kmax)), std::exp(-2.0) * peak, 1e-14 * peak);
}

// |b(k)|^2 is the momentum distribution of the C++ initial wave packet: its Z part g(Z) gives
// g(k) = dz / sqrt(2 pi) sum_Z g(Z) exp(i k Z) = conj b(k) (flux Q9: b has the opposite phase).
// They differ only through the eV -> hartree factors (flux Q5: b precise, the packet legacy; the
// energies differ by 27.2114 * 0.036749308286 - 1 = 1.27e-7, so all packet momenta scale by 1 + eps,
// eps = 6.37e-8). With u = (k0 - k) / w the |b|^2 difference is
// max_u e^(-2u^2) |1 + 4u (k0/w - u)| eps of the peak = 4.61 eps = 2.94e-7 (k0/w = 3.73, u = 0.43)
// and the phase difference z0 k0 eps = 7.3e-6 rad; measured 2026-11-12: 2.93e-7 and 7.3e-6 rad.
// The tolerances are 10 and 14 times these. This fixes the 1 / (2 pi) normalisation of c_n.
TEST(FluxMomentumAmplitude, IsTheMomentumDistributionOfTheInitialWavePacket) {
  const Config config = flux_config(0.2, 0.6, 0.1);
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const PotentialGrid potential = compute_potential_grid(*make_potential(config.input), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());
  const InitialWavePacket packet = make_initial_wave_packet(config, grids, gas, fft);

  // g(Z) = psi(Z, r) / phi(r) at the r point where phi is largest.
  const auto phi = gas.states.wave_function(config.input.v_in, config.input.j_in);
  int ir = 0;
  for (int i = 0; i < config.input.nr; ++i)
    if (std::abs(phi[i]) > std::abs(phi[ir])) ir = i;
  const UniformGrid& z = grids.z_specular;
  auto spectrum = [&](double k) {
    cplx sum(0.0, 0.0);
    for (int iz = 0; iz < z.size(); ++iz) sum += (packet.specular(ir, iz) / phi[ir]) * std::polar(1.0, k * z[iz]);
    return sum * (config.input.dz / std::sqrt(2.0 * pi));
  };

  const auto [kmin, kmax, k0] = packet_momenta(config);
  const double peak = std::norm(hardwall_amplitude(config, k0));
  double worst_density = 0.0;
  double worst_phase = 0.0;
  for (int i = 0; i <= 20; ++i) {
    const double k = kmin + i * (kmax - kmin) / 20.0;
    const cplx b = hardwall_amplitude(config, k);
    const cplx g = spectrum(k);
    worst_density = nan_aware_max(worst_density, std::abs(std::norm(g) - std::norm(b)) / peak);
    worst_phase = nan_aware_max(worst_phase, std::abs(std::arg(b * g)));
  }
  EXPECT_LE(worst_density, 3e-6) << "|g|^2 - |b|^2 relative to the peak";
  EXPECT_LE(worst_phase, 1e-4) << "phase of b g (0 for b = conj g)";
  EXPECT_GT(worst_density, 3e-8) << "flux Q5 should be visible: b uses Constants_eV2Ha, the packet chaev";
}

// ------------------------------------------------------------------------------------------------
// The flux energies
// ------------------------------------------------------------------------------------------------

// INT((E_max - E_min) / dE) + 1 in hartree (flux Q3): for 0.2 - 0.6 eV by 0.1 the quotient rounds
// to just below 4, so 0.6 eV is dropped; by 0.08 it reaches 5; 0.15 does not divide the range; the
// thesis DZ range 0.40 - 0.95 eV by 0.025 gets 22 energies, the last 0.925 eV. The counts are those
// of the reference Fortran's setupFluxAnalyzor (4, 6, 3, 22; a driver run, 2026-11-12).
TEST(FluxEnergies, CountDropsTheEndPointAsTheFortran) {
  struct Case {
    double e_min, e_max, delta_e;
    int count;
    double last_ev;
  };
  for (const Case& c : {Case{0.2, 0.6, 0.1, 4, 0.5}, Case{0.2, 0.6, 0.08, 6, 0.6}, Case{0.2, 0.6, 0.15, 3, 0.5},
                        Case{0.40, 0.95, 0.025, 22, 0.925}}) {
    SCOPED_TRACE(::testing::Message() << c.e_min << " .. " << c.e_max << " by " << c.delta_e);
    const Config config = flux_config(c.e_min, c.e_max, c.delta_e, c.e_min, c.e_max);
    EXPECT_EQ(flux_energy_count(config.input), c.count);
    const FluxEnergies e = make_flux_energies(config, 0.0, 0.0);
    ASSERT_EQ(e.size(), c.count);
    EXPECT_NEAR(units::precise.hartree_to_ev(e.z_energy.back()), c.last_ev, 1e-12);
    const bool end_point_missing = c.last_ev < c.e_max - 1e-9;
    EXPECT_EQ(e.warnings.size(), end_point_missing ? 1U : 0U);
  }
}

TEST(FluxEnergies, NoFluxAnalysisMeansNoEnergies) {
  InputParameters p;
  p.do_flux_analysis = false;
  EXPECT_EQ(flux_energy_count(p), 0);
}

// A DeltaEFluxAna that passes the input checks (> 0 within 1e-12 eV) can still give an absurd
// number of energies, beyond int below about 2e-10 eV for this 0.4 eV range: an input error naming
// DeltaEFluxAna, not an undefined conversion to int.
TEST(FluxEnergies, CountAboveTheLimitIsAnInputError) {
  InputParameters p;
  p.do_flux_analysis = true;
  p.e_min_flux_analysis = 0.2;
  p.e_max_flux_analysis = 0.6;
  for (const double delta_e : {0.4 / (2.0 * max_flux_energies), 1e-11}) {
    SCOPED_TRACE(delta_e);
    p.delta_e_flux_analysis = delta_e;
    try {
      (void)flux_energy_count(p);
      ADD_FAILURE() << "no InputError";
    } catch (const InputError& e) {
      EXPECT_EQ(e.parameter(), "DeltaEFluxAna");
    }
  }
  p.delta_e_flux_analysis = 0.4 / (max_flux_energies - 2);
  EXPECT_LE(flux_energy_count(p), max_flux_energies);
}

// E_n = E_min + n dE with Constants_eV2Ha (flux Q5; the legacy factor would differ), E_tot =
// (E_n + parallel) + rovibrational, k_n = sqrt(2 M E_n), b(k_n) of the packet's precise energies.
TEST(FluxEnergies, FollowTheFortranExpressions) {
  const Config config = flux_config(0.2, 0.6, 0.08);
  const double parallel = 0.0123;
  const double rovibrational = 0.0098765;
  const FluxEnergies e = make_flux_energies(config, parallel, rovibrational);
  const double e_min = 0.2 * 0.036749308286;
  const double step = 0.08 * 0.036749308286;
  const double mass = config.derived.total_mass;
  ASSERT_EQ(e.size(), 6);
  for (int n = 0; n < e.size(); ++n) {
    const double z_energy = e_min + n * step;
    EXPECT_EQ(e.z_energy[n], z_energy);
    EXPECT_EQ(e.total[n], (z_energy + parallel) + rovibrational);
    EXPECT_EQ(e.momentum[n], std::sqrt(2.0 * mass * z_energy));
    EXPECT_EQ(e.amplitude[n], hardwall_amplitude(config, e.momentum[n]));
  }
  EXPECT_NE(e.z_energy[0], units::legacy.ev_to_hartree(0.2));
}

// ------------------------------------------------------------------------------------------------
// The derivative weights D
// ------------------------------------------------------------------------------------------------

namespace {

// The Fortran table, written out with an explicit sign instead of (-1)**d.
double fortran_weight(int nr, double dr, int ir_flux, int ir) {
  const int d = ir_flux - ir;
  if (d == 0) return 0.0;
  const double c = pi / nr;
  const double sign = (d % 2 == 0) ? 1.0 : -1.0;
  const double x = sign / (2.0 * std::tan(c * d));
  return ((2.0 * c) * x) / dr;
}

}  // namespace

// Bitwise the Fortran expression; the two evaluation orders the comment in flux_derivative_weights
// rules out (argument pi d / nr, and 2c / dr applied first) change some weights.
TEST(FluxDerivativeWeights, AreBitwiseTheFortranTable) {
  const int nr = 56;
  const double dr = 0.15;
  const int ir_flux = 8;
  const std::vector<double> w = flux_derivative_weights(nr, dr, ir_flux);
  ASSERT_EQ(static_cast<int>(w.size()), nr);
  EXPECT_EQ(w[ir_flux], 0.0);
  int other_argument = 0;
  int other_scaling = 0;
  for (int ir = 0; ir < nr; ++ir) {
    EXPECT_EQ(w[ir], fortran_weight(nr, dr, ir_flux, ir)) << "ir = " << ir;
    const int d = ir_flux - ir;
    if (d == 0) continue;
    const double sign = (d % 2 == 0) ? 1.0 : -1.0;
    if (((2.0 * (pi / nr)) * (sign / (2.0 * std::tan((pi * d) / nr)))) / dr != w[ir]) ++other_argument;
    if (((2.0 * (pi / nr)) / dr) * (sign / (2.0 * std::tan((pi / nr) * d))) != w[ir]) ++other_scaling;
  }
  EXPECT_GT(other_argument, 0);
  EXPECT_GT(other_scaling, 0);
}

// D(r_c + m dr) = -D(r_c - m dr): the derivative matrix of a periodic grid is antisymmetric.
TEST(FluxDerivativeWeights, AreAntisymmetricAboutTheFluxSurface) {
  const int nr = 56;
  const int ir_flux = 20;
  const std::vector<double> w = flux_derivative_weights(nr, 0.15, ir_flux);
  for (int m = 1; ir_flux + m < nr && ir_flux - m >= 0; ++m) EXPECT_EQ(w[ir_flux + m], -w[ir_flux - m]) << m;
}

// D annihilates a constant and differentiates exp(2 pi i m r / L), L = nr dr, exactly for
// |m| < nr / 2 (the Fourier-grid derivative), to rounding. Rounding scale: one unit roundoff on
// every term is 1.1e-16 sum |D| (sum |D| = 46.1 / bohr), i.e. 6.8e-15 of k for m = 1. Tolerances: 1e-14
// sum |D| for the constant (the weights of d and nr - d cancel up to the rounding of tan; measured
// 2026-11-12: 6.9e-16) and 2e-13 for the waves (about 30 roundoffs; measured 2.1e-14).
TEST(FluxDerivativeWeights, DifferentiatePeriodicWavesExactly) {
  const int nr = 56;
  const double dr = 0.15;
  const double r_start = 0.8;
  const double length = nr * dr;
  for (const int ir_flux : {0, 8, 27, 55}) {
    SCOPED_TRACE(ir_flux);
    const std::vector<double> w = flux_derivative_weights(nr, dr, ir_flux);
    double scale = 0.0;
    for (const double x : w) scale += std::abs(x);

    const std::vector<cplx> constant(nr, cplx(1.0, -2.0));
    EXPECT_LE(std::abs(flux_r_derivative(w, constant)), 1e-14 * scale);

    double worst = 0.0;
    for (const int m : {1, -3, 7, 27}) {
      const double k = 2.0 * pi * m / length;
      std::vector<cplx> wave(nr);
      for (int ir = 0; ir < nr; ++ir) wave[ir] = std::polar(1.0, k * (r_start + ir * dr));
      const cplx expected = cplx(0.0, k) * wave[ir_flux];
      worst = nan_aware_max(worst, std::abs(flux_r_derivative(w, wave) - expected) / std::abs(k));
    }
    EXPECT_LE(worst, 2e-13);
  }
}

// The sum starts at 0 and adds D(ir) psi(ir) for ir = 0, 1, ..., nr - 1, the zero term D(ir_flux)
// psi(ir_flux) included: bitwise that order; the reverse order rounds differently on this column;
// a non-finite psi(r_c) makes the result NaN through 0 * inf, as in the Fortran.
TEST(FluxRDerivative, SumsFromZeroInAscendingOrderWithTheZeroTerm) {
  const int nr = 56;
  const int ir_flux = 8;
  const std::vector<double> w = flux_derivative_weights(nr, 0.15, ir_flux);
  std::vector<cplx> column(nr);
  for (int ir = 0; ir < nr; ++ir) column[ir] = std::polar(1.0 + 0.5 * std::sin(1.3 * ir), 0.37 * ir);

  cplx ascending(0.0, 0.0);
  for (int ir = 0; ir < nr; ++ir) ascending = ascending + w[ir] * column[ir];
  cplx descending(0.0, 0.0);
  for (int ir = nr - 1; ir >= 0; --ir) descending = descending + w[ir] * column[ir];

  const cplx result = flux_r_derivative(w, column);
  EXPECT_EQ(result.real(), ascending.real());
  EXPECT_EQ(result.imag(), ascending.imag());
  EXPECT_NE(descending, ascending) << "this column does not distinguish the summation orders";

  column[ir_flux] = cplx(std::numeric_limits<double>::infinity(), 0.0);
  EXPECT_TRUE(std::isnan(flux_r_derivative(w, column).real()));
}

TEST(FluxRDerivative, RejectsAColumnOfAnotherLength) {
  const std::vector<double> w = flux_derivative_weights(8, 0.1, 3);
  const std::vector<cplx> column(7);
  EXPECT_THROW((void)flux_r_derivative(w, column), std::invalid_argument);
}
