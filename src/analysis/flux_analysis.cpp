// The flux analysis: its energies, the r-derivative at the flux surface, the time-to-energy sums of
// psi and d psi / dr there, and the reaction probability. Always strict floating point (spodvr_strict_objects).
// Fortran: FluxAnalyzor.f90; WaveFunction.f90 getWaveFuncDerivAtZr
#include "analysis/flux_analysis.hpp"

#include <algorithm>
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

using cplx = std::complex<double>;

// FluxAnalyzor.f90 converts its energies and the wave packet's with Constants_eV2Ha (flux Q5).
constexpr const UnitConversion& flux_units = units::precise;

// The last energy "reaches" EMaxFluxAna within this fraction of a step: far above the rounding of
// E_min + n dE (~1e-16 relative), far below a step.
constexpr double end_point_tolerance = 1e-6;

// Blocks per tile of the flush: phi and phi' of one tile (2 x 16 KiB) stay in cache while every
// filled slot is added to them, so they are read once per flush, not once per slot.
constexpr int flush_tile_blocks = 1024;

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
// storeTimeIndependentWF (kZ_initial, bkz, prefactor)
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
    const cplx amplitude = wave_packet_momentum_amplitude(momentum, packet_low, packet_high, mass, p.z0);
    e.z_energy.push_back(z_energy);
    e.total.push_back((z_energy + parallel_energy) + initial_state_energy);
    e.momentum.push_back(momentum);
    e.amplitude.push_back(amplitude);
    e.transform_factor.push_back(flux_transform_factor(momentum, mass, amplitude));
  }
  if (count > 0) {
    std::string warning = end_point_warning(p, e);
    if (!warning.empty()) e.warnings.push_back(std::move(warning));
  }
  return e;
}

// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF, prefactor = kZ_initial / (Constants_2PI * rmass * bkz)
cplx flux_transform_factor(double momentum, double mass, cplx amplitude) {
  const double two_pi_mass = Constants::two_pi * mass;
  // A real over a complex number is a complex division (k, 0) / ((2 pi M) b): Smith's algorithm under
  // -fcx-fortran-rules, as gfortran's. Must stay a division: conj(b) k / |b|^2 rounds differently.
  return momentum / (two_pi_mass * amplitude);
}

// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF,
// coeff = prefactor * EXP(Constants_I * totalScatEnergy * time) * dtAnalysis
cplx flux_step_coefficient(const FluxEnergies& energies, int n, double time, double time_step) {
  // std::exp of the complex argument (i E) t is glibc cexp, as the Fortran's EXP; never std::polar.
  const cplx phase = std::exp(Constants::imaginary_unit * energies.total[n] * time);
  // (A exp) dt in the Fortran's order. gfortran multiplies by the reals E, t, dt (and 2 pi M in A) as
  // by (x, 0): the same values, up to the sign of an exact zero (fortran_deviations.md, signed zeros).
  return (energies.transform_factor[n] * phase) * time_step;
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
// D as the complex number (D, 0): the same values, up to the sign of an exact zero
// (fortran_deviations.md, signed zeros).
cplx flux_r_derivative(std::span<const double> weights, std::span<const cplx> column) {
  if (column.size() != weights.size())
    throw std::invalid_argument(
        std::format("flux_r_derivative: {} values for {} weights", column.size(), weights.size()));
  cplx derivative(0.0, 0.0);
  for (std::size_t ir = 0; ir < weights.size(); ++ir) derivative = derivative + weights[ir] * column[ir];
  return derivative;
}

FluxSurfaceWave::FluxSurfaceWave(int n_blocks, int n_planes)
    : index(n_blocks, n_planes), psi(index.size()), derivative(index.size()) {
  // Zero, written statically over iz: the partition of the flush (first touch).
#pragma omp parallel for schedule(static)
  for (int iz = 0; iz < n_planes; ++iz) {
    std::ranges::fill(std::span(psi).subspan(index(0, iz), index.n0), cplx(0.0, 0.0));
    std::ranges::fill(std::span(derivative).subspan(index(0, iz), index.n0), cplx(0.0, 0.0));
  }
}

FluxCoefficients::FluxCoefficients(int n_blocks, int n_planes, int n_energies)
    : index(n_blocks, n_planes, n_energies), phi(index.size()), phi_derivative(index.size()) {
  // Zero, written statically over iz for every energy: the partition of the flush (first touch).
#pragma omp parallel for schedule(static)
  for (int iz = 0; iz < n_planes; ++iz)
    for (int n = 0; n < n_energies; ++n) {
      std::ranges::fill(std::span(phi).subspan(index(0, iz, n), index.n0), cplx(0.0, 0.0));
      std::ranges::fill(std::span(phi_derivative).subspan(index(0, iz, n), index.n0), cplx(0.0, 0.0));
    }
}

FluxGrid make_flux_grid(const Config& config, const SpatialGrids& grids, int n_rot) {
  return FluxGrid{.n_rot = n_rot,
                  .nx = grids.x.size(),
                  .ny = grids.y.size(),
                  .dz = grids.z_scattering.spacing(),
                  .dx = grids.x.spacing(),
                  .dy = grids.y.spacing(),
                  .skewing_angle = config.derived.system.skewing_angle,
                  .total_mass = config.derived.total_mass,
                  .reduced_mass = config.derived.reduced_mass};
}

namespace {

// (ix, iy) of the block (irot, ixy) of psi, ixy = iy + ny ix. Fortran: getOuterIndex, getInnerIndex
struct GridXY {
  int ix;
  int iy;
};
GridXY block_xy(int block, const FluxGrid& grid) {
  const int ixy = block / grid.n_rot;
  return {ixy / grid.ny, ixy % grid.ny};
}

// P(x, y; E_n) and P(E_n) of energy n into r. Fortran: computeFluxProbability, one ien
void add_energy_results(const FluxCoefficients& c, const FluxEnergies& e, const FluxGrid& grid, int n, FluxResults& r) {
  // ((2 pi M) / (k mu)) dz in the Fortran's order (its prefac), applied to every term as there.
  const double prefactor = ((Constants::two_pi * grid.total_mass) / (e.momentum[n] * grid.reduced_mass)) * grid.dz;
  // iz outer, the blocks in memory order inner: the Fortran's summation order per (x, y).
  for (int iz = 0; iz < c.n_planes(); ++iz)
    for (int block = 0; block < c.n_blocks(); ++block) {
      const auto [ix, iy] = block_xy(block, grid);
      const std::size_t k = c.index(block, iz, n);
      const double flux = std::imag(std::conj(c.phi[k]) * c.phi_derivative[k]);  // Fortran: AIMAG(contrb)
      double& density = r.density[r.density_index(ix, iy, n)];
      density = density + prefactor * flux;
    }

  double sum = 0.0;
  for (int iy = 0; iy < grid.ny; ++iy)
    for (int ix = 0; ix < grid.nx; ++ix) sum = sum + r.density[r.density_index(ix, iy, n)];
  // ((sum dx) dy) sin(angle) in the Fortran's order; not xy_volume_element, which is (dx dy) sin(angle).
  r.probability[n] = ((sum * grid.dx) * grid.dy) * std::sin(grid.skewing_angle);
}

}  // namespace

// Fortran: FluxAnalyzor.f90 computeFluxProbability
FluxResults flux_results(const FluxCoefficients& coefficients, const FluxEnergies& energies, const FluxGrid& grid) {
  if (coefficients.n_blocks() != grid.n_rot * grid.nx * grid.ny || coefficients.n_energies() != energies.size())
    throw std::invalid_argument(
        std::format("flux_results: coefficients of {} blocks and {} energies for {} x {} x {} blocks and {} energies",
                    coefficients.n_blocks(), coefficients.n_energies(), grid.n_rot, grid.nx, grid.ny, energies.size()));
  const int n_energies = energies.size();
  FluxResults r;
  for (const double z_energy : energies.z_energy) r.energy.push_back(flux_units.hartree_to_ev(z_energy));
  r.probability.assign(n_energies, 0.0);
  r.density_index = Index3D(grid.nx, grid.ny, n_energies);
  r.density.assign(r.density_index.size(), 0.0);
  // One energy per iteration, independent of the others (the Fortran's parallel loop over ien);
  // static: every energy costs the same.
#pragma omp parallel for schedule(static)
  for (int n = 0; n < n_energies; ++n) add_energy_results(coefficients, energies, grid, n, r);
  return r;
}

namespace {

// slot_count, checked before the slots and their coefficients are allocated.
int checked_slot_count(int slot_count) {
  if (slot_count < 1)
    throw std::invalid_argument(std::format("FluxAccumulator: {} slots, at least 1 needed", slot_count));
  return slot_count;
}

}  // namespace

FluxAccumulator::FluxAccumulator(int n_blocks, int n_planes, FluxEnergies energies, double time_step, int slot_count)
    : energies_(std::move(energies)),
      time_step_(time_step),
      coefficients_(n_blocks, n_planes, energies_.size()),
      slot_coefficient_index_(energies_.size(), checked_slot_count(slot_count)),
      slot_coefficients_(slot_coefficient_index_.size(), cplx(0.0, 0.0)) {
  slots_.reserve(slot_count);
  for (int k = 0; k < slot_count; ++k) slots_.emplace_back(n_blocks, n_planes);
}

// Fortran: FluxAnalyzor.f90 storeFluxCoeffs (the coeff of each energy and the wave function at r_c)
void FluxAccumulator::add_step(double time, const std::function<void(FluxSurfaceWave&)>& fill) {
  fill(slots_[pending_]);
  for (int n = 0; n < energies_.size(); ++n)
    slot_coefficients_[slot_coefficient_index_(n, pending_)] = flux_step_coefficient(energies_, n, time, time_step_);
  pending_ = pending_ + 1;
  if (pending_ == slot_count()) flush();
}

// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF (the sums). Each element gets the filled slots
// in their time order, ((phi + c_1 psi_1) + c_2 psi_2) + ..., one addition per step as in the Fortran:
// so the result is bitwise the same for any slot count and any thread count.
void FluxAccumulator::flush() {
  const int n_blocks = coefficients_.n_blocks();
  const int n_planes = coefficients_.n_planes();
  const int n_energies = coefficients_.n_energies();
  // Static over iz: the partition that zero-filled the coefficients and the slots (first touch).
#pragma omp parallel for schedule(static)
  for (int iz = 0; iz < n_planes; ++iz)
    for (int first_block = 0; first_block < n_blocks; first_block += flush_tile_blocks)
      for (int n = 0; n < n_energies; ++n)
        add_tile(iz, first_block, std::min(flush_tile_blocks, n_blocks - first_block), n);
  pending_ = 0;
}

// phi += c_n psi and phi' += c_n d psi / dr for tile_size blocks from first_block, on plane iz for
// energy n, slot by slot.
void FluxAccumulator::add_tile(int iz, int first_block, int tile_size, int n) {
  const std::size_t start = coefficients_.index(first_block, iz, n);
  const std::span<cplx> phi = std::span(coefficients_.phi).subspan(start, tile_size);
  const std::span<cplx> phi_derivative = std::span(coefficients_.phi_derivative).subspan(start, tile_size);
  for (int k = 0; k < pending_; ++k) {
    const cplx coefficient = slot_coefficients_[slot_coefficient_index_(n, k)];
    const FluxSurfaceWave& slot = slots_[k];
    const std::span<const cplx> psi = std::span(slot.psi).subspan(slot.index(first_block, iz), tile_size);
    const std::span<const cplx> derivative = std::span(slot.derivative).subspan(slot.index(first_block, iz), tile_size);
    for (int b = 0; b < tile_size; ++b) {
      phi[b] = phi[b] + coefficient * psi[b];  // Fortran: phi = phi + coeff * wfAtZr
      phi_derivative[b] = phi_derivative[b] + coefficient * derivative[b];
    }
  }
}

const FluxCoefficients& FluxAccumulator::coefficients() const {
  if (pending_ > 0)
    throw std::logic_error(
        std::format("FluxAccumulator: {} analysis steps not yet added to the coefficients; "
                    "flush() first",
                    pending_));
  return coefficients_;
}

}  // namespace spodvr
