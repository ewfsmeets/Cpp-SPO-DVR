// On-the-fly S-matrix: accumulation of the time-to-energy transform of the BK coefficients and
// the final S-matrix elements, probabilities and averages.
// Fortran: Analyzor.f90 SetupAnalysis, BalintKurtiCoeffs, smatrix_calc, smatrix_combine
#include "analysis/s_matrix.hpp"

#include <cmath>
#include <span>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

}  // namespace

AnalysisEnergies make_analysis_energies(const Config& config, const GasPhase& gas) {
  const InputParameters& p = config.input;
  const auto& units = units::legacy;
  double start = config.derived.analysis_energy_start;  // Fortran: enrgi
  double step = config.derived.analysis_energy_step;    // Fortran: denrgy
  int count = p.n_energies;                             // Fortran: nenrgy
  if (count == Constants::nen95_sentinel) {
    step = (Constants::pi / p.t_end) * units.energy_factor();
    start = static_cast<int>((p.e_min_wave_packet + p.e_parallel) / step) * step + step;
    count = static_cast<int>((p.e_max_wave_packet + p.e_parallel - start) / step) + 1;
  }

  // Fortran: ensctr = enrgi + dfloat(ien-1) * denrgy is evaluated in quad precision; plain
  // double here (see CLAUDE.md, quad-precision sites).
  AnalysisEnergies e;
  for (int ien = 0; ien < count; ++ien) {
    const double collision = start + ien * step;
    const double total = collision + units.energy_factor() * gas.initial.energy;
    e.collision.push_back(collision);
    e.total.push_back(total);
    e.total_hartree.push_back(total / units.energy_factor());
  }
  return e;
}

SMatrixAccumulator::SMatrixAccumulator(const Config& config, const AnalysisBasis& basis, const GasPhase& gas)
    : config_(config),
      basis_(basis),
      gas_(gas),
      energies_(make_analysis_energies(config, gas)),
      analysis_time_step_(config.derived.steps_per_analysis * config.input.time_step),
      sums_(basis.n_states() * basis.n_v() * energies_.size(), cplx(0.0, 0.0)) {}

// Fortran: BalintKurtiCoeffs, "direct accumulation for some scattering energies"
void SMatrixAccumulator::add(const BkCoefficients& bk, double time) {
  const int n = basis_.n_states() * basis_.n_v();
  for (int ien = 0; ien < energies_.size(); ++ien) {
    const cplx phase_factor = std::exp(Constants::imaginary_unit * energies_.total_hartree[ien] * time);
    const std::span<cplx> sums = std::span(sums_).subspan(n * ien, n);
    for (int k = 0; k < n; ++k) sums[k] = sums[k] + bk.values[k] * phase_factor * analysis_time_step_;
  }
}

ScatteringResults SMatrixAccumulator::results() const {
  ScatteringResults r;
  for (int ien = 0; ien < energies_.size(); ++ien) r.energies.push_back(result_at(ien));
  return r;
}

// Fortran: smatrix_calc (nread = 0, 1) followed by smatrix_combine, for one energy.
EnergyResult SMatrixAccumulator::result_at(int ien) const {
  const InputParameters& p = config_.input;
  const DerivedParameters& d = config_.derived;
  const auto& units = units::legacy;
  const AnalysisBasis& a = basis_;
  const double chaev = units.energy_factor();
  const double mass = d.total_mass;
  const double pi = Constants::pi;
  const int n_states = a.n_states();
  const int n_v = a.n_v();
  const int n_rot = a.n_rotation();
  const int n_orders = a.max_order() + 1;

  EnergyResult e;
  e.total_energy = energies_.total[ien];
  e.collision_energy = energies_.collision[ien];
  e.z_energy = e.collision_energy - a.diffraction_energy(0) * chaev;  // channel 0 is (0, 0): edifin

  // Momentum of the analysis energy perpendicular to the surface (sk0) and the weight of that
  // momentum in the initial wave packet (sumanl).
  const double kx = d.kx;
  const double ky = d.ky;
  const double k_parallel_sq = kx * kx - (2.0 * std::cos(d.system.skewing_angle) * kx * ky) + ky * ky;  // rkxy2
  e.total_momentum = std::sqrt(2.0 * mass * e.collision_energy / chaev);                               // rkt
  const double k0 = std::sqrt(e.total_momentum * e.total_momentum - k_parallel_sq);                     // sk0
  const double weight = std::exp(-((d.kz - k0) * (d.kz - k0)) / (2.0 * d.wave_packet_alpha));          // sumanl
  const double prefactor = std::pow(0.5 * d.wave_packet_alpha / pi, 0.25) / mass;

  // S-matrix elements of the open channels. Fortran: smatrix_calc
  e.s_matrix.assign(n_states * n_v, cplx(0.0, 0.0));
  e.open.assign(n_states * n_v, 0);
  const std::span<const cplx> sums = std::span(sums_).subspan(n_states * n_v * ien, n_states * n_v);
  for (int v = 0; v < n_v; ++v)
    for (int is = 0; is < n_states; ++is) {
      const int k = is + n_states * v;
      const double kz_sq = 2.0 * mass * (e.total_energy / chaev - a.channel_energy(is, v));  // skzsq
      if (!(kz_sq > 0.0)) continue;  // closed channel: S = 0
      const double kz = std::sqrt(kz_sq);
      const double arg = -kz * p.z_analysis + (d.kz - k0) * p.z0 - pi;
      e.s_matrix[k] = sums[k] * prefactor * std::exp(Constants::imaginary_unit * arg) * std::sqrt(kz * k0 / weight);
      e.open[k] = 1;
    }

  // Probabilities, phases, sums and averages, in the regular state order. Fortran: smatrix_combine
  e.probability.assign(n_states * n_v, 0.0);
  e.phase.assign(n_states * n_v, 0.0);
  e.p_vibration.assign(n_v, 0.0);
  e.p_rotation.assign(n_rot * n_v, 0.0);
  e.p_order.assign(n_orders * n_rot * n_v, 0.0);
  e.e_rotation_v.assign(n_v, 0.0);
  e.e_diffraction_v.assign(n_v, 0.0);
  double sum = 0.0;
  double e_vib = 0.0;
  double e_rot = 0.0;
  double e_dif = 0.0;
  const std::vector<int> regular = a.regular_state_order();
  for (int v = 0; v < n_v; ++v)
    for (const int is : regular) {
      const int k = is + n_states * v;
      const cplx s = e.s_matrix[k];
      const double prob = s.real() * s.real() + s.imag() * s.imag();  // Fortran: selmnt * conjg(selmnt)
      const int irotana = a.state_rotation(is);
      const double aerot = gas_.rotational_energies[a.j(is)];
      const double edif = a.diffraction_energy(a.state_diffraction(is));
      e_vib = e_vib + gas_.vibrational_energies[v] * prob;
      e_rot = e_rot + aerot * prob;
      e_dif = e_dif + edif * prob;
      e.e_rotation_v[v] = e.e_rotation_v[v] + aerot * prob;
      e.e_diffraction_v[v] = e.e_diffraction_v[v] + edif * prob;

      e.probability[k] = prob;
      e.phase[k] = s.real() == 0.0 ? 0.0 : std::atan2(s.imag(), s.real());

      e.p_rotation[irotana + n_rot * v] = e.p_rotation[irotana + n_rot * v] + prob;
      const int io = a.channel(is).order + n_orders * (irotana + n_rot * v);
      e.p_order[io] = e.p_order[io] + prob;
      e.p_vibration[v] = e.p_vibration[v] + prob;
      sum = sum + prob;
    }

  for (int v = 0; v < n_v; ++v) {
    if (e.p_vibration[v] != 0.0) {
      e.e_rotation_v[v] = e.e_rotation_v[v] * chaev / e.p_vibration[v];
      e.e_diffraction_v[v] = e.e_diffraction_v[v] * chaev / e.p_vibration[v];
    } else {
      e.e_rotation_v[v] = 0.0;
      e.e_diffraction_v[v] = 0.0;
    }
  }
  e.p_total = sum;
  e.p_reaction = 1.0 - sum;
  e.e_vibration = sum != 0.0 ? e_vib * chaev / sum : 0.0;
  e.e_rotation = sum != 0.0 ? e_rot * chaev / sum : 0.0;
  e.e_diffraction = sum != 0.0 ? e_dif * chaev / sum : 0.0;
  return e;
}

}  // namespace spodvr
