// Fine analysis (Fortran nread = 2): S-matrix on the spectral grid of the BK time series and
// its interpolation to the analysis energies.
// Fortran: AsympInfo.f90 ReadOrWriteAsympInfo, Analyzor.f90 gtaife, getbkc, hlftrf, shbkar,
// probaf, StateToStateProbs, smatrix_calc (nread = 2)
#include "analysis/fine_analysis.hpp"

#include <algorithm>
#include <cmath>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fft_sizes.hpp"
#include "io/binary_file.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

// Number of series transformed together by one batched FFT.
constexpr int fft_batch = 64;

}  // namespace

// Fortran: ReadOrWriteAsympInfo, nread = 2
FineAnalysisGrid make_fine_analysis_grid(const Config& config, int n_records_in_file) {
  const double chaev = units::legacy.energy_factor();
  FineAnalysisGrid g{};
  g.n_records_in_file = n_records_in_file;
  g.time_step = config.derived.analysis_time_step;

  // Records up to tEnd. Fortran: ndatus = nint(tottym / deltat); more than in the file is an error.
  const int n_wanted = static_cast<int>(std::lround(config.input.t_end / g.time_step));
  if (n_wanted > n_records_in_file)
    throw InputError("tEnd", "tEnd needs " + std::to_string(n_wanted) + " BK records but the file has only " +
                                 std::to_string(n_records_in_file));
  g.n_records_available = n_wanted;
  if (g.n_records_available < 2)
    throw InputError("tEnd", "the fine analysis needs at least 2 BK records (tEnd / dtAna >= 2)");

  g.n_records_used = lowest_fft_size(g.n_records_available);
  g.transform_length = 2 * g.n_records_used;
  g.effective_time = g.n_records_used * g.time_step;
  g.spectral_resolution = Constants::pi / g.effective_time * chaev;

  // Grid points k = 1 .. ndused with E_z = k spres >= 0 (Fortran: imeans = 1, nwould = ndused;
  // echeck = 0 because the energy zero is the initial channel energy).
  const int lowest = 1;
  const int highest = g.n_records_used;
  g.meaningful_min = lowest * g.spectral_resolution;
  g.meaningful_max = highest * g.spectral_resolution;

  // Keep the points inside the 95% confidence interval [ectmin, ectmax] (strict comparisons).
  g.last_index = highest;
  for (int k = lowest; k <= highest; ++k)
    if (config.derived.confidence_e_max < k * g.spectral_resolution) {
      g.last_index = k - 1;
      break;
    }
  g.first_index = lowest;
  for (int k = highest; k >= lowest; --k)
    if (config.derived.confidence_e_min > k * g.spectral_resolution) {
      g.first_index = k + 1;
      break;
    }
  return g;
}

std::vector<double> FineAnalysisGrid::z_energies() const {
  std::vector<double> e;
  for (int i = 0; i < size(); ++i) e.push_back(z_energy(i));
  return e;
}

BkSeries::BkSeries(const FineAnalysisGrid& grid, const AnalysisBasis& basis, double initial_channel_energy)
    : n_used_(grid.n_records_used),
      n_series_(basis.n_states() * basis.n_v()),
      time_step_(grid.time_step),
      energy_(initial_channel_energy),
      values_(static_cast<std::size_t>(n_series_) * (n_used_ + 1)) {
  for (int k = 0; k < n_series_; ++k) values_[static_cast<std::size_t>(k) * (n_used_ + 1)] = cplx(0.0, 0.0);
}

// Fortran: getbkc, bkar(idat, is, iv) = balcf(is, iv) * exp(ima * ev0jl * idat * deltat)
void BkSeries::add_record(double time, const BkCoefficients& bk) {
  const int j = n_added_ + 1;
  n_added_ = j;
  if (j > n_used_) return;

  // The record times of native files are t_end + (step - n_steps) dt, equal to j deltat up to
  // rounding; the transform uses j deltat (as Fortran).
  const double expected = j * time_step_;
  if (std::abs(time - expected) > 1e-9 * std::max(1.0, expected))
    throw CheckpointError("BK record " + std::to_string(j) + " is at t = " + std::to_string(time) + ", expected " +
                          std::to_string(expected) + " (the fine analysis needs records every dtAna from t = 0)");

  const cplx shift = std::exp(cplx(0.0, energy_ * j * time_step_));
  for (int k = 0; k < n_series_; ++k) values_[static_cast<std::size_t>(k) * (n_used_ + 1) + j] = bk.values[k] * shift;
}

void BkSeries::check_complete() const {
  if (n_added_ < n_used_)
    throw CheckpointError("the BK file has " + std::to_string(n_added_) + " records, the fine analysis needs " +
                          std::to_string(n_used_));
}

std::span<const cplx> BkSeries::series(int k) const {
  return std::span<const cplx>(values_).subspan(static_cast<std::size_t>(k) * (n_used_ + 1), n_used_ + 1);
}

// Fortran: hlftrf + shbkar (here as one zero-padded transform with trapezoid weights)
std::vector<cplx> half_fourier_transform(const BkSeries& series, const FineAnalysisGrid& grid, const Fft& fft) {
  const int n = grid.n_records_used;
  const int length = grid.transform_length;
  const int n_points = grid.size();
  const int n_series = series.n_series();
  const int batch = std::min(fft_batch, std::max(n_series, 1));
  const auto plan = fft.plan_complex(FftLayout{length, batch, 1, length}, FftSign::Plus, 1.0);

  std::vector<cplx> amplitudes(static_cast<std::size_t>(n_series) * n_points);
  std::vector<cplx> work(static_cast<std::size_t>(length) * batch);
  for (int first = 0; first < n_series; first += batch) {
    // Weighted series, zero-padded to 2 ndused (unused columns of the last batch stay zero).
    std::fill(work.begin(), work.end(), cplx(0.0, 0.0));
    for (int b = 0; b < batch && first + b < n_series; ++b) {
      const std::span<const cplx> c = series.series(first + b);
      const std::span<cplx> w = std::span<cplx>(work).subspan(static_cast<std::size_t>(b) * length, length);
      w[0] = 0.5 * c[0];
      for (int j = 1; j < n; ++j) w[j] = c[j];
      w[n] = 0.5 * c[n];
    }
    plan->execute(work);

    // Output index k is E_k = k pi / T; keep the grid points first_index .. last_index.
    for (int b = 0; b < batch && first + b < n_series; ++b)
      for (int i = 0; i < n_points; ++i)
        amplitudes[(first + b) + static_cast<std::size_t>(n_series) * i] =
            grid.time_step * work[static_cast<std::size_t>(b) * length + grid.first_index + i];
  }
  return amplitudes;
}

// Fortran: probaf (energies, channel momenta, open channels, prefactor), then the sums of
// smatrix_combine (identical to those of probaf)
ScatteringResults fine_spectrum(const Config& config, std::span<const double> z_energies, const AnalysisBasis& basis,
                                const GasPhase& gas, std::span<const cplx> amplitudes) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const double chaev = units::legacy.energy_factor();
  const double mass = d.total_mass;
  const double pi = Constants::pi;
  const int n_states = basis.n_states();
  const int n_v = basis.n_v();
  const double initial_energy = basis.channel_energy(basis.initial_state(), d.v_in_index - 1);  // ev0jl
  const double prefactor = std::pow(0.5 * d.wave_packet_alpha / pi, 0.25) / mass;               // prefac

  ScatteringResults spectrum;
  for (std::size_t i = 0; i < z_energies.size(); ++i) {
    EnergyResult e{};
    e.z_energy = z_energies[i];  // fetrz = fetot (eV)
    e.total_energy = e.z_energy + initial_energy * chaev;
    e.collision_energy = e.z_energy + basis.diffraction_energy(0) * chaev;
    e.total_momentum = std::sqrt(2.0 * mass * e.collision_energy / chaev);

    // Momentum of the elastic channel (fsk0) and its weight in the initial wave packet (fsuman).
    const double k0 = std::sqrt(2.0 * mass * e.z_energy / chaev);
    const double weight = std::exp(-((d.kz - k0) * (d.kz - k0)) / (2.0 * d.wave_packet_alpha));

    e.s_matrix.assign(n_states * n_v, cplx(0.0, 0.0));
    e.open.assign(n_states * n_v, 0);
    for (int v = 0; v < n_v; ++v)
      for (int is = 0; is < n_states; ++is) {
        const int k = is + n_states * v;
        // Final momentum in Z (fskz). Open if kz^2 >= 0 (probaf; the on-the-fly analysis uses > 0).
        const double kz_sq = 2.0 * mass * (e.z_energy / chaev - (basis.channel_energy(is, v) - initial_energy));
        if (!(kz_sq >= 0.0)) continue;
        const double kz = std::sqrt(kz_sq);
        const double arg = -kz * p.z_analysis + (d.kz - k0) * p.z0 - pi;
        e.s_matrix[k] = amplitudes[k + static_cast<std::size_t>(n_states * n_v) * i] * prefactor *
                        std::exp(Constants::imaginary_unit * arg) * std::sqrt(kz * k0 / weight);
        e.open[k] = 1;
      }
    combine_s_matrix(e, basis, gas);
    spectrum.energies.push_back(std::move(e));
  }
  return spectrum;
}

// Fortran: smatrix_calc (nread = 2, non-NAG), search for idat
GridInterval locate_on_grid(const FineAnalysisGrid& grid, double z_energy) {
  const int n = grid.size();
  int upper = 0;
  while (upper < n && grid.z_energy(upper) < z_energy) upper = upper + 1;
  upper = std::clamp(upper, 1, n - 1);
  const double e1 = grid.z_energy(upper - 1);
  return {upper, (z_energy - e1) / (grid.z_energy(upper) - e1)};
}

cplx interpolate_s_matrix(cplx s1, cplx s2, double f) {
  const double p1 = s1.real() * s1.real() + s1.imag() * s1.imag();
  const double p2 = s2.real() * s2.real() + s2.imag() * s2.imag();
  const double p = p1 + (p2 - p1) * f;  // as Fortran: prob_interpolated
  const cplx chord = s1 + (s2 - s1) * f;
  if (p <= 0.0 || chord == cplx(0.0, 0.0)) return {0.0, 0.0};
  return std::polar(std::sqrt(p), std::arg(chord));
}

int count_large_phase_steps(const ScatteringResults& spectrum) {
  double p_max = 0.0;
  for (const EnergyResult& e : spectrum.energies)
    for (const double p : e.probability) p_max = std::max(p_max, p);
  const double floor = 1e-8 * p_max;

  int count = 0;
  for (std::size_t i = 1; i < spectrum.energies.size(); ++i) {
    const EnergyResult& a = spectrum.energies[i - 1];
    const EnergyResult& b = spectrum.energies[i];
    for (std::size_t k = 0; k < a.s_matrix.size(); ++k) {
      if (a.probability[k] <= floor || b.probability[k] <= floor) continue;
      if (std::abs(std::arg(b.s_matrix[k] * std::conj(a.s_matrix[k]))) > Constants::pi / 2.0) count = count + 1;
    }
  }
  return count;
}

// Fortran: StateToStateProbs (nread = 2: range checks) and smatrix_calc (nread = 2, non-NAG:
// linear interpolation between the neighbouring grid points)
FineAnalysisResults interpolate_to_analysis_energies(const Config& config, const FineAnalysisGrid& grid,
                                                     const AnalysisBasis& basis, const GasPhase& gas,
                                                     const AnalysisEnergies& energies, ScatteringResults spectrum) {
  using fortran_compare::rge;
  using fortran_compare::rle;
  const double chaev = units::legacy.energy_factor();
  const DerivedParameters& d = config.derived;
  const int n = grid.size();

  FineAnalysisResults r;
  r.grid = grid;
  r.large_phase_steps = count_large_phase_steps(spectrum);
  for (int ien = 0; ien < energies.size(); ++ien) {
    const double collision = energies.collision[ien];
    const double ez = collision - basis.diffraction_energy(0) * chaev;  // etransZ
    if (n < 2) {
      r.skipped.push_back({ez, "fewer than 2 points on the fine-analysis grid"});
      continue;
    }
    if (!(rle(d.confidence_e_min, ez) && rge(d.confidence_e_max, ez))) {
      r.skipped.push_back({ez, "outside the 95% confidence interval"});
      continue;
    }
    if (!(rge(ez, grid.z_energy(0)) && rle(ez, grid.z_energy(n - 1)))) {
      r.skipped.push_back({ez, "outside the fine-analysis grid"});
      continue;
    }

    const auto [upper, f] = locate_on_grid(grid, ez);
    const EnergyResult& lo = spectrum.energies[upper - 1];
    const EnergyResult& hi = spectrum.energies[upper];

    EnergyResult e{};
    e.total_energy = energies.total[ien];
    e.collision_energy = collision;
    e.z_energy = ez;
    e.total_momentum = std::sqrt(2.0 * d.total_mass * collision / chaev);
    e.s_matrix.resize(lo.s_matrix.size());
    e.open.resize(lo.s_matrix.size());
    for (std::size_t k = 0; k < lo.s_matrix.size(); ++k) {
      e.s_matrix[k] = interpolate_s_matrix(lo.s_matrix[k], hi.s_matrix[k], f);
      e.open[k] = (lo.open[k] != 0 || hi.open[k] != 0) ? 1 : 0;
    }
    combine_s_matrix(e, basis, gas);
    r.interpolated.energies.push_back(std::move(e));
  }
  r.spectrum = std::move(spectrum);
  return r;
}

}  // namespace spodvr
