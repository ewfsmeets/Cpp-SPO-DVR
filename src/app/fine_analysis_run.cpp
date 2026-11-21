// The fine analysis of a BK coefficient file (run mode 2 and alsoAnalysis): reading, transform,
// spectrum and interpolation. Fortran: Main.f90 (nread = 2 part), AsympInfo.f90, Analyzor.f90
#include "app/fine_analysis_run.hpp"

#include <complex>

#include "grids/spatial_grids.hpp"
#include "io/checkpoint.hpp"
#include "io/fortran_files.hpp"
#include "propagation/propagation_tables.hpp"

namespace spodvr {

AnalysisSetup::AnalysisSetup(const Config& config, const AngularGrid& angular, const GasPhase& gas)
    : rotations(angular, config.input.j_in, config.input.mj_in),
      basis(config, angular, rotations, gas, parallel_kinetic_energies(config)) {}

MemoryEstimate estimate_fine_analysis_memory(const FineAnalysisGrid& grid, const AnalysisBasis& basis) {
  constexpr std::size_t complex_bytes = sizeof(std::complex<double>);
  const std::size_t n_series = std::size_t(basis.n_states()) * basis.n_v();
  const std::size_t n_points = grid.size();
  MemoryEstimate e;
  const std::size_t series = n_series * (grid.n_records_used + 1);
  e.items.push_back({"BK coefficient series", series, series * complex_bytes, "", true, true});
  const std::size_t amplitudes = n_series * n_points;
  e.items.push_back({"Transformed series", amplitudes, amplitudes * complex_bytes, "", true, true});
  // Per grid point and state: S (complex), probability, phase (double), open flag (byte).
  const std::size_t spectrum_bytes = amplitudes * (complex_bytes + 2 * sizeof(double) + 1);
  e.items.push_back({"Spectrum (S, P, phase)", amplitudes, spectrum_bytes, "", true, true});
  const std::size_t work = std::size_t(grid.transform_length) * 64;
  e.items.push_back({"FFT work array", work, work * complex_bytes, "64 series per FFT", true, true});
  for (const MemoryItem& item : e.items) e.setup_bytes = e.setup_bytes + item.bytes;
  e.peak_bytes = e.setup_bytes;
  return e;
}

namespace {

// Calls on_record for every record of the file.
int read_records(const Config& config, const AnalysisSetup& setup, const BkInput& input,
                 const std::function<void(double, const BkCoefficients&)>& on_record) {
  if (input.fortran_format) return read_fortran_bk_file(input.path, config, setup.rotations, setup.basis, on_record);
  int n = 0;
  (void)read_bk_file(input.path, config, setup.basis, [&](double time, const BkCoefficients& bk) {
    n = n + 1;
    on_record(time, bk);
  });
  return n;
}

// The fixed specular analysis line of the run that wrote the BK file, from the analysis surface in its
// checkpoint; nothing without a checkpoint or when the line never became fixed.
std::optional<double> specular_line(const Config& config, const BkInput& input) {
  if (!input.checkpoint) return std::nullopt;
  const CheckpointInput& checkpoint = *input.checkpoint;
  const AnalysisSurfaceState surface =
      checkpoint.fortran_format
          ? read_fortran_analysis_surface(checkpoint.path, config.input.nr, config.input.nz_specular)
          : read_checkpoint_surface(checkpoint.path, config);
  if (!surface.fixed) return std::nullopt;
  return make_spatial_grids(config).z_specular[surface.index];
}

}  // namespace

FineAnalysisResults run_fine_analysis(const Config& config, const AnalysisSetup& setup, const GasPhase& gas,
                                      const BkInput& input, const Fft& fft,
                                      const std::function<void(const FineAnalysisGrid&)>& on_grid) {
  const AnalysisBasis& basis = setup.basis;
  // First pass: the number of records (and all header checks). Fortran: ReadOrWriteAsympInfo
  const int n_records = read_records(config, setup, input, [](double, const BkCoefficients&) {});
  const FineAnalysisGrid grid = make_fine_analysis_grid(config, n_records);
  // After the BK file's checks, so that a BK file of another input is reported as such.
  const std::optional<double> line = specular_line(config, input);
  if (on_grid) on_grid(grid);

  // Second pass: the series. Fortran: getbkc (ev0jl = enymo(inchna, nvin))
  const double initial_energy = basis.channel_energy(basis.initial_state(), config.derived.v_in_index - 1);
  BkSeries series(grid, basis, initial_energy);
  (void)read_records(config, setup, input, [&](double time, const BkCoefficients& bk) { series.add_record(time, bk); });
  series.check_complete();

  const std::vector<std::complex<double>> amplitudes = half_fourier_transform(series, grid, fft);
  ScatteringResults spectrum = fine_spectrum(config, grid.z_energies(), basis, gas, amplitudes, line);
  return interpolate_to_analysis_energies(config, grid, basis, gas, make_analysis_energies(config, gas),
                                          std::move(spectrum));
}

}  // namespace spodvr
