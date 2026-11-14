#pragma once

#include <filesystem>
#include <functional>
#include <optional>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/flux_analysis.hpp"
#include "input/config.hpp"
#include "io/binary_file.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Native checkpoint, BK coefficient and flux coefficient files (versioned binary formats, see
// docs/output_files.md). All start with the input parameters that a restart must not change
// (grids, masses, time step, wave packet, absorbers, ...); reading them with a different input
// is a CheckpointError naming the parameter.
// Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc (fort.9 / fort.8), AsympInfo.f90 and
// Analyzor.f90 (fort.21 / fort.22), FluxAnalyzor.f90 readWriteFluxCoeffs (flux_coeffs*.dat)

// The wave function at the end of a run plus the analysis surface state. Written to
// `path`.tmp first and then renamed, so an interrupted write leaves an existing file intact.
void write_checkpoint(const std::filesystem::path& path, const Config& config, double time, const WavePacket& psi,
                      const AnalysisSurfaceState& surface);

struct CheckpointInfo {
  double time;
  AnalysisSurfaceState surface;
};

// Reads a checkpoint into psi (whose shape must match the file).
[[nodiscard]] CheckpointInfo read_checkpoint(const std::filesystem::path& path, const Config& config, WavePacket& psi);

// The BK coefficient time series, one record (time, coefficients) per analysis step, appended
// and flushed as the propagation runs.
class BkFileWriter {
 public:
  BkFileWriter(const std::filesystem::path& path, const Config& config, const AnalysisBasis& basis);
  void append(double time, const BkCoefficients& bk);

 private:
  BinaryWriter file_;
};

// Calls on_record(time, coefficients) for every record of a BK file written for the same
// input and analysis basis. Returns the time of the last record (none: empty file).
std::optional<double> read_bk_file(const std::filesystem::path& path, const Config& config, const AnalysisBasis& basis,
                                   const std::function<void(double, const BkCoefficients&)>& on_record);

// The flux coefficients phi and phi' at `time`, the end of a run with the flux analysis: the
// header, the flux parameters (rFluxAnalysis, EMinFluxAna, EMaxFluxAna, DeltaEFluxAna), the
// time, the sizes and the two arrays. Written to `path`.tmp first and then renamed.
// Fortran: FluxAnalyzor.f90 readWriteFluxCoeffs ("write")
void write_flux_coefficients(const std::filesystem::path& path, const Config& config, double time,
                             const FluxCoefficients& coefficients);

// Reads the file in place into `coefficients` (sizes from the input) and returns its time. The
// flux parameters must equal the input's exactly, as the restart parameters of the header.
// Fortran: FluxAnalyzor.f90 readWriteFluxCoeffs ("read")
[[nodiscard]] double read_flux_coefficients(const std::filesystem::path& path, const Config& config,
                                            FluxCoefficients& coefficients);

}  // namespace spodvr
