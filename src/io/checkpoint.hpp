#pragma once

#include <filesystem>
#include <functional>
#include <optional>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "input/config.hpp"
#include "io/binary_file.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Native checkpoint and BK coefficient files (versioned binary formats, see
// docs/output_files.md). Both start with the input parameters that a restart must not change
// (grids, masses, time step, wave packet, absorbers, ...); reading them with a different input
// is a CheckpointError naming the parameter.
// Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc (fort.9 / fort.8), AsympInfo.f90 and
// Analyzor.f90 (fort.21 / fort.22)

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

}  // namespace spodvr
