#pragma once

#include <filesystem>
#include <functional>

#include "analysis/analysis_basis.hpp"
#include "analysis/fine_analysis.hpp"
#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "input/config.hpp"
#include "propagation/memory_estimate.hpp"

namespace spodvr {

// The rotational and analysis bases of a run, built without the propagation (fine analysis).
// Not movable: the analysis basis refers to the rotational basis.
struct AnalysisSetup {
  AnalysisSetup(const Config& config, const AngularGrid& angular, const GasPhase& gas);
  AnalysisSetup(const AnalysisSetup&) = delete;
  AnalysisSetup& operator=(const AnalysisSetup&) = delete;

  RotationalBasis rotations;
  AnalysisBasis basis;
};

// A BK coefficient file: native (bk_coefficients.bin) or Fortran (fort.21 / fort.22).
struct BkInput {
  std::filesystem::path path;
  bool fortran_format = false;
};

// Memory of the large arrays of the fine analysis: the BK series, the transformed amplitudes,
// the spectrum (S, probabilities, phases at every grid point) and the FFT work array.
[[nodiscard]] MemoryEstimate estimate_fine_analysis_memory(const FineAnalysisGrid& grid, const AnalysisBasis& basis);

// The complete fine analysis of the BK file: grid (on_grid is called once it is known, before
// the large arrays are allocated), series, transform, spectrum and interpolation to the
// NrOfEs analysis energies. Fortran: Main.f90 (ReadOrWriteAsympInfo, gtaife, probaf,
// StateToStateProbs for nread = 2)
[[nodiscard]] FineAnalysisResults run_fine_analysis(const Config& config, const AnalysisSetup& setup,
                                                    const GasPhase& gas, const BkInput& input, const Fft& fft,
                                                    const std::function<void(const FineAnalysisGrid&)>& on_grid = {});

}  // namespace spodvr
