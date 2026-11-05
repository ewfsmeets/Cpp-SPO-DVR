#pragma once

#include <filesystem>
#include <ostream>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/fine_analysis.hpp"
#include "input/config.hpp"
#include "report/format.hpp"

namespace spodvr {

// The spectral grid of the fine analysis: records, transform length, resolution, energy range.
// Fortran: the stdout of AsympInfo.f90 ReadOrWriteAsympInfo (nread = 2)
void print_fine_analysis_grid(std::ostream& out, const Style& style, const Config& config,
                              const FineAnalysisGrid& grid);

// The spectrum (thinned to about 50 rows: E_z, summed and reaction probability, phase of the
// initial channel), the skipped analysis energies, warnings, and the result tables of the
// interpolated energies. Fortran: the stdout of probaf and StateToStateProbs (nread = 2)
void print_fine_analysis_results(std::ostream& out, const Style& style, const Config& config,
                                 const AnalysisBasis& basis, const FineAnalysisResults& results);

// Files in dir (created): spectrum/ (scattering tables on the grid; s_matrix.tsv only with
// writeSMatrix), the tables of the interpolated energies, and grid.json (grid and skipped
// energies). Returns the files written. Format: docs/output_files.md
std::vector<std::filesystem::path> write_fine_analysis_files(const std::filesystem::path& dir, const Config& config,
                                                             const AnalysisBasis& basis,
                                                             const FineAnalysisResults& results);

}  // namespace spodvr
