#pragma once

#include <filesystem>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/s_matrix.hpp"
#include "input/config.hpp"

namespace spodvr {

// Writes the scattering results as tab-separated tables (one '#' comment block with units and
// column meanings, then a header line and one row per record) plus a JSON metadata file:
//   scattering_summary.tsv          one row per energy: total probabilities and averages
//   s_matrix.tsv                    one row per (energy, v, j, mJ, n, m) open channel
//   probabilities_rovib.tsv         per (energy, v, j, mJ), summed over diffraction (writeProbs)
//   probabilities_diffraction.tsv   per (energy, v, j, mJ, order) (writeRDProbs)
//   analysis.json                   run and basis metadata, file and column descriptions
// Returns the files written. Format description: docs/output_files.md
std::vector<std::filesystem::path> write_scattering_files(const std::filesystem::path& dir, const Config& config,
                                                          const AnalysisBasis& basis,
                                                          const ScatteringResults& results);

}  // namespace spodvr
