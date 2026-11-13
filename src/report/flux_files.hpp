#pragma once

#include <filesystem>
#include <vector>

#include "analysis/flux_analysis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"

namespace spodvr {

// Writes the results of the flux analysis as tab-separated tables (a '#' comment block with units and
// column meanings, a header line, one row per record) plus a JSON metadata file:
//   flux_reaction.tsv    one row per flux energy: E_z, P_flux
//   flux_density.tsv     one row per (energy, ix, iy), energy outer and iy fastest as in the Fortran:
//                        E_z, ix, iy, x, y (skewed), X, Y (Cartesian), density
//   flux_analysis.json   run metadata, flux surface, energies, surface grid, file list
// Returns the files written. Format description: docs/output_files.md. Throws std::invalid_argument
// without a flux surface (doFluxAnalysis off) or for results of another surface grid.
// Fortran: FluxAnalyzor.f90 computeFluxProbability (flux_reaction.dat, flux_density.dat)
std::vector<std::filesystem::path> write_flux_files(const std::filesystem::path& dir, const Config& config,
                                                    const SpatialGrids& grids, const FluxResults& results);

}  // namespace spodvr
