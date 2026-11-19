#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "propagation/averages.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

// densities/ (writeDensities): the probability per (Z, r) cell of the scattering and the specular grid,
// the analysis-surface index and the norm, at every averages time (docs/output_files.md, "densities/").
class DensitiesFiles {
 public:
  DensitiesFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids);
  void append(const Averages& averages, const CellProbabilities& cells, int surface_index);

 private:
  SeriesWriter series_;
};

// The size of densities/ (records times record bytes), or nothing without writeDensities.
[[nodiscard]] std::optional<PlotDataSize> densities_data_size(const Config& config, const SpatialGrids& grids);

}  // namespace spodvr
