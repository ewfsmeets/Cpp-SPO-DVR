#pragma once

#include <cstdint>
#include <filesystem>
#include <ostream>
#include <string>
#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "propagation/averages.hpp"
#include "report/format.hpp"

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

// The plot data a run will write (directory, the key that asks for it, projected bytes).
struct PlotDataSize {
  std::string directory;
  std::string key;
  std::uintmax_t bytes = 0;
};

// The plot data of the run's switches (empty when none is on): records times record bytes.
[[nodiscard]] std::vector<PlotDataSize> plot_data_sizes(const Config& config, const SpatialGrids& grids);

// The free space of the file system of `output_dir` (of its nearest existing ancestor when it does not exist yet).
[[nodiscard]] std::uintmax_t free_space(const std::filesystem::path& output_dir);

// The bytes of the plot-data directories in `output_dir` that the run will replace (removed when it starts
// writing them, so their space counts as free).
[[nodiscard]] std::uintmax_t replaced_plot_data_bytes(const std::filesystem::path& output_dir,
                                                      const std::vector<PlotDataSize>& sizes);

// Throws InputError naming the key that adds the most bytes when the plot data exceed `available` bytes.
void check_free_space(const std::vector<PlotDataSize>& sizes, std::uintmax_t available);

// The set-up lines of the plot data: per series its size and the cost in the time loop.
void print_plot_data(std::ostream& out, const Style& style, const std::vector<PlotDataSize>& sizes);

}  // namespace spodvr
