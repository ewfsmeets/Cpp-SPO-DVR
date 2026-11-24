#pragma once

#include <cstdint>
#include <filesystem>
#include <ostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "report/format.hpp"

namespace spodvr {

// The sites of PlotSites as the plot-data manifests record them: name, fractions x, y, grid point ix, iy and
// ixy = iy + nY ix.
[[nodiscard]] nlohmann::ordered_json plot_sites_json(const InputParameters& p);

// Measured on the thesis DZ grid (8 threads; the GPU an RTX 4070 Ti), the plot plan's measurements: seconds per
// complex value of psi (nZ nR nX nY nrot) for one time step, the scale of the plot data's projected costs, and the
// share of the propagation above which a PLOT DATA line warns.
inline constexpr double step_seconds_per_value_cpu = 10.5e-9;
inline constexpr double step_seconds_per_value_gpu = 1.84e-9;
inline constexpr double cost_warning_fraction = 0.02;

// The plot data a run will write: the directory, the key that asks for it, the projected bytes and what it
// costs (docs/output_files.md, "Series directories").
struct PlotDataSize {
  std::string directory;
  std::string key;
  std::uintmax_t bytes = 0;
  std::string cost;  // when it is written and what it adds to the run
};

// The plot data of the run's switches (empty when none is on): densities/ (writeDensities), potential/
// (writePotential), frames/ (dtFrames), norm_budget/ (NormBudget) and energy_resolved/ (EnergyResolvedSites).
[[nodiscard]] std::vector<PlotDataSize> plot_data_sizes(const Config& config, const SpatialGrids& grids,
                                                        const AngularGrid& angular);

// The free space of the file system of `output_dir` (of its nearest existing ancestor when it does not exist yet).
[[nodiscard]] std::uintmax_t free_space(const std::filesystem::path& output_dir);

// The bytes of the plot-data directories in `output_dir` that the run will replace (removed when it starts
// writing them, so their space counts as free).
[[nodiscard]] std::uintmax_t replaced_plot_data_bytes(const std::filesystem::path& output_dir,
                                                      const std::vector<PlotDataSize>& sizes);

// Throws InputError naming its key when a plot-data directory exists in `output_dir` but is not plot data this
// program wrote (no manifest.json of that kind): it would be removed when the run writes its own.
void check_replaceable(const std::filesystem::path& output_dir, const std::vector<PlotDataSize>& sizes);

// Throws InputError naming the key that adds the most bytes when the plot data exceed `available` bytes.
void check_free_space(const std::vector<PlotDataSize>& sizes, std::uintmax_t available);

// The set-up lines of the plot data: per directory its size and its cost.
void print_plot_data(std::ostream& out, const Style& style, const std::vector<PlotDataSize>& sizes);

}  // namespace spodvr
