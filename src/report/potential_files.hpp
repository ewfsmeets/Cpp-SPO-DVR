#pragma once

#include <filesystem>
#include <optional>

#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

// potential/ (writePotential): the run's potential on its grid as the propagation uses it (switching and the
// VMin/VMax limits included), written once at set-up (docs/output_files.md, "potential/"): the minimum over the
// DVR orientations and its orientation at every grid point of the interaction region, and V at every DVR
// orientation at the PlotSites. Written point row by point row: no large buffer. C++ only.
void write_potential_files(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
                           const AngularGrid& angular, const PotentialGrid& potential);

// The size of potential/, or nothing without writePotential.
[[nodiscard]] std::optional<PlotDataSize> potential_data_size(const Config& config, const SpatialGrids& grids,
                                                              const AngularGrid& angular);

}  // namespace spodvr
