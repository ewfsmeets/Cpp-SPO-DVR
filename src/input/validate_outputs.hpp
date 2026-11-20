#pragma once

#include <string>
#include <vector>

#include "input/config.hpp"

namespace spodvr {

// The checks of the plot-data keys (writeDensities, writePotential, PlotSites; docs/input_schema.md), which need
// the run mode and the derived parameters: called after miscellaneous_setup. C++ only (no Fortran counterpart);
// kept out of validate_input, whose checks mirror InputNew.f90 one by one.
void validate_output_keys(const Config& config);

// The names FrameExtras accepts: the optional datasets of frames/ (docs/output_files.md, "frames/").
[[nodiscard]] const std::vector<std::string>& frame_extra_names();

// The fields of a PlotSites entry ("name", "x", "y"), read by json_input.cpp in this order.
[[nodiscard]] const std::vector<std::string>& plot_site_fields();
[[nodiscard]] bool is_plot_site_field(const std::string& field);

// The most sites PlotSites may name.
constexpr int max_plot_sites = 16;

// The (x, y) grid point of a site that passed validate_output_keys: ixy = iy + ny * ix (Fortran order).
[[nodiscard]] int plot_site_ixy(const PlotSite& site, int nx, int ny);

}  // namespace spodvr
