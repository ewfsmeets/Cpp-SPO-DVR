#pragma once
// frames/ (dtFrames, C++ only): the frame records, a series directory (io/series_writer.hpp), their projected size
// and cost, and the lines that report their measured cost (docs/output_files.md, "frames/").

#include <filesystem>
#include <optional>
#include <ostream>

#include "frames/frame_extras.hpp"
#include "frames/frame_reducers.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "report/format.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

// Writes frames/: the manifest at construction, then one record per frame.
class FrameFiles {
 public:
  FrameFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
             const AngularGrid& angular, const RotationalBasis& basis, const FrameLayout& layout);
  // The frame at `time`: its record and the Z index of the specular analysis surface then.
  void append(double time, const FrameRecord& record, int surface_index);

 private:
  SeriesWriter series_;
};

// The size of frames/ (records times record bytes) and its projected cost, or nothing without dtFrames. The cost
// scales the per-value times measured on the thesis DZ grid with 8 threads (0.52 ns per complex value of psi for a
// frame on the CPU, 2.45 ns on the GPU with its plane copies; a time step 10.5 and 1.84 ns), plus the FrameExtras
// that are on (host work, measured on the same grid).
[[nodiscard]] std::optional<PlotDataSize> frames_data_size(const Config& config, const SpatialGrids& grids,
                                                           const AngularGrid& angular, const RotationalBasis& basis);

// The measured cost of the frames at the end of the run: the frames written and the wall-clock time of all frames
// computed.
void print_frame_cost(std::ostream& out, int frames_written, double seconds, double propagation_seconds);

// FrameExtras: the measured cost of each extra that is on, per frame (part of the frames' cost above).
void print_frame_extras_cost(std::ostream& out, int frames_written, const std::vector<double>& seconds);

}  // namespace spodvr
