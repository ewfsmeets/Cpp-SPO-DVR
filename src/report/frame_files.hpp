#pragma once
// frames/ (dtFrames, C++ only): the frame records, a series directory (io/series_writer.hpp), their projected size
// and cost, and the lines that report their measured cost (docs/output_files.md, "frames/").

#include <filesystem>
#include <optional>
#include <ostream>

#include "frames/frame_reducers.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "report/format.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

class FrameFiles {
 public:
  FrameFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
             const RotationalBasis& basis, const FrameLayout& layout);
  void append(double time, const FrameRecord& record, int surface_index);

 private:
  SeriesWriter series_;
};

// The size of frames/ (records times record bytes) and its projected cost, or nothing without dtFrames. The cost
// scales the per-value times measured on the thesis DZ grid with 8 threads (0.52 ns per complex value of psi for a
// frame on the CPU, 2.45 ns on the GPU with its plane copies; a time step 10.5 and 1.84 ns).
[[nodiscard]] std::optional<PlotDataSize> frames_data_size(const Config& config, const SpatialGrids& grids,
                                                           const RotationalBasis& basis);

// The measured cost of the frames at the end of the run: the frames written and the wall-clock time of all frames
// computed.
void print_frame_cost(std::ostream& out, int frames_written, double seconds, double propagation_seconds);

}  // namespace spodvr
