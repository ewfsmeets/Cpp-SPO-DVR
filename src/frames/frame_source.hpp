#pragma once
// A frame of the wave packet (dtFrames, C++ only): one pass over the Z planes the backend hands out (z_plane for
// the scattering grid, the specular wave on the host), reduced by frames/frame_reducers.hpp, and timed.

#include "frames/frame_reducers.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/propagation_backend.hpp"
#include "propagation/propagation_tables.hpp"

namespace spodvr {

class FrameSource {
 public:
  // The tables must outlive the source (it keeps a view of their specular_to_scattering).
  FrameSource(const Config& config, const SpatialGrids& grids, const RotationalBasis& basis,
              const PropagationTables& tables);

  // The frame of the backend's wave packet now (valid until the next call).
  [[nodiscard]] const FrameRecord& frame(const PropagationBackend& backend);

  [[nodiscard]] const FrameLayout& layout() const { return layout_; }
  [[nodiscard]] int frames() const { return frames_; }           // frames made so far
  [[nodiscard]] double seconds() const { return seconds_; }      // their wall-clock time in all
  // The values a frame holds besides the wave packet (its record and the per-r sums): the memory estimate's item.
  [[nodiscard]] static std::size_t memory_values(const FrameLayout& layout);

 private:
  FrameLayout layout_;
  FrameRecord record_;
  std::vector<double> scratch_;
  int frames_ = 0;
  double seconds_ = 0.0;
};

// The layout of the frames of a run (sizes, cell volumes, s2s).
[[nodiscard]] FrameLayout frame_layout(const Config& config, const SpatialGrids& grids, const RotationalBasis& basis,
                                       std::span<const std::complex<double>> specular_to_scattering);

}  // namespace spodvr
