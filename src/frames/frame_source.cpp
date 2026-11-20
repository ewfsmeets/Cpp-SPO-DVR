// A frame of the wave packet (dtFrames): one pass over the backend's Z planes, reduced on the host. C++ only.
#include "frames/frame_source.hpp"

#include <chrono>

namespace spodvr {

FrameLayout frame_layout(const Config& config, const SpatialGrids& grids, const RotationalBasis& basis,
                         std::span<const std::complex<double>> specular_to_scattering) {
  const double dz_dr = grids.z_scattering.spacing() * grids.r.spacing();
  return {grids.z_scattering.size(),
          grids.z_specular.size(),
          grids.r.size(),
          config.input.nx * config.input.ny,
          basis.size(),
          basis.initial_index(),
          dz_dr * config.derived.xy_volume_element,
          dz_dr,
          specular_to_scattering};
}

FrameSource::FrameSource(const Config& config, const SpatialGrids& grids, const RotationalBasis& basis,
                         const PropagationTables& tables)
    : layout_(frame_layout(config, grids, basis, tables.specular_to_scattering)),
      record_(make_frame_record(layout_)),
      scratch_(std::size_t(layout_.nr) * layout_.n_rot) {}

const FrameRecord& FrameSource::frame(const PropagationBackend& backend) {
  const auto start = std::chrono::steady_clock::now();
  const SpecularWave& specular = backend.specular();
  for (int iz = 0; iz < layout_.nz; ++iz)
    reduce_scattering_plane(layout_, iz, backend.z_plane(iz), specular, scratch_, record_);
  reduce_specular_tail(layout_, specular, record_);
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  seconds_ += elapsed.count();
  ++frames_;
  return record_;
}

std::size_t FrameSource::memory_values(const FrameLayout& layout) {
  return layout.zrxy().size() + layout.tail_zr().size() + layout.jmj_z().size() + std::size_t(layout.nr) * layout.n_rot;
}

}  // namespace spodvr
