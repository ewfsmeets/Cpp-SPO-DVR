// A frame of the wave packet (dtFrames): one pass over the backend's Z planes, reduced on the host. C++ only.
#include "frames/frame_source.hpp"

#include <algorithm>
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

FrameSource::FrameSource(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                         const RotationalBasis& basis, const GasPhase& gas, const PropagationTables& tables,
                         const Fft& fft)
    : layout_(frame_layout(config, grids, basis, tables.specular_to_scattering)),
      record_(make_frame_record(layout_)),
      scratch_(std::size_t(layout_.nr) * layout_.n_rot) {
  if (!FrameExtraSwitches::of(config.input).any()) return;
  extras_.emplace(config, grids, angular, basis, gas, tables, fft, layout_);
  const FrameExtraSizes sizes = frame_extra_sizes(config, grids, angular, basis);
  record_.prob_angle_z.assign(sizes.angle_z, 0.0);
  record_.pop_nm_z.assign(sizes.nm_z, 0.0);
  record_.pop_vj_z.assign(sizes.vj_z, 0.0);
  record_.p2_numerator_z.assign(sizes.moments_z, 0.0);
  record_.rot_numerator_z.assign(sizes.moments_z, 0.0);
  record_.prob_kr.assign(sizes.kr, 0.0);
  record_.prob_kz_spec.assign(sizes.kz_spec, 0.0);
  full_plane_.resize(layout_.plane().size());
}

// full_plane_ = the plane of psi_scat at iz with s2s chi(ir, iz) added in the initial channel (psi_full).
void FrameSource::add_specular_part(int iz, std::span<const std::complex<double>> plane, const SpecularWave& specular) {
  const Index3D in_plane = layout_.plane();
  std::ranges::copy(plane, full_plane_.begin());
  for (int ir = 0; ir < layout_.nr; ++ir)
    for (int ixy = 0; ixy < layout_.nxy; ++ixy)
      full_plane_[in_plane(layout_.initial_channel, ixy, ir)] += layout_.specular_to_scattering[ixy] * specular(ir, iz);
}

const FrameRecord& FrameSource::frame(const PropagationBackend& backend) {
  const auto start = std::chrono::steady_clock::now();
  const SpecularWave& specular = backend.specular();
  for (int iz = 0; iz < layout_.nz; ++iz) {
    const std::span<const std::complex<double>> plane = backend.z_plane(iz);
    reduce_scattering_plane(layout_, iz, plane, specular, scratch_, record_);
    if (!extras_) continue;
    add_specular_part(iz, plane, specular);
    extras_->reduce_plane(iz, full_plane_, record_);
  }
  reduce_specular_tail(layout_, specular, record_);
  if (extras_) extras_->reduce_specular(specular, record_);
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  seconds_ += elapsed.count();
  ++frames_;
  return record_;
}

std::size_t FrameSource::memory_values(const FrameLayout& layout) {
  return layout.zrxy().size() + layout.tail_zr().size() + layout.jmj_z().size() + std::size_t(layout.nr) * layout.n_rot;
}

}  // namespace spodvr
