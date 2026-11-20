// The reductions of one frame (dtFrames): probabilities per (Z, r, x, y) cell and per (Z, rotational state) of
// psi_full = psi_scat + s2s chi, and chi beyond the scattering grid. C++ only (plot data).
#include "frames/frame_reducers.hpp"

namespace spodvr {

namespace {

// |c|^2 as ABS(c)**2, as the averages and densities/ compute it (averages.cpp), so that prob_tail_zr equals the
// densities' specular cells bitwise. Only for the tail: a hypot per value would triple the cost of the scattering
// planes (DZ: 0.61 instead of 0.2 s per frame), whose checks are sums to 1e-13, not bitwise.
double abs_squared(const std::complex<double>& c) {
  const double a = std::abs(c);
  return a * a;
}

}  // namespace

FrameRecord make_frame_record(const FrameLayout& layout) {
  return {std::vector<double>(layout.zrxy().size()), std::vector<double>(layout.tail_zr().size()),
          std::vector<double>(layout.jmj_z().size())};
}

void reduce_scattering_plane(const FrameLayout& layout, int iz, std::span<const std::complex<double>> plane,
                             const SpecularWave& specular, std::span<double> scratch, FrameRecord& record) {
  const Index3D in_plane = layout.plane();
  const Index3D zrxy = layout.zrxy();
  const Index2D jmj_z = layout.jmj_z();
  const Index2D sums = {std::size_t(layout.n_rot), std::size_t(layout.nr)};  // scratch: (irot, ir)
#pragma omp parallel for schedule(static)
  for (int ir = 0; ir < layout.nr; ++ir) {
    const std::complex<double> chi = specular(ir, iz);
    for (int irot = 0; irot < layout.n_rot; ++irot) scratch[sums(irot, ir)] = 0.0;
    for (int ixy = 0; ixy < layout.nxy; ++ixy) {
      double cell = 0.0;
      for (int irot = 0; irot < layout.n_rot; ++irot) {
        std::complex<double> amplitude = plane[in_plane(irot, ixy, ir)];
        if (irot == layout.initial_channel) amplitude += layout.specular_to_scattering[ixy] * chi;  // psi_full
        const double density = std::norm(amplitude);
        cell += density;
        scratch[sums(irot, ir)] += density;
      }
      record.prob_zrxy[zrxy(ixy, ir, iz)] = cell * layout.cell_volume;
    }
  }
  for (int irot = 0; irot < layout.n_rot; ++irot) {
    double sum = 0.0;
    for (int ir = 0; ir < layout.nr; ++ir) sum += scratch[sums(irot, ir)];  // in r order, whatever the threads
    record.pop_jmj_z[jmj_z(irot, iz)] = sum * layout.cell_volume;
  }
}

void reduce_specular_tail(const FrameLayout& layout, const SpecularWave& specular, FrameRecord& record) {
  const Index2D tail = layout.tail_zr();
  for (int iz = layout.nz; iz < layout.nz_specular; ++iz)
    for (int ir = 0; ir < layout.nr; ++ir)
      record.prob_tail_zr[tail(ir, iz - layout.nz)] = abs_squared(specular(ir, iz)) * layout.specular_cell_volume;
}

}  // namespace spodvr
