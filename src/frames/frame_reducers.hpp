#pragma once
// The reductions of one frame (dtFrames, C++ only): the wave function psi_full = psi_scat + s2s chi (chi added in the
// initial rotational channel) reduced on the host, one Z plane at a time, as the backend hands the planes out.

#include <complex>
#include <span>
#include <vector>

#include "core/array_index.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// The sizes and factors of the reductions.
struct FrameLayout {
  int nz;               // points of the scattering Z grid
  int nz_specular;      // points of the specular Z grid (>= nz; both start at ZStart with the same dZ)
  int nr;
  int nxy;              // nX nY
  int n_rot;            // rotational states (J, mJ)
  int initial_channel;  // the state (JIn, MJIn), which holds chi
  double cell_volume;           // dZ dr dx dy sin(skew): |psi|^2 at a (Z, r, x, y) point times it is a probability
  double specular_cell_volume;  // dZ dr: |chi|^2 at a (Z, r) point times it is a probability
  // s2s(ixy) = exp(i k_par . R(ixy)) unit_cell_norm: chi on the scattering grid. Fortran: TwoGrids.f90
  std::span<const std::complex<double>> specular_to_scattering;

  // The layout of one Z plane of the scattering wave as the backend hands it out: (irot, ixy, ir).
  [[nodiscard]] Index3D plane() const { return {std::size_t(n_rot), std::size_t(nxy), std::size_t(nr)}; }
  // prob_zrxy [nZ, nR, nX, nY]: ixy fastest, then ir, then iz.
  [[nodiscard]] Index3D zrxy() const { return {std::size_t(nxy), std::size_t(nr), std::size_t(nz)}; }
  // pop_jmj_z [nZ, nrot]: irot fastest.
  [[nodiscard]] Index2D jmj_z() const { return {std::size_t(n_rot), std::size_t(nz)}; }
  // prob_tail_zr [nZSp - nZ, nR]: ir fastest; the first row is Z index nz.
  [[nodiscard]] Index2D tail_zr() const { return {std::size_t(nr), std::size_t(nz_specular - nz)}; }
};

// One frame's datasets (docs/output_files.md, "frames/").
struct FrameRecord {
  std::vector<double> prob_zrxy;     // probability per (Z, r, x, y) cell of psi_full
  std::vector<double> prob_tail_zr;  // probability per (Z, r) cell of chi beyond the scattering grid
  std::vector<double> pop_jmj_z;     // probability per (Z, rotational state) of psi_full
};

// A record of the layout's sizes (zeros).
[[nodiscard]] FrameRecord make_frame_record(const FrameLayout& layout);

// The prob_zrxy and pop_jmj_z entries of Z index iz < nZ from that plane of the scattering wave and chi.
// `scratch` holds nr n_rot values (the per-r sums of pop_jmj_z). Every output element is written by one thread and
// summed in a fixed order: the result does not depend on the number of threads.
void reduce_scattering_plane(const FrameLayout& layout, int iz, std::span<const std::complex<double>> plane,
                             const SpecularWave& specular, std::span<double> scratch, FrameRecord& record);

// prob_tail_zr: |chi|^2 dZ dr at the specular grid's points beyond the scattering grid (nothing counted twice).
void reduce_specular_tail(const FrameLayout& layout, const SpecularWave& specular, FrameRecord& record);

}  // namespace spodvr
