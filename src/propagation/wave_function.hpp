#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "core/aligned_vector.hpp"
#include "core/array_index.hpp"

namespace spodvr {

// Wave function on the scattering grid, Fortran layout (irot, ixy, izr) with ixy = iy + ny*ix
// and izr = ir + nr*iz (first index fastest). One "block" holds all (irot, ixy) values of one
// (Z, r) point. Fortran: WaveFunction.f90 ScatWave
struct ScatteringWave {
  Index3D index;
  ComplexVector values;  // 64-byte aligned (fast SIMD FFTs)

  // Zero wave function. Initialised in parallel per group of `group_points` consecutive (Z, r)
  // points, with the schedule of the propagation loops over those groups (first touch, see
  // AlignedAllocator and core/loop_schedule.hpp; the propagator's groups: xy_group_points).
  ScatteringWave(int n_rot, int n_xy, int n_zr, int group_points = 1) : index(n_rot, n_xy, n_zr), values(index.size()) {
    const std::size_t group = index.n0 * index.n1 * group_points;
#pragma omp parallel for schedule(runtime)
    for (int g = 0; g < n_zr / group_points; ++g)
      std::fill_n(values.begin() + static_cast<std::ptrdiff_t>(group * g), group, std::complex<double>(0.0, 0.0));
  }

  [[nodiscard]] int block_size() const { return static_cast<int>(index.n0 * index.n1); }
  [[nodiscard]] int n_blocks() const { return static_cast<int>(index.n2); }
  [[nodiscard]] std::span<std::complex<double>> block(int izr) {
    return std::span(values).subspan(izr * index.n0 * index.n1, index.n0 * index.n1);
  }
  [[nodiscard]] std::span<const std::complex<double>> block(int izr) const {
    return std::span(values).subspan(izr * index.n0 * index.n1, index.n0 * index.n1);
  }
};

// Wave function of the incoming channel on the specular (Z, r) grid, layout (ir, iz).
// Fortran: WaveFunction.f90 SpecWave
struct SpecularWave {
  Index2D index;
  ComplexVector values;  // 64-byte aligned (fast SIMD FFTs)

  SpecularWave(int n_r, int n_z) : index(n_r, n_z), values(index.size(), std::complex<double>(0.0, 0.0)) {}

  [[nodiscard]] int n_r() const { return static_cast<int>(index.n0); }
  [[nodiscard]] int n_z() const { return static_cast<int>(index.n1); }
  [[nodiscard]] std::complex<double>& operator()(int ir, int iz) { return values[index(ir, iz)]; }
  [[nodiscard]] const std::complex<double>& operator()(int ir, int iz) const { return values[index(ir, iz)]; }
};

// The complete wave packet: the scattering grid part plus the specular (incoming channel) part.
struct WavePacket {
  ScatteringWave scattering;
  SpecularWave specular;
};

}  // namespace spodvr
