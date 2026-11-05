#pragma once

#include "cuda/device_buffer.hpp"
#include "cuda/device_complex.hpp"

namespace spodvr::cuda {

// The wave packet in GPU memory: the scattering wave in the device layout (izr, irot, ixy) (see
// DeviceLayout) and the specular wave in its host layout (ir, iz).
struct DeviceWavePacket {
  DeviceLayout layout;
  int nr;
  int nz_specular;
  DeviceBuffer<Complex> scattering;
  DeviceBuffer<Complex> specular;

  DeviceWavePacket(const DeviceLayout& shape, int n_r, int n_z_specular)
      : layout(shape),
        nr(n_r),
        nz_specular(n_z_specular),
        scattering(shape.size()),
        specular(std::size_t(n_r) * n_z_specular) {}
};

}  // namespace spodvr::cuda
