#pragma once

#include "propagation/propagation_backend.hpp"

namespace spodvr {

// The backend factory for a device. Device::Cuda throws NotSupportedError when the program was
// built without CUDA support (cmake -DSPODVR_USE_CUDA=ON).
[[nodiscard]] BackendFactory backend_factory(Device device);

// Whether this build contains the CUDA backend.
[[nodiscard]] bool cuda_backend_available();

}  // namespace spodvr
