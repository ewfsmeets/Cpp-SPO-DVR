#pragma once

#include <ostream>
#include <string>

#include "propagation/propagation_backend.hpp"
#include "report/format.hpp"

namespace spodvr {

// The backend factory for a device. Device::Cuda throws NotSupportedError when the program was
// built without CUDA support (cmake -DSPODVR_USE_CUDA=ON).
[[nodiscard]] BackendFactory backend_factory(Device device);

// Whether this build contains the CUDA backend.
[[nodiscard]] bool cuda_backend_available();

// Prepares the device (GPU: selects it) and returns a line describing it (empty for the CPU).
[[nodiscard]] std::string prepare_device(Device device);

// GPU: prints the GPU memory estimate and stops with a clear error if the free GPU memory is too
// small for this run (no-op for the CPU). Called before the potential is computed, so that a run
// does not fail after a long set-up.
void check_device_memory(std::ostream& out, const Style& style, Device device, const Config& config,
                         const SpatialGrids& grids, const AngularGrid& angular, int n_representatives);

}  // namespace spodvr
