#pragma once

// The GPU (CUDA) propagation backend: the entry points the rest of the program uses. This header
// needs no CUDA headers.

#include <cstddef>
#include <memory>
#include <string>

#include "propagation/memory_estimate.hpp"
#include "propagation/propagation_backend.hpp"

namespace spodvr::cuda {

// Selects the first GPU and sets CUDA_MODULE_LOADING=LAZY (unless set: loads only the kernels
// that are used, which saves GPU memory). Throws CudaError if there is no usable GPU.
void initialise_cuda();

// One line describing the GPU and the CUDA libraries.
[[nodiscard]] std::string cuda_device_summary();

// Free and total GPU memory (bytes).
[[nodiscard]] std::size_t free_gpu_memory();
[[nodiscard]] std::size_t total_gpu_memory();

// The GPU memory the backend allocates for a run, array by array (the cuFFT work area at its
// limit; the CUDA context itself, a few hundred MiB, is not included).
[[nodiscard]] MemoryEstimate estimate_gpu_memory(const Config& config, const SpatialGrids& grids,
                                                 const AngularGrid& angular, int n_representatives);

// The CUDA backend as a BackendFactory: the wave packet lives in GPU memory for the whole run.
[[nodiscard]] std::unique_ptr<PropagationBackend> make_cuda_backend(const BackendInput& input, PotentialGrid potential,
                                                                    SpecularWave initial);

}  // namespace spodvr::cuda
