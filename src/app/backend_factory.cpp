// Choice of the propagation backend (--device). The only place that knows whether the CUDA
// backend was built.
#include "app/backend_factory.hpp"

#include "core/error.hpp"
#include "propagation/cpu_backend.hpp"
#include "report/sections.hpp"
#ifdef SPODVR_HAVE_CUDA
#include "cuda/cuda_backend.hpp"
#endif

namespace spodvr {

bool cuda_backend_available() {
#ifdef SPODVR_HAVE_CUDA
  return true;
#else
  return false;
#endif
}

BackendFactory backend_factory(Device device) {
  if (device == Device::Cpu) return make_cpu_backend;
#ifdef SPODVR_HAVE_CUDA
  return cuda::make_cuda_backend;
#else
  throw NotSupportedError("--device cuda: this program was built without CUDA support (cmake -DSPODVR_USE_CUDA=ON)");
#endif
}

std::string prepare_device(Device device) {
  if (device == Device::Cpu) return "";
#ifdef SPODVR_HAVE_CUDA
  cuda::initialise_cuda();
  return cuda::cuda_device_summary();
#else
  (void)backend_factory(device);  // throws
  return "";
#endif
}

void check_device_memory([[maybe_unused]] std::ostream& out, [[maybe_unused]] const Style& style,
                         [[maybe_unused]] Device device, [[maybe_unused]] const Config& config,
                         [[maybe_unused]] const SpatialGrids& grids, [[maybe_unused]] const AngularGrid& angular,
                         [[maybe_unused]] int n_representatives) {
#ifdef SPODVR_HAVE_CUDA
  if (device != Device::Cuda) return;
  const MemoryEstimate estimate = cuda::estimate_gpu_memory(config, grids, angular, n_representatives);
  print_memory_estimate(out, style, estimate, cuda::total_gpu_memory(), "GPU MEMORY ESTIMATE");
  const std::size_t free = cuda::free_gpu_memory();
  if (estimate.peak_bytes > free)
    throw NotSupportedError("--device cuda: this run needs about " + format_bytes(estimate.peak_bytes) +
                            " of GPU memory, but only " + format_bytes(free) +
                            " are free (use --device cpu, surface symmetry or a larger GPU)");
#endif
}

}  // namespace spodvr
