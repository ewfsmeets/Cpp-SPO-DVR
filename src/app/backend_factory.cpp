// Choice of the propagation backend (--device). The only place that knows whether the CUDA
// backend was built.
#include "app/backend_factory.hpp"

#include "core/error.hpp"
#include "propagation/cpu_backend.hpp"
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

void check_device_memory([[maybe_unused]] Device device, [[maybe_unused]] const Config& config,
                         [[maybe_unused]] const SpatialGrids& grids, [[maybe_unused]] const AngularGrid& angular,
                         [[maybe_unused]] int n_representatives) {
#ifdef SPODVR_HAVE_CUDA
  if (device != Device::Cuda) return;
  const std::size_t needed = cuda::estimate_gpu_memory(config, grids, angular, n_representatives);
  const std::size_t free = cuda::free_gpu_memory();
  if (needed > free) {
    auto gib = [](std::size_t b) { return std::to_string(double(b) / double(1 << 30)).substr(0, 5) + " GiB"; };
    throw NotSupportedError("--device cuda: this run needs about " + gib(needed) + " of GPU memory, but only " +
                            gib(free) + " are free (use --device cpu, symmetry or a larger GPU)");
  }
#endif
}

}  // namespace spodvr
