// Choice of the propagation backend (--device). The only place that knows whether the CUDA
// backend was built.
#include "app/backend_factory.hpp"

#include "core/error.hpp"
#include "propagation/cpu_backend.hpp"

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
  throw NotSupportedError("--device cuda: this program was built without CUDA support (cmake -DSPODVR_USE_CUDA=ON)");
}

}  // namespace spodvr
