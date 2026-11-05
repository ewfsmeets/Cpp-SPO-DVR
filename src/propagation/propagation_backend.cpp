// Device names of the propagation backends (--device).
#include "propagation/propagation_backend.hpp"

#include "core/error.hpp"

namespace spodvr {

Device parse_device(std::string_view name) {
  if (name == "cpu") return Device::Cpu;
  if (name == "cuda") return Device::Cuda;
  throw InputError("--device", "unknown device '" + std::string(name) + "' (expected cpu or cuda)");
}

std::string to_string(Device device) { return device == Device::Cpu ? "cpu" : "cuda"; }

}  // namespace spodvr
