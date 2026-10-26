#pragma once

#include <cstddef>
#include <optional>

namespace spodvr {

// Physical memory of this machine in bytes (empty if the system does not report it).
[[nodiscard]] std::optional<std::size_t> physical_memory_bytes();

}  // namespace spodvr
