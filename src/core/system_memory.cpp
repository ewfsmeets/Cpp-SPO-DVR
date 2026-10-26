// Physical memory of the machine (POSIX sysconf).
#include "core/system_memory.hpp"

#include <unistd.h>

namespace spodvr {

std::optional<std::size_t> physical_memory_bytes() {
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGE_SIZE);
  if (pages <= 0 || page_size <= 0) return std::nullopt;
  return std::size_t(pages) * std::size_t(page_size);
}

}  // namespace spodvr
