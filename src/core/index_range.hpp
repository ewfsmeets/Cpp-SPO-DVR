#pragma once

namespace spodvr {

// Inclusive 0-based index range [first, last].
struct IndexRange {
  int first;
  int last;
  [[nodiscard]] int size() const { return last - first + 1; }
};

}  // namespace spodvr
