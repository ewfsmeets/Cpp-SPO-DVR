#pragma once

#include <cstddef>

namespace spodvr {

// Flat-index helpers for multi-dimensional arrays stored in Fortran memory order (first index
// fastest). All indices are 0-based.
struct Index2D {
  std::size_t n0;  // extent of the fastest index
  std::size_t n1;

  [[nodiscard]] constexpr std::size_t operator()(std::size_t i0, std::size_t i1) const { return i0 + n0 * i1; }
  [[nodiscard]] constexpr std::size_t size() const { return n0 * n1; }
};

struct Index3D {
  std::size_t n0;  // extent of the fastest index
  std::size_t n1;
  std::size_t n2;

  [[nodiscard]] constexpr std::size_t operator()(std::size_t i0, std::size_t i1, std::size_t i2) const {
    return i0 + n0 * (i1 + n1 * i2);
  }
  [[nodiscard]] constexpr std::size_t size() const { return n0 * n1 * n2; }
};

}  // namespace spodvr
