#pragma once

#include <complex>
#include <cstddef>
#include <new>
#include <vector>

namespace spodvr {

// std::vector storage aligned to 64 bytes (a cache line; enough for the SIMD code paths of FFTW
// and oneMKL, which are faster on aligned data). The pointers exist only inside the standard
// allocator interface used by std::vector (the one place besides the C-API boundaries where
// they appear; see CLAUDE.md Pointers).
template <class T>
struct AlignedAllocator {
  using value_type = T;
  static constexpr std::align_val_t alignment{64};

  AlignedAllocator() = default;
  template <class U>
  explicit AlignedAllocator(const AlignedAllocator<U>&) noexcept {}

  [[nodiscard]] T* allocate(std::size_t n) { return static_cast<T*>(::operator new(n * sizeof(T), alignment)); }
  void deallocate(T* p, std::size_t) noexcept { ::operator delete(p, alignment); }

  friend bool operator==(const AlignedAllocator&, const AlignedAllocator&) { return true; }
};

using ComplexVector = std::vector<std::complex<double>, AlignedAllocator<std::complex<double>>>;

}  // namespace spodvr
