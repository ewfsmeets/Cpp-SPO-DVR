#pragma once

#include <complex>
#include <cstddef>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace spodvr {

// std::vector storage for the large arrays:
//   - aligned to 64 bytes (a cache line; enough for the SIMD code paths of FFTW and oneMKL,
//     which are faster on aligned data);
//   - NOT initialised on allocation: `AlignedVector<double> v(n)` leaves the n values
//     undefined (double and std::complex<double> are implicit-lifetime types, so this is valid
//     C++20). The owner must write every element before reading it, and does so in parallel for
//     the large arrays, with the same static OpenMP schedule as the loops that use them: on
//     multi-socket (NUMA) machines a memory page is placed on the socket of the thread that
//     first writes it ("first touch"), so each thread then works on local memory.
// The pointers exist only inside the standard allocator interface used by std::vector (the one
// place besides the C-API boundaries where they appear; see CLAUDE.md Pointers).
template <class T>
struct AlignedAllocator {
  using value_type = T;
  static constexpr std::align_val_t alignment{64};

  AlignedAllocator() = default;
  template <class U>
  explicit AlignedAllocator(const AlignedAllocator<U>&) noexcept {}

  [[nodiscard]] T* allocate(std::size_t n) { return static_cast<T*>(::operator new(n * sizeof(T), alignment)); }
  void deallocate(T* p, std::size_t) noexcept { ::operator delete(p, alignment); }

  // Default construction leaves the element uninitialised (see above); copies are constructed.
  // Debug builds with SPODVR_POISON_UNINITIALISED fill NaN instead, so that any read of an
  // unwritten element shows up as NaN in the results (cmake -DSPODVR_POISON_UNINITIALISED=ON).
  template <class U>
  void construct([[maybe_unused]] U* p) noexcept {
#ifdef SPODVR_POISON_UNINITIALISED
    ::new (static_cast<void*>(p)) U(poison<U>());
#endif
  }
  template <class U, class Arg, class... Args>
  void construct(U* p, Arg&& arg, Args&&... args) {
    ::new (static_cast<void*>(p)) U(std::forward<Arg>(arg), std::forward<Args>(args)...);
  }

  friend bool operator==(const AlignedAllocator&, const AlignedAllocator&) { return true; }

 private:
  template <class U>
  static U poison() {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    if constexpr (std::is_same_v<U, std::complex<double>>)
      return U(nan, nan);
    else
      return U(nan);
  }
};

template <class T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;
using ComplexVector = AlignedVector<std::complex<double>>;

}  // namespace spodvr
