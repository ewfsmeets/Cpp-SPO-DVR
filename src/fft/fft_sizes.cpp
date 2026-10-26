#include "fft/fft_sizes.hpp"

#include <cstdint>

namespace spodvr {

// Fortran: ffts.f90 GetNearestFFTSize (same loop structure and limits; 2^30 sentinel)
int nearest_fft_size(int n, int dividable) {
  constexpr std::int64_t limit = std::int64_t{1} << 30;
  std::int64_t fft_dividable = dividable;
  if (fft_dividable % 2 == 0) fft_dividable /= 2;

  std::int64_t best = limit;
  for (int h = 1; h <= 25; ++h) {
    for (int i = 0; i <= 2; ++i) {
      for (int j = 0; j <= 1; ++j) {
        for (int k = 0; k <= 1; ++k) {
          for (int m = 0; m <= 1; ++m) {
            std::int64_t size = fft_dividable * (std::int64_t{1} << h);
            for (int c = 0; c < i; ++c) size *= 3;
            if (j == 1) size *= 5;
            if (k == 1) size *= 7;
            if (m == 1) size *= 11;
            if (size <= 0 || size >= limit) continue;
            if (size >= n && size <= best) best = size;
          }
        }
      }
    }
  }
  return static_cast<int>(best);
}

bool is_valid_fft_size(int n) { return nearest_fft_size(n, 1) == n; }

}  // namespace spodvr
