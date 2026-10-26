#pragma once

namespace spodvr {

// Smallest allowed FFT length >= n: dividable' * 2^h * 3^i * 5^j * 7^k * 11^m with h = 1..25,
// i <= 2, j, k, m <= 1, where dividable' = dividable, halved once if it is even (so the
// length is always even). Fortran: ffts.f90 GetNearestFFTSize
[[nodiscard]] int nearest_fft_size(int n, int dividable = 1);

// True if n is an allowed FFT length (nearest_fft_size(n, 1) == n).
[[nodiscard]] bool is_valid_fft_size(int n);

}  // namespace spodvr
