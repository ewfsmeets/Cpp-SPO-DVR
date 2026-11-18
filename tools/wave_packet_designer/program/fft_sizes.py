"""Valid FFT sizes: 2^h 3^i 5^j 7^k 11^m with h >= 1, i <= 2 and j, k, m <= 1 (so every size is even).

Shadows: src/fft/fft_sizes.cpp (nearest_fft_size, is_valid_fft_size)

Examples: 112, 160, 198, 210, 240, 256 are valid; 196, 200, 216, 250, 270 and 1 are not.
"""

LIMIT = 1 << 30  # sizes at or above this are not considered


def nearest_fft_size(n: int, multiple: int = 1) -> int:
    """The smallest size >= n of the form multiple' * 2^h 3^i 5^j 7^k 11^m (the C++ loop, same limits).

    multiple' is `multiple`, halved once if it is even, so the size is a multiple of `multiple` but need not
    be valid itself: nearest_fft_size(50, 6) = 54 = 3 * 2 * 3^2, and 3^3 is not allowed.
    """
    factor = multiple // 2 if multiple % 2 == 0 else multiple  # the factor 2 is always present (h >= 1)
    best = LIMIT
    for h in range(1, 26):
        for i in range(3):
            for j in range(2):
                for k in range(2):
                    for m in range(2):
                        size = factor * 2**h * 3**i * 5**j * 7**k * 11**m
                        if size <= 0 or size >= LIMIT:
                            continue
                        if n <= size <= best:
                            best = size
    return best


def is_valid_fft_size(n: int) -> bool:
    """True if n is itself a valid size."""
    return nearest_fft_size(n) == n
