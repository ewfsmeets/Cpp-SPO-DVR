"""Tests (ctest WavePacketDesigner.FftSizes): the valid FFT sizes and the nearest size with a multiple
(program/fft_sizes.py)."""

from program.fft_sizes import is_valid_fft_size, nearest_fft_size
from tests.helpers import check, check_equal


def test_sizes_of_the_input_guide() -> None:
    """The examples of inputs/README.md section 9: valid and invalid sizes, and 1 is invalid (h >= 1)."""
    for n in (112, 160, 198, 210, 240, 256):
        check(is_valid_fft_size(n), f"{n} should be valid")
    for n in (1, 196, 200, 216, 250, 270):
        check(not is_valid_fft_size(n), f"{n} should be invalid")


def test_nearest_size_with_a_multiple() -> None:
    """nPhi for MJMax 28: 2*28+1 = 57 -> 60, also as a multiple of 6 (C6v) and 4 (C4v)."""
    check_equal(nearest_fft_size(57), 60, "nearest_fft_size(57)")
    check_equal(nearest_fft_size(57, 6), 60, "nearest_fft_size(57, 6)")
    check_equal(nearest_fft_size(57, 4), 60, "nearest_fft_size(57, 4)")
    check_equal(nearest_fft_size(9, 2), 10, "nearest_fft_size(9, 2)")
    check_equal(nearest_fft_size(101), 110, "nearest_fft_size(101)")  # 102, 104, 108 are invalid
