#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "core/constants.hpp"
#include "fft/fft_sizes.hpp"
#include "fft/fftw_fft.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

std::vector<cplx> naive_dft(const std::vector<cplx>& x, double sign, double scale) {
  const std::size_t n = x.size();
  std::vector<cplx> y(n);
  for (std::size_t k = 0; k < n; ++k) {
    cplx s = 0.0;
    for (std::size_t j = 0; j < n; ++j)
      s += x[j] * std::polar(1.0, sign * Constants::two_pi * static_cast<double>(k * j) / static_cast<double>(n));
    y[k] = scale * s;
  }
  return y;
}

std::vector<double> test_signal(std::size_t n, double phase) {
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; ++i) x[i] = std::sin(0.7 * static_cast<double>(i) + phase) + 0.1 * static_cast<double>(i);
  return x;
}

}  // namespace

TEST(FftSizes, NearestSizes) {
  EXPECT_EQ(nearest_fft_size(9), 10);
  EXPECT_EQ(nearest_fft_size(56), 56);
  EXPECT_EQ(nearest_fft_size(57), 60);
  EXPECT_EQ(nearest_fft_size(1), 2);
  EXPECT_EQ(nearest_fft_size(176), 176);
  EXPECT_EQ(nearest_fft_size(280), 280);
  EXPECT_EQ(nearest_fft_size(27), 28);   // 27 = 3^3 is not allowed (i <= 2), 28 = 4*7
  EXPECT_EQ(nearest_fft_size(49), 56);   // 50 = 2*5^2 and 54 = 2*3^3 are not allowed
  EXPECT_EQ(nearest_fft_size(9, 4), 12); // dividable 4: multiples of 4
  EXPECT_EQ(nearest_fft_size(9, 3), 12); // multiples of 3, even
  EXPECT_TRUE(is_valid_fft_size(4));
  EXPECT_FALSE(is_valid_fft_size(3));
  EXPECT_FALSE(is_valid_fft_size(5));
}

TEST(Fftw, ComplexMatchesNaiveDftBothSigns) {
  const FftwFft fft;
  const std::size_t n = 12;
  std::vector<cplx> x(n);
  for (std::size_t i = 0; i < n; ++i) x[i] = cplx(std::cos(0.3 * static_cast<double>(i)), 0.2 * static_cast<double>(i));
  for (const auto& [sign, s] : {std::pair{FftSign::Minus, -1.0}, std::pair{FftSign::Plus, 1.0}}) {
    std::vector<cplx> y = x;
    fft.plan_complex({static_cast<int>(n)}, sign, 0.5)->execute(y);
    const auto ref = naive_dft(x, s, 0.5);
    for (std::size_t k = 0; k < n; ++k) EXPECT_NEAR(std::abs(y[k] - ref[k]), 0.0, 1e-13);
  }
}

TEST(Fftw, StridedBatchedComplex) {
  const FftwFft fft;
  // 3 transforms of length 8 interleaved: element n of transform b at b + 3*n.
  const int n = 8;
  const int howmany = 3;
  std::vector<cplx> data(static_cast<std::size_t>(n * howmany));
  for (std::size_t i = 0; i < data.size(); ++i) data[i] = cplx(static_cast<double>(i % 5), static_cast<double>(i % 3));
  std::vector<cplx> y = data;
  fft.plan_complex({n, howmany, howmany, 1}, FftSign::Plus, 1.0)->execute(y);
  for (int b = 0; b < howmany; ++b) {
    std::vector<cplx> x(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) x[static_cast<std::size_t>(i)] = data[static_cast<std::size_t>(b + howmany * i)];
    const auto ref = naive_dft(x, 1.0, 1.0);
    for (int k = 0; k < n; ++k)
      EXPECT_NEAR(std::abs(y[static_cast<std::size_t>(b + howmany * k)] - ref[static_cast<std::size_t>(k)]), 0.0, 1e-13);
  }
}

TEST(Fftw, RealToComplexSignConventions) {
  const FftwFft fft;
  const int n = 10;
  const auto x = test_signal(n, 0.1);
  std::vector<cplx> xc(x.begin(), x.end());
  for (const auto& [sign, s] : {std::pair{FftSign::Minus, -1.0}, std::pair{FftSign::Plus, 1.0}}) {
    std::vector<cplx> y(n / 2 + 1);
    fft.plan_real_to_complex({n}, {n}, sign, 2.0)->execute(x, y);
    const auto ref = naive_dft(xc, s, 2.0);
    for (int k = 0; k <= n / 2; ++k)
      EXPECT_NEAR(std::abs(y[static_cast<std::size_t>(k)] - ref[static_cast<std::size_t>(k)]), 0.0, 1e-12);
  }
}

TEST(Fftw, RoundTripLikeFortranInverse) {
  // Forward r2c with sign + and scale 1, inverse c2r with sign - and scale 1/n
  // (Fortran: InverseRealToComplexFFT) must reproduce the input; the input array is preserved.
  const FftwFft fft;
  const int n = 14;
  const int howmany = 2;
  std::vector<double> x = test_signal(static_cast<std::size_t>(n * howmany), 0.4);
  std::vector<cplx> y(static_cast<std::size_t>((n / 2 + 1) * howmany));
  const FftLayout real_layout{n, howmany, 1, n};
  const FftLayout half_layout{n, howmany, 1, n / 2 + 1};
  fft.plan_real_to_complex(real_layout, half_layout, FftSign::Plus, 1.0)->execute(x, y);
  const std::vector<cplx> y_copy = y;
  std::vector<double> back(x.size());
  fft.plan_complex_to_real(half_layout, real_layout, FftSign::Minus, 1.0 / n)->execute(y, back);
  for (std::size_t i = 0; i < x.size(); ++i) EXPECT_NEAR(back[i], x[i], 1e-13);
  EXPECT_EQ(y, y_copy);
}
