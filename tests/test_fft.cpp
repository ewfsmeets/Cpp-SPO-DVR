// Tests: FFT sizes and the FFTW backend against a naive DFT (signs, scales, strides).
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "core/constants.hpp"
#include "fft/fft_sizes.hpp"
#include "core/aligned_vector.hpp"
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
  std::vector<cplx> data((n * howmany));
  for (std::size_t i = 0; i < data.size(); ++i) data[i] = cplx(static_cast<double>(i % 5), static_cast<double>(i % 3));
  std::vector<cplx> y = data;
  fft.plan_complex({n, howmany, howmany, 1}, FftSign::Plus, 1.0)->execute(y);
  for (int b = 0; b < howmany; ++b) {
    std::vector<cplx> x(n);
    for (int i = 0; i < n; ++i) x[i] = data[(b + howmany * i)];
    const auto ref = naive_dft(x, 1.0, 1.0);
    for (int k = 0; k < n; ++k)
      EXPECT_NEAR(std::abs(y[(b + howmany * k)] - ref[k]), 0.0, 1e-13);
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
      EXPECT_NEAR(std::abs(y[k] - ref[k]), 0.0, 1e-12);
  }
}

TEST(Fftw, RoundTripLikeFortranInverse) {
  // Forward r2c with sign + and scale 1, inverse c2r with sign - and scale 1/n
  // (Fortran: InverseRealToComplexFFT) must reproduce the input; the input array is preserved.
  const FftwFft fft;
  const int n = 14;
  const int howmany = 2;
  std::vector<double> x = test_signal((n * howmany), 0.4);
  std::vector<cplx> y(((n / 2 + 1) * howmany));
  const FftLayout real_layout{n, howmany, 1, n};
  const FftLayout half_layout{n, howmany, 1, n / 2 + 1};
  fft.plan_real_to_complex(real_layout, half_layout, FftSign::Plus, 1.0)->execute(x, y);
  const std::vector<cplx> y_copy = y;
  std::vector<double> back(x.size());
  fft.plan_complex_to_real(half_layout, real_layout, FftSign::Minus, 1.0 / n)->execute(y, back);
  for (std::size_t i = 0; i < x.size(); ++i) EXPECT_NEAR(back[i], x[i], 1e-13);
  EXPECT_EQ(y, y_copy);
}

TEST(Fftw, ComplexPaddedLayoutLeavesPaddingAndScalesOnlyTransforms) {
  // Two transforms of length 4 with one padding element each (distance 5), as in the angular
  // DVR blocks. The padding must be neither transformed nor scaled.
  const FftwFft fft;
  const int n = 4;
  std::vector<cplx> data(2 * (n + 1));
  for (std::size_t i = 0; i < data.size(); ++i) data[i] = cplx(0.5 * static_cast<double>(i), 1.0);
  const cplx pad0 = data[4];
  const cplx pad1 = data[9];
  std::vector<cplx> y = data;
  fft.plan_complex({n, 2, 1, n + 1}, FftSign::Minus, 3.0)->execute(y);
  EXPECT_EQ(y[4], pad0);
  EXPECT_EQ(y[9], pad1);
  for (int b = 0; b < 2; ++b) {
    const std::vector<cplx> x(data.begin() + b * (n + 1), data.begin() + b * (n + 1) + n);
    const auto ref = naive_dft(x, -1.0, 3.0);
    for (int k = 0; k < n; ++k) EXPECT_NEAR(std::abs(y[b * (n + 1) + k] - ref[k]), 0.0, 1e-13);
  }
}

TEST(Fftw, ComplexPlanIsSafeToUseConcurrently) {
  // One plan applied by several OpenMP threads to different arrays gives the serial results.
  const FftwFft fft;
  const int n = 16;
  const int arrays = 32;
  const auto plan = fft.plan_complex({n}, FftSign::Plus, 0.25);
  std::vector<std::vector<cplx>> serial(arrays, std::vector<cplx>(n));
  for (int a = 0; a < arrays; ++a)
    for (int i = 0; i < n; ++i) serial[a][i] = cplx(std::cos(0.1 * a + i), std::sin(0.3 * i - a));
  std::vector<std::vector<cplx>> parallel = serial;
  for (auto& v : serial) plan->execute(v);
#pragma omp parallel for
  for (int a = 0; a < arrays; ++a) plan->execute(parallel[a]);
  EXPECT_EQ(serial, parallel);
}

TEST(Fftw, Complex2dBatchedStridedMatchesNaiveDft) {
  // 3 interleaved transforms of 4 (outer) x 6 (inner) points: element (o, i) of transform b at
  // b + 3 * (i + 6 * o), i.e. the layout of the gathered (Z, r) blocks of the propagator.
  // Checked on an aligned array and on an unaligned one (offset by one element).
  const FftwFft fft;
  const int howmany = 3, n_outer = 4, n_inner = 6;
  const FftLayout2D layout{n_outer, howmany * n_inner, n_inner, howmany, howmany, 1};
  const auto plan = fft.plan_complex_2d(layout, FftSign::Minus, 0.5);
  for (const int offset : {0, 1}) {
    ComplexVector storage(offset + howmany * n_outer * n_inner);
    const std::span<cplx> data = std::span(storage).subspan(offset);
    for (std::size_t i = 0; i < data.size(); ++i) data[i] = cplx(std::cos(0.7 * double(i)), 0.1 * double(i % 5));
    const std::vector<cplx> x(data.begin(), data.end());
    plan->execute(data);
    for (int b = 0; b < howmany; ++b)
      for (int ko = 0; ko < n_outer; ++ko)
        for (int ki = 0; ki < n_inner; ++ki) {
          cplx ref(0.0, 0.0);
          for (int o = 0; o < n_outer; ++o)
            for (int i = 0; i < n_inner; ++i)
              ref += x[b + howmany * (i + n_inner * o)] *
                     std::exp(cplx(0.0, -2.0 * M_PI * (double(ko * o) / n_outer + double(ki * i) / n_inner)));
          EXPECT_NEAR(std::abs(data[b + howmany * (ki + n_inner * ko)] - 0.5 * ref), 0.0, 1e-13) << offset;
        }
  }
}
