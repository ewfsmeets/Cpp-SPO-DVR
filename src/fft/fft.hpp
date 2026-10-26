#pragma once

#include <complex>
#include <memory>
#include <span>

namespace spodvr {

// Sign of the exponent: Plus computes y_k = scale * sum_n x_n exp(+2 pi i k n / N),
// Minus uses exp(-2 pi i k n / N). The Fortran FFT module uses the SCS-library convention in
// which its "forward" transforms have exponent sign Plus (ffts.f90 / fftw_2_scs.f90).
enum class FftSign { Minus, Plus };

// A batch of `howmany` 1D transforms of `length` points; element n of transform b is at
// offset b * distance + n * stride (in elements of the array it applies to).
struct FftLayout {
  int length;
  int howmany = 1;
  int stride = 1;
  int distance = 0;  // 0 means: length (contiguous transforms)
};

// A batch of `howmany` 2D transforms of outer_length x inner_length points; element (o, i) of
// transform b is at offset b * distance + o * outer_stride + i * inner_stride.
// (Maps onto cufftPlanMany for the layouts used: inner_stride = howmany-interleave,
// outer_stride = inner_length * inner_stride.)
struct FftLayout2D {
  int outer_length;
  int outer_stride;
  int inner_length;
  int inner_stride;
  int howmany = 1;
  int distance = 1;
};

// In-place complex transform. execute() may be called concurrently on different arrays.
class ComplexFftPlan {
 public:
  virtual ~ComplexFftPlan() = default;
  virtual void execute(std::span<std::complex<double>> data) const = 0;
};

// Real -> complex (half spectrum, length/2+1 values per transform).
class RealToComplexFftPlan {
 public:
  virtual ~RealToComplexFftPlan() = default;
  virtual void execute(std::span<const double> in, std::span<std::complex<double>> out) const = 0;
};

// Complex (half spectrum) -> real. The imaginary parts of the k = 0 and k = N/2 inputs are
// ignored (they must vanish for a Hermitian spectrum). The input is not modified.
class ComplexToRealFftPlan {
 public:
  virtual ~ComplexToRealFftPlan() = default;
  virtual void execute(std::span<const std::complex<double>> in, std::span<double> out) const = 0;
};

// Creates FFT plans. Layouts: `in` describes the input array, `out` the output array (for the
// half spectrum the logical length is length/2+1, `length` always refers to the real length).
class Fft {
 public:
  virtual ~Fft() = default;
  [[nodiscard]] virtual std::unique_ptr<ComplexFftPlan> plan_complex(const FftLayout& layout,
                                                                     FftSign sign, double scale) const = 0;
  [[nodiscard]] virtual std::unique_ptr<ComplexFftPlan> plan_complex_2d(const FftLayout2D& layout, FftSign sign,
                                                                        double scale) const = 0;
  [[nodiscard]] virtual std::unique_ptr<RealToComplexFftPlan> plan_real_to_complex(
      const FftLayout& in, const FftLayout& out, FftSign sign, double scale) const = 0;
  [[nodiscard]] virtual std::unique_ptr<ComplexToRealFftPlan> plan_complex_to_real(
      const FftLayout& in, const FftLayout& out, FftSign sign, double scale) const = 0;
};

}  // namespace spodvr
