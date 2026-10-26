#pragma once

#include "fft/fft.hpp"

namespace spodvr {

// FFTW3 implementation (double precision, FFTW_ESTIMATE plans like the Fortran build).
// Fortran: fftw_2_scs.f90 (sign/scale handling), ffts.f90
class FftwFft final : public Fft {
 public:
  [[nodiscard]] std::unique_ptr<ComplexFftPlan> plan_complex(const FftLayout& layout, FftSign sign,
                                                             double scale) const override;
  [[nodiscard]] std::unique_ptr<RealToComplexFftPlan> plan_real_to_complex(
      const FftLayout& in, const FftLayout& out, FftSign sign, double scale) const override;
  [[nodiscard]] std::unique_ptr<ComplexToRealFftPlan> plan_complex_to_real(
      const FftLayout& in, const FftLayout& out, FftSign sign, double scale) const override;
};

}  // namespace spodvr
