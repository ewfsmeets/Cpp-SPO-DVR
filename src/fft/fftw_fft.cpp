// FFTW3 implementation of the Fft interface with the Fortran sign/scale conventions.
// Fortran: fftw_2_scs.f90
#include "fft/fftw_fft.hpp"

#include <fftw3.h>
#include <omp.h>

#include <algorithm>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <vector>

// This file is the C-API boundary to FFTW: raw pointers (fftw_plan, vector.data()) are used
// only here.

namespace spodvr {

namespace {

struct PlanDeleter {
  void operator()(fftw_plan plan) const { fftw_destroy_plan(plan); }
};
using PlanPtr = std::unique_ptr<std::remove_pointer_t<fftw_plan>, PlanDeleter>;

int distance_of(const FftLayout& layout, int logical_length) {
  return layout.distance > 0 ? layout.distance : logical_length;
}

// Number of array elements spanned by a layout with the given logical transform length.
std::size_t extent(const FftLayout& layout, int logical_length) {
  const auto dist = static_cast<std::size_t>(distance_of(layout, logical_length));
  return static_cast<std::size_t>(layout.howmany - 1) * dist +
         static_cast<std::size_t>(logical_length - 1) * static_cast<std::size_t>(layout.stride) + 1;
}

void check_layout(const FftLayout& layout) {
  if (layout.length < 1 || layout.howmany < 1 || layout.stride < 1 || layout.distance < 0)
    throw std::invalid_argument("FftwFft: invalid FFT layout");
}

void check_size(std::size_t have, std::size_t need) {
  if (have < need) throw std::invalid_argument("FftwFft: array too small for the FFT layout");
}

fftw_complex* as_fftw(std::vector<std::complex<double>>& v) { return reinterpret_cast<fftw_complex*>(v.data()); }

// One dimension of a transform: length and stride (in complex elements).
struct Dimension {
  int length;
  int stride;
};

// Memory from fftw_malloc (aligned for FFTW's SIMD code), freed with fftw_free.
struct FftwFree {
  void operator()(void* p) const { fftw_free(p); }
};
using FftwBuffer = std::unique_ptr<void, FftwFree>;

// In-place complex transforms (rank 1 or 2, batched) executed directly on the caller's array
// with FFTW's new-array interface; execute() is safe to call concurrently on different arrays.
//   - Two variants of the plan: one for arrays with FFTW's SIMD alignment (fast) and an
//     FFTW_UNALIGNED one for any other array. Which one runs depends only on the address of
//     the array, so results are reproducible (arrays of the propagator are allocated aligned).
//   - oneMKL's FFTW3 interface stores the arrays of a new-array execute in the plan and
//     serialises concurrent calls on one plan; so every OpenMP thread gets its own (identical)
//     copies. A thread without a copy (more threads than when the plan was made) uses the
//     shared copy under a lock.
//   - As in the Fortran wrapper, the scale is applied afterwards and only when it differs from
//     1, and only to the transformed elements (not to padding between transforms).
// Fortran: fftw_2_scs.f90 fftw_2_zzfftm
class FftwComplexPlan final : public ComplexFftPlan {
 public:
  FftwComplexPlan(std::vector<Dimension> dims, Dimension batch, FftSign sign, double scale)
      : dims_(std::move(dims)), batch_(batch), scale_(scale), size_(1 + std::size_t(batch.length - 1) * batch.stride) {
    for (const Dimension& d : dims_) size_ += std::size_t(d.length - 1) * d.stride;

    std::vector<fftw_iodim> iodims;
    for (const Dimension& d : dims_) iodims.push_back({d.length, d.stride, d.stride});
    const fftw_iodim howmany{batch.length, batch.stride, batch.stride};
    const int fftw_sign = sign == FftSign::Minus ? FFTW_FORWARD : FFTW_BACKWARD;
    const int rank = static_cast<int>(iodims.size());

    const FftwBuffer work(fftw_malloc(size_ * sizeof(fftw_complex)));  // only used for planning
    auto* w = static_cast<fftw_complex*>(work.get());
    for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
      aligned_.emplace_back(fftw_plan_guru_dft(rank, iodims.data(), 1, &howmany, w, w, fftw_sign, FFTW_ESTIMATE));
      unaligned_.emplace_back(
          fftw_plan_guru_dft(rank, iodims.data(), 1, &howmany, w, w, fftw_sign, FFTW_ESTIMATE | FFTW_UNALIGNED));
      if (!aligned_.back() || !unaligned_.back()) throw std::runtime_error("FftwFft: complex plan creation failed");
    }
  }

  void execute(std::span<std::complex<double>> data) const override {
    check_size(data.size(), size_);
    auto* p = reinterpret_cast<fftw_complex*>(data.data());
    const bool aligned = fftw_alignment_of(reinterpret_cast<double*>(data.data())) == 0;
    const std::vector<PlanPtr>& plans = aligned ? aligned_ : unaligned_;
    const auto thread = static_cast<std::size_t>(omp_get_thread_num());
    if (thread < plans.size()) {
      fftw_execute_dft(plans[thread].get(), p, p);
    } else {
      const std::lock_guard<std::mutex> lock(shared_plan_mutex_);
      fftw_execute_dft(plans[0].get(), p, p);
    }
    if (scale_ != 1.0) apply_scale(data);
  }

 private:
  void apply_scale(std::span<std::complex<double>> data) const {
    const Dimension outer = dims_.size() == 2 ? dims_[0] : Dimension{1, 0};
    const Dimension inner = dims_.back();
    for (int b = 0; b < batch_.length; ++b)
      for (int o = 0; o < outer.length; ++o)
        for (int i = 0; i < inner.length; ++i) {
          std::complex<double>& v =
              data[std::size_t(b) * batch_.stride + std::size_t(o) * outer.stride + std::size_t(i) * inner.stride];
          v = scale_ * v;
        }
  }

  std::vector<Dimension> dims_;  // outer first
  Dimension batch_;
  double scale_;
  std::size_t size_;
  std::vector<PlanPtr> aligned_;    // one per OpenMP thread
  std::vector<PlanPtr> unaligned_;  // one per OpenMP thread
  mutable std::mutex shared_plan_mutex_;
};

// Fortran: fftw_2_scs.f90 fftw_2_dzfftm_rc_2D. FFTW computes the exp(-i...) transform; for
// sign Plus the result is conjugated (and scaled), as in the Fortran wrapper.
class FftwRealToComplexPlan final : public RealToComplexFftPlan {
 public:
  FftwRealToComplexPlan(const FftLayout& in, const FftLayout& out, FftSign sign, double scale)
      : sign_(sign),
        scale_(scale),
        in_size_(extent(in, in.length)),
        out_size_(extent(out, in.length / 2 + 1)),
        in_buffer_(in_size_),
        out_buffer_(out_size_) {
    int n = in.length;
    plan_.reset(fftw_plan_many_dft_r2c(1, &n, in.howmany, in_buffer_.data(), nullptr, in.stride,
                                       distance_of(in, in.length), as_fftw(out_buffer_), nullptr, out.stride,
                                       distance_of(out, in.length / 2 + 1), FFTW_ESTIMATE));
    if (!plan_) throw std::runtime_error("FftwFft: r2c plan creation failed");
  }

  void execute(std::span<const double> in, std::span<std::complex<double>> out) const override {
    check_size(in.size(), in_size_);
    check_size(out.size(), out_size_);
    std::copy_n(in.begin(), in_size_, in_buffer_.begin());
    fftw_execute(plan_.get());
    if (sign_ == FftSign::Plus) {
      for (auto& v : out_buffer_) v = scale_ != 1.0 ? scale_ * std::conj(v) : std::conj(v);
    } else if (scale_ != 1.0) {
      for (auto& v : out_buffer_) v = scale_ * v;
    }
    std::copy_n(out_buffer_.begin(), out_size_, out.begin());
  }

 private:
  FftSign sign_;
  double scale_;
  std::size_t in_size_;
  std::size_t out_size_;
  mutable std::vector<double> in_buffer_;
  mutable std::vector<std::complex<double>> out_buffer_;
  PlanPtr plan_;
};

// Fortran: fftw_2_scs.f90 fftw_2_zdfftm_cr_2D. FFTW computes the exp(+i...) transform; for
// sign Minus the input is conjugated first, as in the Fortran wrapper.
class FftwComplexToRealPlan final : public ComplexToRealFftPlan {
 public:
  FftwComplexToRealPlan(const FftLayout& in, const FftLayout& out, FftSign sign, double scale)
      : sign_(sign),
        scale_(scale),
        in_size_(extent(in, out.length / 2 + 1)),
        out_size_(extent(out, out.length)),
        in_buffer_(in_size_),
        out_buffer_(out_size_) {
    int n = out.length;
    plan_.reset(fftw_plan_many_dft_c2r(1, &n, out.howmany, as_fftw(in_buffer_), nullptr, in.stride,
                                       distance_of(in, out.length / 2 + 1), out_buffer_.data(), nullptr, out.stride,
                                       distance_of(out, out.length), FFTW_ESTIMATE));
    if (!plan_) throw std::runtime_error("FftwFft: c2r plan creation failed");
  }

  void execute(std::span<const std::complex<double>> in, std::span<double> out) const override {
    check_size(in.size(), in_size_);
    check_size(out.size(), out_size_);
    if (sign_ == FftSign::Minus) {
      std::transform(in.begin(), in.begin() + static_cast<std::ptrdiff_t>(in_size_), in_buffer_.begin(),
                     [](const std::complex<double>& v) { return std::conj(v); });
    } else {
      std::copy_n(in.begin(), in_size_, in_buffer_.begin());
    }
    fftw_execute(plan_.get());  // destroys in_buffer_ (FFTW c2r), never the caller's input
    if (scale_ != 1.0)
      for (auto& v : out_buffer_) v = scale_ * v;
    std::copy_n(out_buffer_.begin(), out_size_, out.begin());
  }

 private:
  FftSign sign_;
  double scale_;
  std::size_t in_size_;
  std::size_t out_size_;
  mutable std::vector<std::complex<double>> in_buffer_;
  mutable std::vector<double> out_buffer_;
  PlanPtr plan_;
};

}  // namespace

std::unique_ptr<ComplexFftPlan> FftwFft::plan_complex(const FftLayout& layout, FftSign sign, double scale) const {
  check_layout(layout);
  return std::make_unique<FftwComplexPlan>(std::vector<Dimension>{{layout.length, layout.stride}},
                                           Dimension{layout.howmany, distance_of(layout, layout.length)}, sign, scale);
}

std::unique_ptr<ComplexFftPlan> FftwFft::plan_complex_2d(const FftLayout2D& layout, FftSign sign, double scale) const {
  if (layout.outer_length < 1 || layout.inner_length < 1 || layout.howmany < 1 || layout.outer_stride < 1 ||
      layout.inner_stride < 1 || layout.distance < 1)
    throw std::invalid_argument("FftwFft: invalid 2D FFT layout");
  return std::make_unique<FftwComplexPlan>(
      std::vector<Dimension>{{layout.outer_length, layout.outer_stride}, {layout.inner_length, layout.inner_stride}},
      Dimension{layout.howmany, layout.distance}, sign, scale);
}

std::unique_ptr<RealToComplexFftPlan> FftwFft::plan_real_to_complex(const FftLayout& in, const FftLayout& out,
                                                                    FftSign sign, double scale) const {
  check_layout(in);
  if (out.howmany != in.howmany) throw std::invalid_argument("FftwFft: howmany mismatch");
  return std::make_unique<FftwRealToComplexPlan>(in, FftLayout{in.length, out.howmany, out.stride, out.distance}, sign,
                                                 scale);
}

std::unique_ptr<ComplexToRealFftPlan> FftwFft::plan_complex_to_real(const FftLayout& in, const FftLayout& out,
                                                                    FftSign sign, double scale) const {
  check_layout(out);
  if (out.howmany != in.howmany) throw std::invalid_argument("FftwFft: howmany mismatch");
  return std::make_unique<FftwComplexToRealPlan>(FftLayout{out.length, in.howmany, in.stride, in.distance}, out, sign,
                                                 scale);
}

}  // namespace spodvr
