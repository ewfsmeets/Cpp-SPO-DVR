#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "cuda/device_complex.hpp"
#include "fft/fft.hpp"

namespace spodvr::cuda {

// A batch of in-place complex FFTs of rank 1 or 2 on device memory. Element (o, i) of transform
// b (for rank 1 only i) is at b * distance + (o * inner_length + i) * stride: the dimensions of
// one transform are packed, `stride` apart, which covers every layout of the GPU backend.
struct DeviceFftShape {
  std::vector<std::int64_t> lengths;  // outer first (1 or 2 values)
  std::int64_t stride = 1;
  std::int64_t distance = 0;
  std::int64_t batch = 1;
};

// A cuFFT load callback in LTO-IR (cufftXtSetJITCallback): cuFFT calls the function `name` of
// `lto_fatbin` for every value it reads, with `caller_info` (device memory). The offsets it passes
// count from the start of the data of one execute() call, so a plan with a load callback must run
// its whole batch in one call (no slabs).
struct DeviceFftLoadCallback {
  std::span<const unsigned char> lto_fatbin;
  std::string name;
  void* caller_info = nullptr;  // device pointer
};

// One cuFFT plan (cufftMakePlanMany64). cuFFT does not scale: the scale factors of the CPU plans
// are applied by the kernels around the transforms (folded into their factor tables).
// Large batches are executed in slabs of consecutive transforms, so that cuFFT's work area
// stays below a limit; the work area itself is shared by all plans of a backend
// (set_work_area), since they run one after the other on one stream.
// The sign follows FftSign: Minus = exp(-2 pi i k n / N) = CUFFT_FORWARD.
class DeviceFftPlan {
 public:
  DeviceFftPlan(const DeviceFftShape& shape, FftSign sign, std::size_t max_work_bytes,
                std::optional<DeviceFftLoadCallback> load_callback = std::nullopt);
  ~DeviceFftPlan();
  DeviceFftPlan(const DeviceFftPlan&) = delete;
  DeviceFftPlan& operator=(const DeviceFftPlan&) = delete;

  // Work area needed (bytes); set_work_area must be given at least that much before execute().
  [[nodiscard]] std::size_t work_size() const { return work_size_; }
  void set_work_area(void* work);
  [[nodiscard]] std::int64_t batch_per_call() const { return slab_; }

  void execute(Complex* data) const;

 private:
  struct Plan;  // cufftHandle and its batch
  [[nodiscard]] std::unique_ptr<Plan> make_plan(std::int64_t batch, std::size_t& work) const;

  DeviceFftShape shape_;
  int direction_;
  std::optional<DeviceFftLoadCallback> load_callback_;
  std::int64_t slab_;              // transforms per call
  std::unique_ptr<Plan> full_;     // slab_ transforms
  std::unique_ptr<Plan> rest_;     // the last batch % slab_ transforms (if any)
  std::size_t work_size_ = 0;
};

// The cuFFT version as text (e.g. "cuFFT 11.0.1").
[[nodiscard]] std::string cufft_version();

}  // namespace spodvr::cuda
