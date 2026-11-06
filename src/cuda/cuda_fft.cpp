// cuFFT plans on device memory for the GPU backend.
#include "cuda/cuda_fft.hpp"

#include <cufft.h>
#include <cufftXt.h>

#include <stdexcept>
#include <string>
#include <utility>

#include "cuda/cuda_error.hpp"

namespace spodvr::cuda {

struct DeviceFftPlan::Plan {
  cufftHandle handle = 0;
  std::int64_t batch = 0;
  ~Plan() {
    if (handle != 0) cufftDestroy(handle);
  }
};

std::unique_ptr<DeviceFftPlan::Plan> DeviceFftPlan::make_plan(std::int64_t batch, std::size_t& work) const {
  auto plan = std::make_unique<Plan>();
  plan->batch = batch;
  check_cufft_status(cufftCreate(&plan->handle), "cufftCreate");
  check_cufft_status(cufftSetAutoAllocation(plan->handle, 0), "cufftSetAutoAllocation");
  if (load_callback_) {  // must be set between cufftCreate and the plan initialisation
    void* info = load_callback_->caller_info;
    check_cufft_status(
        cufftXtSetJITCallback(plan->handle, load_callback_->name.c_str(), load_callback_->lto_fatbin.data(),
                              load_callback_->lto_fatbin.size(), CUFFT_CB_LD_COMPLEX_DOUBLE, &info),
        "cufftXtSetJITCallback");
  }
  std::vector<long long> n(shape_.lengths.begin(), shape_.lengths.end());
  const int rank = static_cast<int>(n.size());
  const int status = cufftMakePlanMany64(plan->handle, rank, n.data(), n.data(), shape_.stride, shape_.distance,
                                         n.data(), shape_.stride, shape_.distance, CUFFT_Z2Z, batch, &work);
  if (load_callback_ && status != CUFFT_SUCCESS)
    throw CudaError("cuFFT could not link its load callback (cufftMakePlanMany64 returned " + std::to_string(status) +
                    "). It needs libnvJitLink and libnvrtc of the CUDA release of cuFFT on the library path.");
  check_cufft_status(status, "cufftMakePlanMany64");
  return plan;
}

DeviceFftPlan::DeviceFftPlan(const DeviceFftShape& shape, FftSign sign, std::size_t max_work_bytes,
                             std::optional<DeviceFftLoadCallback> load_callback)
    : shape_(shape),
      direction_(sign == FftSign::Minus ? CUFFT_FORWARD : CUFFT_INVERSE),
      load_callback_(std::move(load_callback)),
      slab_(shape.batch) {
  // Halve the slab until cuFFT's work area fits (cuFFT may need as much work memory as data).
  std::size_t work = 0;
  full_ = make_plan(slab_, work);
  while (work > max_work_bytes && slab_ > 1) {
    slab_ = (slab_ + 1) / 2;
    full_ = make_plan(slab_, work);
  }
  work_size_ = work;
  if (load_callback_ && slab_ != shape.batch)
    throw std::invalid_argument("DeviceFftPlan: a plan with a load callback must run its batch in one call");
  if (const std::int64_t rest = shape.batch % slab_; rest != 0) {
    std::size_t rest_work = 0;
    rest_ = make_plan(rest, rest_work);
    work_size_ = std::max(work_size_, rest_work);
  }
}

DeviceFftPlan::~DeviceFftPlan() = default;

void DeviceFftPlan::set_work_area(void* work) {
  check_cufft_status(cufftSetWorkArea(full_->handle, work), "cufftSetWorkArea");
  if (rest_) check_cufft_status(cufftSetWorkArea(rest_->handle, work), "cufftSetWorkArea");
}

void DeviceFftPlan::execute(Complex* data) const {
  std::int64_t first = 0;
  for (; first + slab_ <= shape_.batch; first += slab_) {
    auto* p = reinterpret_cast<cufftDoubleComplex*>(data + first * shape_.distance);
    check_cufft_status(cufftExecZ2Z(full_->handle, p, p, direction_), "cufftExecZ2Z");
  }
  if (rest_) {
    auto* p = reinterpret_cast<cufftDoubleComplex*>(data + first * shape_.distance);
    check_cufft_status(cufftExecZ2Z(rest_->handle, p, p, direction_), "cufftExecZ2Z");
  }
}

std::string cufft_version() {
  int v = 0;
  cufftGetVersion(&v);
  return "cuFFT " + std::to_string(v / 1000) + "." + std::to_string(v / 100 % 10) + "." + std::to_string(v % 100);
}

}  // namespace spodvr::cuda
