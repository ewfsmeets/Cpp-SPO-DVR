#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

namespace spodvr::cuda {

// A failed CUDA runtime or cuFFT call (or a kernel launch).
class CudaError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Throws CudaError naming `what` if a CUDA runtime call failed. The status is passed as an int
// so that this header needs no CUDA header (cudaError_t / cufftResult values).
void check_cuda_status(int status, std::string_view what);
void check_cufft_status(int status, std::string_view what);

// Checks the launch of the kernel(s) just started (cudaGetLastError). With
// SPODVR_CUDA_SYNC_DEBUG it also waits for them, so that an error in a kernel is reported here
// and not at some later call.
void check_kernel_launch(std::string_view kernel);

}  // namespace spodvr::cuda
