// Error checking for the CUDA runtime, cuFFT and kernel launches.
#include "cuda/cuda_error.hpp"

#include <cuda_runtime_api.h>
#include <cufft.h>

namespace spodvr::cuda {

void check_cuda_status(int status, std::string_view what) {
  if (status == cudaSuccess) return;
  throw CudaError("CUDA error in " + std::string(what) + ": " + cudaGetErrorString(static_cast<cudaError_t>(status)));
}

void check_cufft_status(int status, std::string_view what) {
  if (status == CUFFT_SUCCESS) return;
  throw CudaError("cuFFT error " + std::to_string(status) + " in " + std::string(what));
}

void check_kernel_launch(std::string_view kernel) {
  check_cuda_status(cudaGetLastError(), kernel);
#ifdef SPODVR_CUDA_SYNC_DEBUG
  check_cuda_status(cudaDeviceSynchronize(), kernel);
#endif
}

}  // namespace spodvr::cuda
