// Device memory allocation and host <-> device copies (CUDA runtime).
#include "cuda/device_buffer.hpp"

#include <cuda_runtime_api.h>

#include <string>

#include "cuda/cuda_error.hpp"

namespace spodvr::cuda {

void free_device_memory(void* pointer) noexcept { cudaFree(pointer); }

void* allocate_device_memory(std::size_t bytes) {
  void* p = nullptr;
  if (bytes == 0) return p;
  check_cuda_status(cudaMalloc(&p, bytes), "cudaMalloc of " + std::to_string(bytes) + " bytes");
  return p;
}

void copy_to_device(void* device, const void* host, std::size_t bytes) {
  check_cuda_status(cudaMemcpy(device, host, bytes, cudaMemcpyHostToDevice), "copy to the GPU");
}

void copy_to_host(void* host, const void* device, std::size_t bytes) {
  check_cuda_status(cudaMemcpy(host, device, bytes, cudaMemcpyDeviceToHost), "copy from the GPU");
}

void copy_on_device(void* to, const void* from, std::size_t bytes) {
  check_cuda_status(cudaMemcpy(to, from, bytes, cudaMemcpyDeviceToDevice), "copy on the GPU");
}

void zero_device_memory(void* device, std::size_t bytes) {
  check_cuda_status(cudaMemset(device, 0, bytes), "cudaMemset");
}

}  // namespace spodvr::cuda
