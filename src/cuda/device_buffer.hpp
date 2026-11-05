#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <span>

#include "cuda/device_complex.hpp"

namespace spodvr::cuda {

// Raw allocation and copies (cuda_runtime calls, in device_buffer.cpp).
void free_device_memory(void* pointer) noexcept;
[[nodiscard]] void* allocate_device_memory(std::size_t bytes);
void copy_to_device(void* device, const void* host, std::size_t bytes);
void copy_to_host(void* host, const void* device, std::size_t bytes);
void copy_on_device(void* to, const void* from, std::size_t bytes);
void zero_device_memory(void* device, std::size_t bytes);

// An array in GPU memory (RAII, like std::vector but without element access from the host).
// Device pointers are plain pointers by nature: they are the documented exception to the
// pointer rule inside src/cuda/ and never leave it.
template <class T>
class DeviceBuffer {
 public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t size)
      : data_(static_cast<T*>(allocate_device_memory(size * sizeof(T)))), size_(size) {}
  // Allocates and copies the host values.
  explicit DeviceBuffer(std::span<const T> host) : DeviceBuffer(host.size()) { upload(host); }

  [[nodiscard]] T* data() { return data_.get(); }
  [[nodiscard]] const T* data() const { return data_.get(); }
  [[nodiscard]] std::size_t size() const { return size_; }
  [[nodiscard]] std::size_t bytes() const { return size_ * sizeof(T); }

  // host -> device / device -> host, starting at element `offset` of this buffer.
  void upload(std::span<const T> host, std::size_t offset = 0) {
    copy_to_device(data() + offset, host.data(), host.size_bytes());
  }
  void download(std::span<T> host, std::size_t offset = 0) const {
    copy_to_host(host.data(), data() + offset, host.size_bytes());
  }
  void set_zero() { zero_device_memory(data(), bytes()); }

 private:
  struct Free {
    void operator()(T* p) const noexcept { free_device_memory(p); }
  };
  std::unique_ptr<T, Free> data_;
  std::size_t size_ = 0;
};

// std::complex<double> host arrays as device Complex arrays (same layout, see Complex).
[[nodiscard]] inline std::span<const Complex> as_device_complex(std::span<const std::complex<double>> values) {
  return {reinterpret_cast<const Complex*>(values.data()), values.size()};
}
[[nodiscard]] inline std::span<Complex> as_device_complex(std::span<std::complex<double>> values) {
  return {reinterpret_cast<Complex*>(values.data()), values.size()};
}

// A device copy of a host table of complex numbers.
[[nodiscard]] inline DeviceBuffer<Complex> upload_complex(std::span<const std::complex<double>> values) {
  return DeviceBuffer<Complex>(as_device_complex(values));
}

}  // namespace spodvr::cuda
