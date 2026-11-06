// The embedded LTO-IR of the cuFFT load callback (fft_callbacks.cu). The header is generated at
// build time (CMakeLists.txt: nvcc -dc -fatbin, then cmake/embed_binary.cmake).
#include "cuda/fft_callbacks.hpp"

#include "fft_callbacks_fatbin.h"  // spodvr_fft_callbacks_fatbin[]

namespace spodvr::cuda {

std::span<const unsigned char> exp_potential_callback_fatbin() { return spodvr_fft_callbacks_fatbin; }

}  // namespace spodvr::cuda
