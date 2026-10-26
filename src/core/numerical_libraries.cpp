// Set-up and description of the external numerical libraries (oneMKL or system BLAS/LAPACK,
// oneMKL's FFTW3 interface or FFTW).
#include "core/numerical_libraries.hpp"

#include <cstdlib>
#include <mutex>

#ifdef SPODVR_HAVE_MKL
#include <mkl_service.h>
#endif
#ifdef SPODVR_FFT_FROM_MKL
#include <fftw3_mkl.h>
#include <omp.h>
#endif

namespace spodvr {

namespace {

#ifdef SPODVR_HAVE_MKL
std::string cnr_branch_name(int branch) {
  switch (branch) {
    case MKL_CBWR_OFF: return "OFF";
    case MKL_CBWR_BRANCH_OFF: return "OFF";
    case MKL_CBWR_AUTO: return "AUTO";
    case MKL_CBWR_COMPATIBLE: return "COMPATIBLE";
    case MKL_CBWR_SSE2: return "SSE2";
    case MKL_CBWR_SSE4_2: return "SSE4_2";
    case MKL_CBWR_AVX2: return "AVX2";
    case MKL_CBWR_AVX512: return "AVX512";
    case MKL_CBWR_AVX512_E1: return "AVX512_E1";
    case MKL_CBWR_AVX10: return "AVX10";
    default: return std::to_string(branch);
  }
}
#endif

}  // namespace

void configure_numerical_libraries() {
  static std::once_flag once;
  std::call_once(once, [] {
#ifdef SPODVR_HAVE_MKL
    if (std::getenv("MKL_CBWR") == nullptr) mkl_cbwr_set(MKL_CBWR_AUTO);
#endif
#ifdef SPODVR_FFT_FROM_MKL
    fftw3_mkl.number_of_user_threads = omp_get_max_threads();  // plans executed from all threads
#endif
  });
}

std::string numerical_libraries_summary() {
  configure_numerical_libraries();
#ifdef SPODVR_HAVE_MKL
  const std::string source = std::getenv("MKL_CBWR") != nullptr ? "from MKL_CBWR" : "default";
  std::string text = std::string("BLAS/LAPACK: oneMKL ") + SPODVR_MKL_VERSION + " (sequential), reproducibility branch " + cnr_branch_name(mkl_cbwr_get(MKL_CBWR_BRANCH)) +
                     " (" + source + ")";
#else
  std::string text = "BLAS/LAPACK: system libraries";
#endif
#ifdef SPODVR_FFT_FROM_MKL
  text += "; FFT: oneMKL (FFTW3 interface)";
#else
  text += "; FFT: FFTW";
#endif
  return text;
}

}  // namespace spodvr
