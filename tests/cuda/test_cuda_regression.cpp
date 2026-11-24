// Regression gates 3 and 4 (docs/validation.md) with the GPU backend: the same Fortran references and
// tolerances as test_regression_propagation.cpp. Skipped without a GPU or without the references,
// except the SCM case and cu111_dailight_oblique, which fail without their references.
#include <cstdlib>

#include "cuda/cuda_backend.hpp"
#include "cuda_test_setup.hpp"
#include "regression_gate.hpp"

using namespace spodvr;
using namespace spodvr::testing;

TEST(CudaRegressionPropagation, Gate3ShortRun) {
  SPODVR_REQUIRE_GPU();
  run_gate("hardwall_short", cuda::make_cuda_backend);
}

TEST(CudaRegressionPropagation, Gate3DaiLightShort) {
  SPODVR_REQUIRE_GPU();
  run_gate("cu111_dailight_short", cuda::make_cuda_backend);
}

TEST(CudaRegressionPropagation, Gate3DaiLightShortSymmetric) {
  SPODVR_REQUIRE_GPU();
  run_gate("cu111_dailight_short_sym", cuda::make_cuda_backend);
}

TEST(CudaRegressionPropagation, Gate3DaiLightShortLowZ) {
  SPODVR_REQUIRE_GPU();
  run_gate("cu111_dailight_short_lowz", cuda::make_cuda_backend);
}

// Oblique incidence with J > 0 (test_regression_propagation.cpp Gate3DaiLightOblique); a missing reference file
// fails it.
TEST(CudaRegressionPropagation, Gate3DaiLightOblique) {
  SPODVR_REQUIRE_GPU();
  run_gate("cu111_dailight_oblique", cuda::make_cuda_backend, References::required);
}

// The SCM case (test_regression_propagation.cpp Gate3ScmShort); a missing reference file fails it.
TEST(CudaRegressionPropagation, Gate3ScmShort) {
  SPODVR_REQUIRE_GPU();
  run_gate("scm_cu111_srp48_small", cuda::make_cuda_backend, References::required);
}

TEST(CudaRegressionPropagation, Gate3FullRunLong) {
  SPODVR_REQUIRE_GPU();
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("hardwall_small", cuda::make_cuda_backend);
}

TEST(CudaRegressionPropagation, Gate3DaiLightFullRunLong) {
  SPODVR_REQUIRE_GPU();
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("cu111_dailight_small", cuda::make_cuda_backend);
}
