// Tests: the decision of SPODVR_REQUIRE_GPU (gpu_check) for what the CUDA runtime can report. Pure
// function: runs without a GPU.
#include "cuda_test_setup.hpp"

using spodvr::testing::gpu_check;
using spodvr::testing::GpuCheck;

// A machine without a GPU skips (no driver, no device); a driver with any other runtime error fails,
// so that a broken CUDA set-up is not reported as skipped tests.
TEST(CudaTestSetup, GpuCheckDecision) {
  const int driver = 13040;  // a driver is installed (cudaDriverGetVersion != 0)

  EXPECT_EQ(gpu_check(driver, cudaSuccess, 1), GpuCheck::run);
  EXPECT_EQ(gpu_check(driver, cudaSuccess, 2), GpuCheck::run);

  // No GPU: no driver installed (whatever cudaGetDeviceCount says), no device, or none counted.
  EXPECT_EQ(gpu_check(0, cudaErrorInsufficientDriver, 0), GpuCheck::skip);
  EXPECT_EQ(gpu_check(0, cudaErrorNoDevice, 0), GpuCheck::skip);
  EXPECT_EQ(gpu_check(0, cudaSuccess, 0), GpuCheck::skip);
  EXPECT_EQ(gpu_check(driver, cudaErrorNoDevice, 0), GpuCheck::skip);
  EXPECT_EQ(gpu_check(driver, cudaSuccess, 0), GpuCheck::skip);

  // A driver is present but the runtime cannot use it: a broken set-up, not a missing GPU.
  EXPECT_EQ(gpu_check(driver, cudaErrorInsufficientDriver, 0), GpuCheck::fail);
  EXPECT_EQ(gpu_check(driver, cudaErrorInitializationError, 0), GpuCheck::fail);
  EXPECT_EQ(gpu_check(driver, cudaErrorSystemDriverMismatch, 0), GpuCheck::fail);
  EXPECT_EQ(gpu_check(driver, cudaErrorUnknown, 1), GpuCheck::fail);  // the count does not matter
}
