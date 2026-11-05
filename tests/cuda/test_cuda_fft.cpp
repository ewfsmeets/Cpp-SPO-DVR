// Tests: cuFFT plans of the GPU backend (DeviceFftPlan) against the CPU FFT (FftwFft) for the
// layouts the backend uses, including execution in slabs.
#include <vector>

#include "cuda/cuda_fft.hpp"
#include "cuda/device_buffer.hpp"
#include "cuda_test_setup.hpp"
#include "fft/fftw_fft.hpp"

using namespace spodvr;
using namespace spodvr::cuda;
using namespace spodvr::testing;

namespace {

using cplx = std::complex<double>;

// Transforms `values` with the GPU plan and the CPU plan; returns the relative difference.
double compare(const DeviceFftShape& gpu_shape, const FftLayout2D* cpu_2d, const FftLayout* cpu_1d, FftSign sign,
               std::vector<cplx> values, std::size_t max_work = std::size_t(1) << 30) {
  DeviceFftPlan plan(gpu_shape, sign, max_work);
  DeviceBuffer<std::byte> work(std::max<std::size_t>(plan.work_size(), 1));
  plan.set_work_area(work.data());
  DeviceBuffer<Complex> data(as_device_complex(std::span<const cplx>(values)));
  plan.execute(data.data());
  std::vector<cplx> gpu(values.size());
  data.download(as_device_complex(std::span<cplx>(gpu)));

  const FftwFft fft;
  const auto cpu = cpu_2d ? fft.plan_complex_2d(*cpu_2d, sign, 1.0) : fft.plan_complex(*cpu_1d, sign, 1.0);
  cpu->execute(values);
  return relative_difference(values, gpu);
}

}  // namespace

// T_Z + T_r: contiguous 2D (Z, r) blocks, batch at distance n_zr, both signs.
TEST(CudaFft, ContiguousBatchOf2DTransformsMatchesTheCpu) {
  SPODVR_REQUIRE_GPU();
  const int nz = 24;
  const int nr = 10;
  const int batch = 7;
  std::vector<cplx> v(std::size_t(nz) * nr * batch);
  fill_pseudo_random(v);
  const DeviceFftShape gpu{{nz, nr}, 1, nz * nr, batch};
  const FftLayout2D cpu{nz, nr, nr, 1, batch, nz * nr};
  EXPECT_LT(compare(gpu, &cpu, nullptr, FftSign::Minus, v), 1e-15);
  EXPECT_LT(compare(gpu, &cpu, nullptr, FftSign::Plus, v), 1e-15);
}

// T_xy: 2D (x, y) with the points n_k apart and the batches at distance 1.
TEST(CudaFft, StridedBatchOf2DTransformsMatchesTheCpu) {
  SPODVR_REQUIRE_GPU();
  const int nx = 6;
  const int ny = 4;
  const int n_k = 15;
  std::vector<cplx> v(std::size_t(nx) * ny * n_k);
  fill_pseudo_random(v, 7);
  const DeviceFftShape gpu{{nx, ny}, n_k, 1, n_k};
  const FftLayout2D cpu{nx, n_k * ny, ny, n_k, n_k, 1};
  EXPECT_LT(compare(gpu, &cpu, nullptr, FftSign::Minus, v), 1e-15);
  EXPECT_LT(compare(gpu, &cpu, nullptr, FftSign::Plus, v), 1e-15);
}

// The phi transforms of the angular transform: 1D, values n apart, batches at distance 1; a tiny
// work-area limit forces execution in slabs (with a remainder).
TEST(CudaFft, StridedBatchOf1DTransformsMatchesTheCpuAlsoInSlabs) {
  SPODVR_REQUIRE_GPU();
  const int n_phi = 12;
  const int n = 1001;  // columns x theta
  std::vector<cplx> v(std::size_t(n_phi) * n);
  fill_pseudo_random(v, 3);
  const DeviceFftShape gpu{{n_phi}, n, 1, n};
  const FftLayout cpu{n_phi, n, n, 1};
  EXPECT_LT(compare(gpu, nullptr, &cpu, FftSign::Plus, v), 1e-15);
  DeviceFftPlan small(gpu, FftSign::Plus, 1);  // work limit 1 byte: smallest slabs cuFFT allows
  EXPECT_LE(small.batch_per_call(), n);
  EXPECT_LT(compare(gpu, nullptr, &cpu, FftSign::Plus, v, 1), 1e-15);
}
