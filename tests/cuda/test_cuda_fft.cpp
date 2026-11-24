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

// The phi transforms of the angular transform: 1D, values n apart, batches at distance 1. cuFFT
// needs no work area for them (CUDA 13.3, 1001 columns: none for the 22 lengths from 10 to 8192
// tried), so they run in one call whatever the limit. A length of 16384 needs 16 bytes of work per
// value: a limit of 10 such transforms forces execution in slabs (of 8 transforms) with a remainder.
TEST(CudaFft, StridedBatchOf1DTransformsMatchesTheCpuAlsoInSlabs) {
  SPODVR_REQUIRE_GPU();
  const int n_phi = 12;
  const int n = 1001;  // columns x theta
  std::vector<cplx> v(std::size_t(n_phi) * n);
  fill_pseudo_random(v, 3);
  const DeviceFftShape gpu{{n_phi}, n, 1, n};
  const FftLayout cpu{n_phi, n, n, 1};
  EXPECT_LT(compare(gpu, nullptr, &cpu, FftSign::Plus, v), 1e-15);

  const int n_long = 16384;
  const int n_columns = 61;
  std::vector<cplx> w(std::size_t(n_long) * n_columns);
  fill_pseudo_random(w, 4);
  const DeviceFftShape gpu_long{{n_long}, n_columns, 1, n_columns};
  const FftLayout cpu_long{n_long, n_columns, n_columns, 1};
  const std::size_t limit = std::size_t(10) * n_long * sizeof(cplx);
  DeviceFftPlan sliced(gpu_long, FftSign::Plus, limit);
  EXPECT_LT(sliced.batch_per_call(), n_columns);              // more than one slab
  EXPECT_NE(n_columns % sliced.batch_per_call(), 0);          // and a remainder slab
  EXPECT_LT(compare(gpu_long, nullptr, &cpu_long, FftSign::Plus, w), 1e-15);         // one call
  EXPECT_LT(compare(gpu_long, nullptr, &cpu_long, FftSign::Plus, w, limit), 1e-15);  // in slabs
}
