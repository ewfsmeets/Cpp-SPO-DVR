// Tests: the flux surface (psi and d psi / dr at r = r_c on the planes below Z_asymptotic): the CPU
// function against the column-by-column reference flux_r_derivative, and the CPU backend's seam.
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <bit>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "backend_setup.hpp"
#include "propagation/flux_surface.hpp"
#include "propagation/propagation_backend.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

// A synthetic wave: one full tile of the derivative sum plus a partial one of 96 blocks per (Z, r)
// point (1120 blocks with 1024-block tiles), 12 points in r, 5 planes in Z of which the surface uses
// the first 4.
constexpr int partial_tile_blocks = 96;
constexpr int n_blocks = flux_surface_tile_blocks + partial_tile_blocks;
constexpr int n_rot = 8;
static_assert(partial_tile_blocks < flux_surface_tile_blocks && n_blocks % n_rot == 0);
constexpr int n_xy = n_blocks / n_rot;
constexpr int nr = 12;
constexpr int nz = 5;
constexpr int n_planes = 4;
constexpr double dr = 0.15;  // bohr

const cplx nan_value(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

// Distinct values of varying magnitude and phase at every point.
void fill_pseudo_random(std::span<cplx> values) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    const double x = static_cast<double>(i);
    values[i] = cplx(std::sin(0.37 * x + 0.1), std::cos(1.31 * x)) * (1.0 + 0.5 * std::sin(0.011 * x));
  }
}

ScatteringWave synthetic_wave() {
  ScatteringWave psi(n_rot, n_xy, nr * nz);
  fill_pseudo_random(psi.values);
  return psi;
}

// A surface whose every element is NaN, as an unwritten element of a reused slot would be garbage.
FluxSurfaceWave nan_surface(int blocks, int planes) {
  FluxSurfaceWave out(blocks, planes);
  std::ranges::fill(out.psi, nan_value);
  std::ranges::fill(out.derivative, nan_value);
  return out;
}

// The (Z, r) point izr = ir + nr iz (the index of ScatteringWave::block).
int zr_point(int ir, int iz, int n_r) { return ir + n_r * iz; }

// psi(block, ir, iz) for ir = 0 .. nr - 1: the column the Fortran sums over.
std::vector<cplx> r_column(const ScatteringWave& psi, int n_r, int block, int iz) {
  std::vector<cplx> column(n_r);
  for (int ir = 0; ir < n_r; ++ir) column[ir] = psi.block(zr_point(ir, iz, n_r))[block];
  return column;
}

// flux_surface.cpp has the build's floating-point flags, flux_r_derivative always strict ones: bitwise
// equal in strict builds. SPODVR_FAST_MATH (__FAST_MATH__ here: the same flags) may fuse each multiply-add
// into an FMA, moving each of the nr roundings by <= eps sum |D psi|: bound 2 nr eps sum |D psi| (build-fast,
// 2026-11-12: 35546 of 40560 values differed, the largest by 0.06 of the bound). That path may miss a NaN
// (-ffinite-math-only); the strict path's EXPECT_EQ fails on it.
void expect_reference_derivative(cplx actual, std::span<const double> weights, std::span<const cplx> column) {
  const cplx expected = flux_r_derivative(weights, column);
#ifdef __FAST_MATH__
  double magnitude = 0.0;
  for (std::size_t ir = 0; ir < weights.size(); ++ir) magnitude = magnitude + std::abs(weights[ir] * column[ir]);
  EXPECT_LE(std::abs(actual - expected), 2.0 * double(weights.size()) * DBL_EPSILON * magnitude);
#else
  EXPECT_EQ(actual.real(), expected.real());
  EXPECT_EQ(actual.imag(), expected.imag());
#endif
}

// Every (block, iz) of out against the column reference of psi: psi_c copied, d psi / dr as
// flux_r_derivative sums it. Counts the elements checked.
int expect_reference_surface(const FluxSurfaceWave& out, const ScatteringWave& psi, int ir_flux,
                             std::span<const double> weights) {
  const int n_r = static_cast<int>(weights.size());
  int checked = 0;
  for (int iz = 0; iz < out.n_planes(); ++iz)
    for (int block = 0; block < out.n_blocks(); ++block) {
      SCOPED_TRACE(::testing::Message() << "block " << block << ", iz " << iz);
      const std::vector<cplx> column = r_column(psi, n_r, block, iz);
      const cplx psi_c = out.psi[out.index(block, iz)];
      EXPECT_EQ(psi_c.real(), column[ir_flux].real());
      EXPECT_EQ(psi_c.imag(), column[ir_flux].imag());
      expect_reference_derivative(out.derivative[out.index(block, iz)], weights, column);
      if (::testing::Test::HasFailure()) return checked;  // one failing element is enough output
      checked = checked + 1;
    }
  return checked;
}

// Number of doubles whose bit patterns differ: equal patterns also mean that neither array kept a
// NaN where the other has a number, in every build (no floating-point comparison involved).
int bit_differences(std::span<const cplx> a, std::span<const cplx> b) {
  int differences = 0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::bit_cast<std::uint64_t>(a[i].real()) != std::bit_cast<std::uint64_t>(b[i].real())) ++differences;
    if (std::bit_cast<std::uint64_t>(a[i].imag()) != std::bit_cast<std::uint64_t>(b[i].imag())) ++differences;
  }
  return differences;
}

}  // namespace

// psi_c = psi(:, ir_c, iz) and d psi / dr = sum_ir D(ir) psi(:, ir, iz) summed as flux_r_derivative sums
// one column (bitwise in strict builds), on every block of both derivative tiles and every plane below
// nz_potential, for r_c at either end and inside the r grid. Every element is written (NaN before).
TEST(FluxSurface, EqualsTheColumnReferenceAtEveryBlockAndPlane) {
  const ScatteringWave psi = synthetic_wave();
  for (const int ir_flux : {0, 5, nr - 1}) {
    SCOPED_TRACE(::testing::Message() << "ir_flux " << ir_flux);
    const std::vector<double> weights = flux_derivative_weights(nr, dr, ir_flux);
    FluxSurfaceWave out = nan_surface(n_blocks, n_planes);
    flux_surface(psi, ir_flux, weights, out);
    EXPECT_EQ(expect_reference_surface(out, psi, ir_flux, weights), n_blocks * n_planes);
  }
}

// The planes beyond out.n_planes() are not read: a NaN on plane n_planes (Z >= Z_asymptotic) leaves
// the surface finite and unchanged.
TEST(FluxSurface, ReadsOnlyThePlanesOfTheSurface) {
  ScatteringWave psi = synthetic_wave();
  const std::vector<double> weights = flux_derivative_weights(nr, dr, 5);
  FluxSurfaceWave before(n_blocks, n_planes);
  flux_surface(psi, 5, weights, before);
  for (int ir = 0; ir < nr; ++ir) std::ranges::fill(psi.block(zr_point(ir, n_planes, nr)), nan_value);
  FluxSurfaceWave after(n_blocks, n_planes);
  flux_surface(psi, 5, weights, after);
  EXPECT_EQ(bit_differences(before.psi, after.psi), 0);
  EXPECT_EQ(bit_differences(before.derivative, after.derivative), 0);
}

// The zero term D(ir_c) psi(ir_c) is added, as in the Fortran: a non-finite psi at r_c makes
// d psi / dr NaN (0 * inf) at that block only.
TEST(FluxSurface, AddsTheZeroTermAtTheFluxSurface) {
#ifdef __FAST_MATH__
  GTEST_SKIP() << "-ffinite-math-only: inf and NaN have no defined arithmetic in this build";
#else
  ScatteringWave psi = synthetic_wave();
  const int ir_flux = 5;
  const int iz = 2;
  const int block = flux_surface_tile_blocks + 6;  // in the second tile
  const std::vector<double> weights = flux_derivative_weights(nr, dr, ir_flux);
  ASSERT_EQ(weights[ir_flux], 0.0);
  psi.block(zr_point(ir_flux, iz, nr))[block] = cplx(std::numeric_limits<double>::infinity(), 0.0);
  FluxSurfaceWave out(n_blocks, n_planes);
  flux_surface(psi, ir_flux, weights, out);
  EXPECT_TRUE(std::isnan(out.derivative[out.index(block, iz)].real()));
  EXPECT_FALSE(std::isnan(out.derivative[out.index(block - 1, iz)].real()));
  EXPECT_FALSE(std::isnan(out.derivative[out.index(block, iz - 1)].real()));
#endif
}

// Each element is computed by one thread in a fixed order: 1 and 4 threads give the same bits (any
// build). The 1-thread surface starts as NaN, the 4-thread one as zero: equal bits also show that
// every element was written.
TEST(FluxSurface, DoesNotDependOnTheThreadCount) {
  const ScatteringWave psi = synthetic_wave();
  const std::vector<double> weights = flux_derivative_weights(nr, dr, 5);
  const int threads = omp_get_max_threads();
  FluxSurfaceWave serial = nan_surface(n_blocks, n_planes);
  omp_set_num_threads(1);
  flux_surface(psi, 5, weights, serial);
  FluxSurfaceWave parallel(n_blocks, n_planes);
  omp_set_num_threads(4);
  flux_surface(psi, 5, weights, parallel);
  omp_set_num_threads(threads);
  EXPECT_EQ(bit_differences(serial.psi, parallel.psi), 0);
  EXPECT_EQ(bit_differences(serial.derivative, parallel.derivative), 0);
}

TEST(FluxSurface, RejectsShapesThatDoNotFit) {
  const ScatteringWave psi = synthetic_wave();
  const std::vector<double> weights = flux_derivative_weights(nr, dr, 5);
  FluxSurfaceWave fits(n_blocks, n_planes);
  EXPECT_NO_THROW(flux_surface(psi, 5, weights, fits));
  FluxSurfaceWave other_blocks(n_blocks - 1, n_planes);
  EXPECT_THROW(flux_surface(psi, 5, weights, other_blocks), std::invalid_argument);
  FluxSurfaceWave too_many_planes(n_blocks, nz + 1);
  EXPECT_THROW(flux_surface(psi, 5, weights, too_many_planes), std::invalid_argument);
  EXPECT_THROW(flux_surface(psi, nr, weights, fits), std::invalid_argument);
  EXPECT_THROW(flux_surface(psi, -1, weights, fits), std::invalid_argument);
  const std::vector<double> seven_weights = flux_derivative_weights(7, dr, 5);  // 60 (Z, r) points: no whole planes
  EXPECT_THROW(flux_surface(psi, 5, seven_weights, fits), std::invalid_argument);
  EXPECT_THROW(flux_surface(psi, 5, std::span<const double>(), fits), std::invalid_argument);
}

// ------------------------------------------------------------------------------------------------
// The seam: CpuBackend::flux_surface through PropagationBackend
// ------------------------------------------------------------------------------------------------

namespace {

using spodvr::testing::BackendSetup;

// hardwall_small with the flux analysis at r_c = 2.0 bohr (ir 8).
nlohmann::json hardwall_small_flux_keys() { return spodvr::testing::flux_analysis_keys("hardwall_small", 2.0, 0.15); }

}  // namespace

// The backend reads its own psi at r_c = rFluxAnalysis with the weights D of that r_c: equal to the
// column reference built here from the grid alone (bitwise in strict builds).
TEST(FluxSurface, CpuBackendReadsItsWaveAtRFluxAnalysis) {
  const BackendSetup s("hardwall_small", hardwall_small_flux_keys());
  ASSERT_TRUE(s.grids.ir_flux_analysis.has_value());
  const int ir_flux = *s.grids.ir_flux_analysis;
  ASSERT_EQ(ir_flux, 8);
  const std::unique_ptr<PropagationBackend> backend = s.cpu();
  backend->read_wave_packet([](WavePacket& psi) { fill_pseudo_random(psi.scattering.values); });

  const ScatteringWave& psi = backend->wave_packet_on_host().scattering;
  FluxSurfaceWave out = nan_surface(psi.block_size(), s.grids.nz_potential);
  std::as_const(*backend).flux_surface(out);
  const std::vector<double> weights = flux_derivative_weights(s.config.input.nr, s.grids.r.spacing(), ir_flux);
  EXPECT_EQ(expect_reference_surface(out, psi, ir_flux, weights), psi.block_size() * s.grids.nz_potential);
}

TEST(FluxSurface, CpuBackendRefusesWithoutFluxAnalysis) {
  const BackendSetup s("hardwall_small");
  const std::unique_ptr<PropagationBackend> backend = s.cpu();
  FluxSurfaceWave out(backend->wave_packet_on_host().scattering.block_size(), s.grids.nz_potential);
  EXPECT_THROW(backend->flux_surface(out), std::logic_error);
}

// The seam's contract: a slot of nz_potential planes of n_rot nx ny blocks, nothing else (also not
// fewer planes, which the free function would accept).
TEST(FluxSurface, CpuBackendRefusesASlotOfAnotherShape) {
  const BackendSetup s("hardwall_small", hardwall_small_flux_keys());
  const std::unique_ptr<PropagationBackend> backend = s.cpu();
  const int blocks = backend->wave_packet_on_host().scattering.block_size();
  const int planes = s.grids.nz_potential;
  FluxSurfaceWave fits(blocks, planes);
  EXPECT_NO_THROW(backend->flux_surface(fits));
  FluxSurfaceWave fewer_planes(blocks, planes - 1);
  EXPECT_THROW(backend->flux_surface(fewer_planes), std::logic_error);
  FluxSurfaceWave more_planes(blocks, planes + 1);
  EXPECT_THROW(backend->flux_surface(more_planes), std::logic_error);
  FluxSurfaceWave other_blocks(blocks + 1, planes);
  EXPECT_THROW(backend->flux_surface(other_blocks), std::logic_error);
}
