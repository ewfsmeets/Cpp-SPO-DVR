// Kinetic energy propagators on the scattering grid: 2D FFT in (r, Z) on gathered work units,
// 2D FFT in (x, y) per (Z, r) point. Fortran: KineticOp.f90, SplitOperator.f90
#include "propagation/kinetic_operators.hpp"

#include <algorithm>

#include "core/aligned_vector.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

const std::vector<cplx>& factors_for(const StepFactors& f, TimeStep step) {
  return step == TimeStep::Full ? f.full : f.half;
}

bool is_one(const cplx& c) { return c.real() == 1.0 && c.imag() == 0.0; }

}  // namespace

int zr_chunk_size(int n_roty) {
  for (int c = std::min(4, n_roty); c > 1; --c)
    if (n_roty % c == 0) return c;
  return 1;
}

int xy_group_points(const ScatteringShape& shape) {
  constexpr std::size_t max_bytes = std::size_t(1) << 20;
  const std::size_t block_bytes = sizeof(cplx) * shape.n_rot * shape.nx * shape.ny;
  for (int g = shape.nr; g > 1; --g)
    if (shape.nr % g == 0 && g * block_bytes <= max_bytes) return g;
  return 1;
}

// Fortran: KineticOp.f90 SetupFFTs, SetupScaleFactors
ScatteringKineticOperator::ScatteringKineticOperator(const ScatteringShape& shape, const PropagationTables& tables,
                                                     const Fft& fft)
    : shape_(shape),
      tables_(tables),
      chunk_(zr_chunk_size(shape.n_rot * shape.ny)),
      shift_(!std::all_of(tables.shift_forward.begin(), tables.shift_forward.end(), is_one)),
      group_points_(xy_group_points(shape)) {
  // (r, Z) unit buffer: (e, ir, iz) with e = 0 .. chunk-1 fastest.
  const FftLayout2D zr_layout{shape.nz, chunk_ * shape.nr, shape.nr, chunk_, chunk_, 1};
  const double scale_z = tables.fft_scale_z;
  const double scale_r = tables.fft_scale_r;
  zr_to_momentum_ = fft.plan_complex_2d(zr_layout, FftSign::Minus, scale_z * scale_r);
  zr_to_coordinate_ =
      fft.plan_complex_2d(zr_layout, FftSign::Plus, (1.0 / (shape.nz * scale_z)) * (1.0 / (shape.nr * scale_r)));

  // (x, y) for one group, gathered as (k, iy, ix) with k = irot + n_rot * (point in group)
  // fastest. A single grid point in x or y is not transformed (nor scaled) in the Fortran code.
  const int k_count = shape.n_rot * group_points_;
  const FftLayout2D xy_layout{shape.nx, k_count * shape.ny, shape.ny, k_count, k_count, 1};
  const double scale_x = shape.nx > 1 ? tables.fft_scale_x : 1.0;
  const double scale_y = shape.ny > 1 ? tables.fft_scale_y : 1.0;
  xy_to_momentum_ = fft.plan_complex_2d(xy_layout, FftSign::Minus, scale_x * scale_y);
  xy_to_coordinate_ =
      fft.plan_complex_2d(xy_layout, FftSign::Plus, (1.0 / (shape.nx * scale_x)) * (1.0 / (shape.ny * scale_y)));
}

// Every work unit (a chunk of adjacent (irot, iy) values at one ix) is gathered into a
// contiguous buffer, transformed in (r, Z), multiplied by exp(-i (T_Z + T_r) dt), transformed
// back and scattered. Fixed units, so the result does not depend on the number of threads.
// Fortran: KineticOp.f90 PerformEXPKinOpZScat, PerformEXPKinOpRScat
void ScatteringKineticOperator::apply_zr(ScatteringWave& psi, TimeStep step) const {
  const int n_roty = shape_.n_rot * shape_.ny;
  const std::size_t r_stride = std::size_t(n_roty) * shape_.nx;
  const std::size_t z_stride = r_stride * shape_.nr;
  const int nr = shape_.nr;
  const int nz = shape_.nz;
  const int c = chunk_;
  const int chunks_per_x = n_roty / c;
  const std::vector<cplx>& factors = factors_for(tables_.zr_kinetic, step);

#pragma omp parallel
  {
    ComplexVector unit(std::size_t(c) * nr * nz);
#pragma omp for schedule(static)
    for (int u = 0; u < chunks_per_x * shape_.nx; ++u) {
      const int ix = u / chunks_per_x;
      const std::size_t first = std::size_t(c) * (u % chunks_per_x) + std::size_t(n_roty) * ix;

      for (int iz = 0; iz < nz; ++iz)
        for (int ir = 0; ir < nr; ++ir)
          for (int e = 0; e < c; ++e)
            unit[e + c * (ir + nr * iz)] = psi.values[first + e + r_stride * ir + z_stride * iz];

      zr_to_momentum_->execute(unit);
      for (int k = 0; k < nr * nz; ++k)
        for (int e = 0; e < c; ++e) unit[e + c * k] = factors[k] * unit[e + c * k];
      zr_to_coordinate_->execute(unit);

      for (int iz = 0; iz < nz; ++iz)
        for (int ir = 0; ir < nr; ++ir)
          for (int e = 0; e < c; ++e)
            psi.values[first + e + r_stride * ir + z_stride * iz] = unit[e + c * (ir + nr * iz)];
    }
  }
}

void ScatteringKineticOperator::apply_xy(ScatteringWave& psi, TimeStep step) const {
  const std::size_t group_size = std::size_t(psi.block_size()) * group_points_;
  const int n_groups = psi.n_blocks() / group_points_;
#pragma omp parallel
  {
    ComplexVector work(group_work_size());
#pragma omp for schedule(static)
    for (int g = 0; g < n_groups; ++g)
      apply_xy_group(std::span(psi.values).subspan(group_size * g, group_size), step, work);
  }
}

std::size_t ScatteringKineticOperator::group_work_size() const {
  return std::size_t(shape_.n_rot) * shape_.nx * shape_.ny * group_points_;
}

// A group of points: gather their blocks (irot, ixy) into the work array as (k, ixy) with
// k = irot + n_rot * (point in group), shift to the moving frame, FFT (x, y), multiply by
// exp(-i T_xy dt), back to coordinates, shift back, scatter. The shift is exactly 1 without
// parallel momentum and then skipped. Fortran: SplitOperator.f90 PerformEXPKineticOp
void ScatteringKineticOperator::apply_xy_group(std::span<cplx> group, TimeStep step, std::span<cplx> work) const {
  const std::vector<cplx>& factors = factors_for(tables_.xy_kinetic, step);
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  const int k_count = n_rot * group_points_;

  for (int ir = 0; ir < group_points_; ++ir)
    for (int ixy = 0; ixy < n_xy; ++ixy)
      for (int irot = 0; irot < n_rot; ++irot)
        work[irot + n_rot * ir + k_count * ixy] = group[irot + n_rot * (ixy + n_xy * ir)];

  if (shift_)
    for (int ixy = 0; ixy < n_xy; ++ixy)
      for (int k = 0; k < k_count; ++k) work[k + k_count * ixy] = work[k + k_count * ixy] * tables_.shift_forward[ixy];

  xy_to_momentum_->execute(work);  // Fortran: TransToMomenta
  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int k = 0; k < k_count; ++k) work[k + k_count * ixy] = factors[ixy] * work[k + k_count * ixy];
  xy_to_coordinate_->execute(work);  // Fortran: TransToCoords

  if (shift_)
    for (int ixy = 0; ixy < n_xy; ++ixy)
      for (int k = 0; k < k_count; ++k) work[k + k_count * ixy] = work[k + k_count * ixy] * tables_.shift_backward[ixy];

  for (int ir = 0; ir < group_points_; ++ir)
    for (int ixy = 0; ixy < n_xy; ++ixy)
      for (int irot = 0; irot < n_rot; ++irot)
        group[irot + n_rot * (ixy + n_xy * ir)] = work[irot + n_rot * ir + k_count * ixy];
}

}  // namespace spodvr
