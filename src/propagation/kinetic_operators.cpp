// Kinetic energy propagators on the scattering grid (FFT in Z and r over the whole wave
// function, FFT in x and y per (Z, r) point). Fortran: KineticOp.f90, SplitOperator.f90
#include "propagation/kinetic_operators.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

const std::vector<cplx>& factors_for(const StepFactors& f, TimeStep step) {
  return step == TimeStep::Full ? f.full : f.half;
}

}  // namespace

// Fortran: KineticOp.f90 SetupFFTs. The FFT along one coordinate is planned for a fixed work
// unit (all (irot, iy) values of one (ix, ir), resp. (ix, iz), resp. (Z, r) point), so results
// do not depend on how the units are distributed over threads.
ScatteringKineticOperator::ScatteringKineticOperator(const ScatteringShape& shape, const PropagationTables& tables,
                                                     const Fft& fft)
    : shape_(shape), tables_(tables) {
  const int n_roty = shape.n_rot * shape.ny;
  const int n_rotxy = n_roty * shape.nx;
  const int z_stride = n_rotxy * shape.nr;

  const FftLayout z_layout{shape.nz, n_roty, z_stride, 1};
  z_to_momentum_ = fft.plan_complex(z_layout, FftSign::Minus, tables.fft_scale_z);
  z_to_coordinate_ = fft.plan_complex(z_layout, FftSign::Plus, 1.0 / (shape.nz * tables.fft_scale_z));

  const FftLayout r_layout{shape.nr, n_roty, n_rotxy, 1};
  r_to_momentum_ = fft.plan_complex(r_layout, FftSign::Minus, tables.fft_scale_r);
  r_to_coordinate_ = fft.plan_complex(r_layout, FftSign::Plus, 1.0 / (shape.nr * tables.fft_scale_r));

  const FftLayout x_layout{shape.nx, n_roty, n_roty, 1};
  x_to_momentum_ = fft.plan_complex(x_layout, FftSign::Minus, tables.fft_scale_x);
  x_to_coordinate_ = fft.plan_complex(x_layout, FftSign::Plus, 1.0 / (shape.nx * tables.fft_scale_x));

  const FftLayout y_layout{shape.ny, shape.n_rot, shape.n_rot, 1};
  y_to_momentum_ = fft.plan_complex(y_layout, FftSign::Minus, tables.fft_scale_y);
  y_to_coordinate_ = fft.plan_complex(y_layout, FftSign::Plus, 1.0 / (shape.ny * tables.fft_scale_y));
}

// Fortran: KineticOp.f90 PerformEXPKinOpZScat, then PerformEXPKinOpRScat
void ScatteringKineticOperator::apply_zr(ScatteringWave& psi, TimeStep step) const {
  const int n_roty = shape_.n_rot * shape_.ny;
  const int n_rotxy = n_roty * shape_.nx;
  const int z_stride = n_rotxy * shape_.nr;
  const std::span<cplx> values(psi.values);

  // Along Z: one unit per (ix, ir), holding all (irot, iy) values at every iz.
  const std::vector<cplx>& ez = factors_for(tables_.z_kinetic, step);
#pragma omp parallel for schedule(static)
  for (int unit = 0; unit < shape_.nx * shape_.nr; ++unit) {
    const std::span<cplx> column = values.subspan(n_roty * unit);
    z_to_momentum_->execute(column);
    for (int iz = 0; iz < shape_.nz; ++iz)
      for (int e = 0; e < n_roty; ++e) column[e + z_stride * iz] = ez[iz] * column[e + z_stride * iz];
    z_to_coordinate_->execute(column);
  }

  // Along r: one unit per (ix, iz).
  const std::vector<cplx>& er = factors_for(tables_.r_kinetic, step);
#pragma omp parallel for schedule(static)
  for (int unit = 0; unit < shape_.nx * shape_.nz; ++unit) {
    const int ix = unit % shape_.nx;
    const int iz = unit / shape_.nx;
    const std::span<cplx> column = values.subspan(n_roty * ix + z_stride * iz);
    r_to_momentum_->execute(column);
    for (int ir = 0; ir < shape_.nr; ++ir)
      for (int e = 0; e < n_roty; ++e) column[e + n_rotxy * ir] = er[ir] * column[e + n_rotxy * ir];
    r_to_coordinate_->execute(column);
  }
}

void ScatteringKineticOperator::apply_xy(ScatteringWave& psi, TimeStep step) const {
  const std::vector<cplx>& exy = factors_for(tables_.xy_kinetic, step);
#pragma omp parallel for schedule(static)
  for (int izr = 0; izr < psi.n_blocks(); ++izr) apply_xy_block(psi.block(izr), exy);
}

// One (Z, r) point: shift to the moving frame, FFT x then y, multiply by exp(-i T_xy dt),
// back to coordinates, shift back. Fortran: SplitOperator.f90 PerformEXPKineticOp
void ScatteringKineticOperator::apply_xy_block(std::span<cplx> block, const std::vector<cplx>& factors) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;

  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot) block[irot + n_rot * ixy] = block[irot + n_rot * ixy] * tables_.shift_forward[ixy];

  // Fortran: TransToMomenta (FFTAlongX, FFTAlongY); no transform for a single point.
  if (shape_.nx > 1) x_to_momentum_->execute(block);
  if (shape_.ny > 1)
    for (int ix = 0; ix < shape_.nx; ++ix) y_to_momentum_->execute(block.subspan(n_rot * shape_.ny * ix));

  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot) block[irot + n_rot * ixy] = factors[ixy] * block[irot + n_rot * ixy];

  // Fortran: TransToCoords (FFTAlongX, FFTAlongY)
  if (shape_.nx > 1) x_to_coordinate_->execute(block);
  if (shape_.ny > 1)
    for (int ix = 0; ix < shape_.nx; ++ix) y_to_coordinate_->execute(block.subspan(n_rot * shape_.ny * ix));

  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot) block[irot + n_rot * ixy] = block[irot + n_rot * ixy] * tables_.shift_backward[ixy];
}

}  // namespace spodvr
