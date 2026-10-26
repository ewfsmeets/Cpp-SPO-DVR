// Balint-Kurti coefficients: projection of the wave packet at the analysis surface onto the
// asymptotic channels. Fortran: Analyzor.f90 BalintKurtiCoeffs
#include "analysis/balint_kurti.hpp"

#include <cmath>

namespace spodvr {

namespace {

using cplx = std::complex<double>;

}  // namespace

BalintKurtiAnalyzer::BalintKurtiAnalyzer(const Config& config, const SpatialGrids& grids, const AnalysisBasis& basis,
                                         const GasPhase& gas, const PropagationTables& tables, const Fft& fft)
    : basis_(basis),
      gas_(gas),
      tables_(tables),
      nx_(config.input.nx),
      ny_(config.input.ny),
      nr_(config.input.nr),
      iz_analysis_(grids.iz_analysis),
      j_in_(config.input.j_in),
      dr_(config.input.dr) {
  const DerivedParameters& d = config.derived;
  xy_scale_ = std::sqrt(d.dx * d.dy * std::sin(d.system.skewing_angle) / (nx_ * ny_));
  // One (x, y) array in the layout iy + ny*ix, to momentum space (exponent sign -, scale 1).
  const FftLayout2D layout{.outer_length = nx_, .outer_stride = ny_, .inner_length = ny_, .inner_stride = 1};
  xy_to_momentum_ = fft.plan_complex_2d(layout, FftSign::Minus, 1.0);
}

BkCoefficients BalintKurtiAnalyzer::coefficients(const WavePacket& psi, const AnalysisSurface& surface) const {
  const AnalysisBasis& a = basis_;
  const int n_xy = nx_ * ny_;
  const int n_rot = a.n_rotation();
  const int n_rot_basis = static_cast<int>(psi.scattering.index.n0);  // stride of ixy in a (Z, r) block
  BkCoefficients bk(a.n_states(), a.n_v());

  // Scattering grid: every (v, irotana) is independent.
#pragma omp parallel for schedule(static)
  for (int task = 0; task < a.n_v() * n_rot; ++task) {
    const int v = task / n_rot;
    const int irotana = task % n_rot;
    const int irot = a.rotation_index(irotana);
    const int j = a.j(irotana * a.n_diffraction());
    const auto phi = gas_.states.wave_function(v, j);

    // Projection onto phi_vj(r) at Zana, with the parallel momentum shifted out.
    std::vector<cplx> overlap(n_xy, cplx(0.0, 0.0));  // Fortran: pft2d(iy, ix)
    for (int ir = 0; ir < nr_; ++ir) {
      const double r_weight = phi[ir] * dr_;  // Fortran: rwave
      const auto block = psi.scattering.block(ir + nr_ * iz_analysis_);
      for (int ixy = 0; ixy < n_xy; ++ixy)
        overlap[ixy] = overlap[ixy] + tables_.shift_forward[ixy] * block[irot + n_rot_basis * ixy] * r_weight;
    }
    xy_to_momentum_->execute(overlap);

    for (int idif = 0; idif < a.n_diffraction(); ++idif) {
      const DiffractionChannel& c = a.diffraction_channels()[idif];
      bk(irotana * a.n_diffraction() + idif, v) = overlap[c.iy + ny_ * c.ix] * xy_scale_;
    }
  }

  // Specular grid: the initial channel at the fixed specular analysis surface.
  if (surface.fixed()) {
    const int iz = surface.index();
    for (int v = 0; v < a.n_v(); ++v) {
      const auto phi = gas_.states.wave_function(v, j_in_);
      cplx sum(0.0, 0.0);
      for (int ir = 0; ir < nr_; ++ir) sum = sum + psi.specular(ir, iz) * phi[ir];
      bk(a.initial_state(), v) = sum * dr_;
    }
  }
  return bk;
}

}  // namespace spodvr
