// Analysis surface on the specular grid. Fortran: Analyzor.f90 SetupAnalysis, FindNewAnalysisPoint
#include "propagation/analysis_surface.hpp"

#include "core/real_compare.hpp"

namespace spodvr {

namespace {

// sum over ir of |psi(ir, iz)|^2, evaluated like the Fortran conjg(psi)*psi sum.
double add_row_norm(double sum, const SpecularWave& psi, int iz) {
  for (int ir = 0; ir < psi.n_r(); ++ir) {
    const std::complex<double> c = psi(ir, iz);
    sum = sum + (c.real() * c.real() + c.imag() * c.imag());
  }
  return sum;
}

}  // namespace

// Fortran: Analyzor.f90 SetupAnalysis (nread = 0, nproj = 2)
AnalysisSurface::AnalysisSurface(const Config& config, const SpatialGrids& grids, const SpecularWave& initial)
    : z_(grids.z_specular),
      iz_analysis_(grids.iz_analysis),
      last_before_absorber_(grids.z_specular_optical->first - 1),
      norm_tolerance_(0.0),
      tolerance_(0.0),
      fixed_(config.derived.fixed_z_analysis_specular),
      index_(last_before_absorber_),
      previous_index_(last_before_absorber_) {
  const InputParameters& p = config.input;
  const double volume = p.dz * p.dr;

  if (fixed_) {
    // Fixed by the input: just before ZAnaSp if given, else at the start of the absorber.
    const bool at_z_ana_sp = p.z_analysis_specular && *p.z_analysis_specular > 0.0;
    const double z_fix = at_z_ana_sp ? *p.z_analysis_specular : p.zsp_opt_start;
    for (int iz = 0; iz < z_.size(); ++iz)
      if (fortran_compare::rge(z_[iz], z_fix)) {
        index_ = at_z_ana_sp ? iz - 1 : iz;
        break;
      }
    return;
  }

  // Dynamic: the tolerance is the initial norm inside the absorber, but at least tolsp.
  double norm = 0.0;
  for (int iz = grids.z_specular_optical->first; iz < z_.size(); ++iz) norm = add_row_norm(norm, initial, iz);
  norm = norm * volume;
  norm_tolerance_ = (norm > p.tol_specular ? norm : p.tol_specular) / volume;
  tolerance_ = norm_tolerance_ * volume;
}

// Fortran: Analyzor.f90 FindNewAnalysisPoint. Walks inwards from the absorber, accumulating
// the norm; the surface is put just beyond the first point where it exceeds the tolerance.
// Fortran quirks kept: the clamp to nschmi is overwritten on the next line, and reaching
// izanly fixes the surface without leaving the loop.
std::optional<double> AnalysisSurface::update(const SpecularWave& specular) {
  if (fixed_) return std::nullopt;
  std::optional<double> fixed_at;
  double norm = 0.0;
  for (int iz = last_before_absorber_; iz >= 0; --iz) {
    norm = add_row_norm(norm, specular, iz);
    if (norm > norm_tolerance_) {
      if (iz + 1 > previous_index_) {
        fixed_ = true;
        index_ = iz + 1;
        if (index_ >= last_before_absorber_) index_ = last_before_absorber_;
        fixed_at = z_[index_];
      }
      index_ = iz + 1;
      previous_index_ = index_;
      break;
    }
    if (iz == iz_analysis_) {
      fixed_ = true;
      index_ = iz_analysis_;
      fixed_at = z_[index_];
    }
  }
  return fixed_at;
}

}  // namespace spodvr
