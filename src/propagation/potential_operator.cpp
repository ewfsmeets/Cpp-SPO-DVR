// Potential energy propagator on the scattering grid, incl. absorption and the coupling to the
// specular grid. Fortran: Potential.f90, Absorption.f90, TwoGrids.f90
#include "propagation/potential_operator.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace spodvr {

namespace {
using cplx = std::complex<double>;
}  // namespace

PotentialOperator::PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                                     const PropagationTables& tables, const AngularTransform& transform,
                                     std::span<const double> premultiplied_potential)
    : shape_(shape),
      nz_potential_(nz_potential),
      initial_channel_(initial_channel),
      tables_(tables),
      transform_(transform),
      potential_(premultiplied_potential) {
  const std::size_t expected = std::size_t(transform.n_theta()) * transform.n_phi() * shape.nx * shape.ny *
                               std::size_t(nz_potential) * shape.nr;
  if (potential_.size() != expected) throw std::invalid_argument("PotentialOperator: potential has the wrong size");
}

// Fortran: Potential.f90 PerformEXPPotentialOp (called per izr in KineticAndPotentialOp)
void PotentialOperator::apply(WavePacket& psi) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  const int ichn = initial_channel_;
  ScatteringWave& scat = psi.scattering;
  SpecularWave& spec = psi.specular;

#pragma omp parallel
  {
    std::vector<cplx> dvr(std::size_t(transform_.dvr_block_size()) * n_xy);  // work space per thread
#pragma omp for schedule(static)
    for (int izr = 0; izr < scat.n_blocks(); ++izr) {
      const int ir = izr % shape_.nr;
      const int iz = izr / shape_.nr;
      const std::span<cplx> block = scat.block(izr);
      const bool interaction_region = iz < nz_potential_;

      // Absorbing potentials. Fortran: PerformEXPOpticalZ, PerformEXPOpticalR
      for (cplx& v : block) v = tables_.z_optical[iz] * v;
      if (shape_.nr != 1)
        for (cplx& v : block) v = tables_.r_optical[ir] * v;

      // Add the specular wave to the initial channel. Fortran: AddSpecToScat
      if (interaction_region)
        for (int ixy = 0; ixy < n_xy; ++ixy)
          block[ichn + n_rot * ixy] = block[ichn + n_rot * ixy] + tables_.specular_to_scattering[ixy] * spec(ir, iz);

      // Fortran: PerformEXPInteractionOp
      apply_rotation_half_step(block, ir);
      if (interaction_region) {
        apply_interaction(block, izr, dvr);
      } else {
        for (cplx& v : block) v = v * tables_.gas_phase[ir];  // Fortran: PerformEXPGasPhaseOnScatGrid
      }
      apply_rotation_half_step(block, ir);

      // Project the initial channel back onto the specular grid. Fortran: SubSpecFromScat
      if (interaction_region) {
        cplx projection(0.0, 0.0);
        for (int ixy = 0; ixy < n_xy; ++ixy)
          projection = projection + tables_.scattering_to_specular[ixy] * block[ichn + n_rot * ixy];
        spec(ir, iz) = projection;
        for (int ixy = 0; ixy < n_xy; ++ixy)
          block[ichn + n_rot * ixy] = block[ichn + n_rot * ixy] - tables_.specular_to_scattering[ixy] * projection;
      }
    }
  }
}

// exp(-i E_rot dt/2) for every rotational state. Fortran: EXPKinEn%halfStep%RotChan
void PotentialOperator::apply_rotation_half_step(std::span<cplx> block, int ir) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot)
      block[irot + n_rot * ixy] = tables_.rotation.half[irot + n_rot * ir] * block[irot + n_rot * ixy];
}

// FBR -> DVR, multiply by exp(i A) with A = -V dt, DVR -> FBR, for every (x, y) point.
// Fortran: PerformEXPInteractionOp (B = COS(A) + ima*SIN(A))
void PotentialOperator::apply_interaction(std::span<cplx> block, int izr, std::span<cplx> dvr) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  const int n_theta = transform_.n_theta();
  const int n_phi = transform_.n_phi();
  const int ld = transform_.leading_dimension();
  const int dvr_size = transform_.dvr_block_size();
  const std::size_t n_tp = std::size_t(n_theta) * n_phi;

  for (int ixy = 0; ixy < n_xy; ++ixy)
    transform_.fbr_to_dvr(block.subspan(n_rot * ixy, n_rot), dvr.subspan(dvr_size * ixy, dvr_size));

  for (int ixy = 0; ixy < n_xy; ++ixy) {
    const std::size_t v0 = n_tp * (ixy + std::size_t(n_xy) * izr);
    for (int itheta = 0; itheta < n_theta; ++itheta)
      for (int iphi = 0; iphi < n_phi; ++iphi) {
        const double a = potential_[v0 + iphi + n_phi * itheta];
        const cplx b(std::cos(a), std::sin(a));
        cplx& v = dvr[dvr_size * ixy + iphi + ld * itheta];
        v = b * v;
      }
  }

  for (int ixy = 0; ixy < n_xy; ++ixy)
    transform_.dvr_to_fbr(dvr.subspan(dvr_size * ixy, dvr_size), block.subspan(n_rot * ixy, n_rot));
}

}  // namespace spodvr
