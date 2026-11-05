// Potential energy propagator on the scattering grid, incl. absorption and the coupling to the
// specular grid. Fortran: Potential.f90, Absorption.f90, TwoGrids.f90
#include "propagation/potential_operator.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace spodvr {

namespace {
using cplx = std::complex<double>;

// exp(i A) as the Fortran code forms it: B = COS(A) + ima*SIN(A).
cplx exp_i(double a) { return cplx(std::cos(a), std::sin(a)); }
}  // namespace

PotentialOperator::PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                                     const PropagationTables& tables, const AngularTransform& transform,
                                     std::span<const double> premultiplied_potential, std::size_t table_limit)
    : PotentialOperator(shape, nz_potential, initial_channel, tables, transform, premultiplied_potential,
                        GridSymmetry(shape.nx, shape.ny, transform.n_phi()), table_limit) {}

PotentialOperator::PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                                     const PropagationTables& tables, const AngularTransform& transform,
                                     std::span<const double> premultiplied_potential, GridSymmetry symmetry,
                                     std::size_t table_limit)
    : shape_(shape),
      nz_potential_(nz_potential),
      initial_channel_(initial_channel),
      tables_(tables),
      transform_(transform),
      potential_(premultiplied_potential),
      symmetry_(std::move(symmetry)) {
  if (symmetry_.nx() != shape.nx || symmetry_.ny() != shape.ny || symmetry_.nphi() != transform.n_phi())
    throw std::invalid_argument("PotentialOperator: the symmetry does not match the grid");
  const std::size_t expected = std::size_t(transform.n_theta()) * transform.n_phi() * symmetry_.n_representatives() *
                               std::size_t(nz_potential) * shape.nr;
  if (potential_.size() != expected) throw std::invalid_argument("PotentialOperator: potential has the wrong size");

  if ((shape.nx * shape.ny) % transform.n_columns() != 0)
    throw std::invalid_argument("PotentialOperator: nx*ny must be a multiple of the transform's columns");

  // Merged per-(irot, ir) factors (one multiplication instead of several):
  //   interaction region: exp(-V_opt,r dt) rot(dt/2)  (before the DVR step)
  //   beyond it:          exp(-V_opt,r dt) rot(dt/2) exp(-i vasy dt) rot(dt/2)
  const int n_rot = shape.n_rot;
  for (int ir = 0; ir < shape.nr; ++ir)
    for (int irot = 0; irot < n_rot; ++irot) {
      const cplx rot = tables.rotation.half[irot + n_rot * ir];
      const double optical_r = shape.nr != 1 ? tables.r_optical[ir] : 1.0;
      before_interaction_.push_back(optical_r * rot);
      gas_phase_step_.push_back(optical_r * rot * tables.gas_phase[ir] * rot);
    }

  if (expected * sizeof(cplx) <= table_limit) {
    exp_potential_.resize(expected);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < expected; ++i) exp_potential_[i] = exp_i(potential_[i]);
    potential_ = {};
  }
}

std::size_t PotentialOperator::dvr_work_size() const {
  return std::size_t(transform_.dvr_block_size()) * shape_.nx * shape_.ny;
}

// The DVR blocks of all (x, y) points, plus (on-the-fly path with symmetry) exp(i A) of the
// representatives of one (Z, r) point.
std::size_t PotentialOperator::work_size() const {
  const bool scratch = !uses_precomputed_exponentials() && reduces();
  const std::size_t n_tp = std::size_t(transform_.n_theta()) * transform_.n_phi();
  return dvr_work_size() + (scratch ? n_tp * symmetry_.n_representatives() : 0);
}

void PotentialOperator::apply(WavePacket& psi) const {
#pragma omp parallel
  {
    ComplexVector work(work_size());
#pragma omp for schedule(static)
    for (int izr = 0; izr < psi.scattering.n_blocks(); ++izr)
      apply_block(izr, psi.scattering.block(izr), psi.specular, work);
  }
}

// Fortran: Potential.f90 PerformEXPPotentialOp (called per izr in KineticAndPotentialOp), with
// the scalar factors merged: optical Z (real) times the per-(irot, ir) tables of the constructor.
void PotentialOperator::apply_block(int izr, std::span<cplx> block, SpecularWave& spec, std::span<cplx> work) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  const int ichn = initial_channel_;
  const int ir = izr % shape_.nr;
  const int iz = izr / shape_.nr;
  const double optical_z = tables_.z_optical[iz];
  const auto rot_half = std::span(tables_.rotation.half).subspan(std::size_t(n_rot) * ir, n_rot);

  if (iz >= nz_potential_) {
    // Beyond the interaction region: absorption and the free molecule.
    // Fortran: PerformEXPOpticalZ/R, PerformEXPGasPhaseOnScatGrid between rot(dt/2)
    const auto factor = std::span(gas_phase_step_).subspan(std::size_t(n_rot) * ir, n_rot);
    for (int ixy = 0; ixy < n_xy; ++ixy)
      for (int irot = 0; irot < n_rot; ++irot)
        block[irot + n_rot * ixy] = optical_z * (factor[irot] * block[irot + n_rot * ixy]);
    return;
  }

  // Absorption, rot(dt/2), and the specular wave added to the initial channel (it enters after
  // the absorption, before rot(dt/2)). Fortran: PerformEXPOpticalZ/R, AddSpecToScat
  const auto factor = std::span(before_interaction_).subspan(std::size_t(n_rot) * ir, n_rot);
  const cplx incoming = spec(ir, iz);
  for (int ixy = 0; ixy < n_xy; ++ixy) {
    for (int irot = 0; irot < n_rot; ++irot)
      block[irot + n_rot * ixy] = optical_z * (factor[irot] * block[irot + n_rot * ixy]);
    block[ichn + n_rot * ixy] =
        block[ichn + n_rot * ixy] + rot_half[ichn] * (tables_.specular_to_scattering[ixy] * incoming);
  }

  // Fortran: PerformEXPInteractionOp
  apply_interaction(block, izr, work);
  apply_rotation_half_step(block, ir);

  // Project the initial channel back onto the specular grid. Fortran: SubSpecFromScat
  cplx projection(0.0, 0.0);
  for (int ixy = 0; ixy < n_xy; ++ixy)
    projection = projection + tables_.scattering_to_specular[ixy] * block[ichn + n_rot * ixy];
  spec(ir, iz) = projection;
  for (int ixy = 0; ixy < n_xy; ++ixy)
    block[ichn + n_rot * ixy] = block[ichn + n_rot * ixy] - tables_.specular_to_scattering[ixy] * projection;
}

// exp(-i E_rot dt/2) for every rotational state. Fortran: EXPKinEn%halfStep%RotChan
void PotentialOperator::apply_rotation_half_step(std::span<cplx> block, int ir) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot)
      block[irot + n_rot * ixy] = tables_.rotation.half[irot + n_rot * ir] * block[irot + n_rot * ixy];
}

// exp(i A) for every (theta, phi) point of every representative of izr (layout (itp, irep)).
void PotentialOperator::exponentials_of_block(int izr, std::span<cplx> scratch) const {
  const std::size_t block_size = scratch.size();
  const std::size_t first = block_size * izr;
  for (std::size_t i = 0; i < block_size; ++i) scratch[i] = exp_i(potential_[first + i]);
}

// FBR -> DVR, multiply by exp(i A) with A = -V dt, DVR -> FBR, for every (x, y) point. A of
// point (ixy, theta, phi) is that of its representative at the mapped phi (surface symmetry;
// without symmetry the map is the identity). Fortran: PerformEXPInteractionOp
// (B = COS(A) + ima*SIN(A); with symmetry A = pV(symmetryIndex(index1)))
void PotentialOperator::apply_interaction(std::span<cplx> block, int izr, std::span<cplx> work) const {
  const int n_rot = shape_.n_rot;
  const int n_xy = shape_.nx * shape_.ny;
  const int n_theta = transform_.n_theta();
  const int n_phi = transform_.n_phi();
  const int ld = transform_.leading_dimension();
  const int dvr_size = transform_.dvr_block_size();
  const std::size_t n_tp = std::size_t(n_theta) * n_phi;
  const std::size_t n_rep = symmetry_.n_representatives();
  const std::span<cplx> dvr = work.first(dvr_work_size());

  const int n_col = transform_.n_columns();
  for (int c0 = 0; c0 < n_xy; c0 += n_col)
    transform_.fbr_to_dvr(block.subspan(n_rot * c0, n_rot * n_col), dvr.subspan(dvr_size * c0, dvr_size * n_col));

  // Where exp(i A) of this (Z, r) point comes from: the precomputed table, the scratch array
  // (on the fly, with symmetry: each representative's value is used by several (x, y) points),
  // or (on the fly, without symmetry) evaluated at the point itself.
  const bool precomputed = uses_precomputed_exponentials();
  const bool scratch = !precomputed && reduces();
  std::span<const cplx> exponentials;
  if (precomputed) exponentials = std::span<const cplx>(exp_potential_).subspan(n_tp * n_rep * izr, n_tp * n_rep);
  if (scratch) {
    const std::span<cplx> s = work.subspan(dvr_work_size(), n_tp * n_rep);
    exponentials_of_block(izr, s);
    exponentials = s;
  }
  const std::span<const double> potential =
      precomputed ? std::span<const double>() : potential_.subspan(n_tp * n_rep * izr, n_tp * n_rep);

  for (int ixy = 0; ixy < n_xy; ++ixy) {
    const std::size_t v0 = n_tp * symmetry_.representative_of(ixy);
    const std::vector<int>& phi_index = symmetry_.phi_index(symmetry_.operation_of(ixy));
    for (int itheta = 0; itheta < n_theta; ++itheta)
      for (int iphi = 0; iphi < n_phi; ++iphi) {
        const std::size_t i = v0 + phi_index[iphi] + n_phi * itheta;
        const cplx b = (precomputed || scratch) ? exponentials[i] : exp_i(potential[i]);
        cplx& v = dvr[dvr_size * ixy + iphi + ld * itheta];
        v = b * v;
      }
  }

  for (int c0 = 0; c0 < n_xy; c0 += n_col)
    transform_.dvr_to_fbr(dvr.subspan(dvr_size * c0, dvr_size * n_col), block.subspan(n_rot * c0, n_rot * n_col));
}

}  // namespace spodvr
