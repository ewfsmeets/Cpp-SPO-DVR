// Kernels of the potential stage: absorption, the coupling to the specular grid, the angular
// transform FBR <-> DVR and exp(-i V dt). CPU twins: propagation/potential_operator.cpp
// (apply_block, apply_interaction) and propagation/angular_transform.cpp.
//
// The interaction region is processed in slabs of (Z, r) points; within a slab every (Z, r, x, y)
// point is one DVR "column" (the CPU transforms the columns of one (Z, r) point together).
// Every product and sum is formed in the same order as on the CPU, so for the same input these
// kernels give bitwise the same result as the CPU stage (the phi FFTs are cuFFT's).
#include "cuda/cuda_error.hpp"
#include "cuda/kernels.hpp"
#include "cuda/launch.cuh"

namespace spodvr::cuda {

namespace {

// Index of DVR value (theta, phi) of column col in the scratch array (col, theta, phi).
__device__ inline std::int64_t dvr_index(const PotentialArgs& a, std::int64_t col, int itheta, int iphi) {
  const std::int64_t n_columns = a.n_slab * a.psi.n_xy;
  return col + n_columns * (itheta + std::int64_t(a.n_theta) * iphi);
}

// One thread per element of the gas-phase region izr >= n_zr_potential; i runs izr fastest.
// CPU: PotentialOperator::apply_block, branch iz >= nz_potential.
__global__ void gas_phase_step_kernel(PotentialArgs a, Complex* psi) {
  const std::int64_t n_points = a.psi.n_zr - a.n_zr_potential;
  const std::int64_t i = thread_index();
  if (i >= n_points * a.psi.n_rot * a.psi.n_xy) return;
  const std::int64_t izr = a.n_zr_potential + i % n_points;
  const std::int64_t k = i / n_points;  // irot + n_rot * ixy
  const std::int64_t irot = k % a.psi.n_rot;
  const std::int64_t ir = izr % a.nr;
  const std::int64_t iz = izr / a.nr;
  Complex& v = psi[izr + a.psi.n_zr * k];
  v = a.z_optical[iz] * (a.gas_phase_step[irot + a.psi.n_rot * ir] * v);
}

// One thread per DVR column (izr, ixy) of the slab; neighbouring threads take neighbouring izr,
// so the psi reads (n_zr apart between rotational states) and the DVR writes are coalesced.
// For every mJ block: load its FBR coefficients with absorption and rot(dt/2) applied (and the
// specular wave added to the initial channel), then psi(theta_t, phi slot of mJ) =
// sum_J M(t, J) c_J. The matrix is read at the same address by all threads of a warp (one
// cached load, broadcast).
// CPU: PotentialOperator::apply_block (up to apply_interaction) and AngularTransform::fbr_to_dvr.
__global__ void absorb_and_transform_to_dvr_kernel(PotentialArgs a, const Complex* psi, const Complex* spec,
                                                   Complex* dvr) {
  const std::int64_t col = thread_index();
  if (col >= a.n_slab * a.psi.n_xy) return;
  const std::int64_t izr = a.izr_first + col % a.n_slab;
  const std::int64_t ixy = col / a.n_slab;
  const std::int64_t ir = izr % a.nr;
  const double optical_z = a.z_optical[izr / a.nr];
  const std::int64_t n_rot = a.psi.n_rot;

  // The specular wave of this (Z, r) point (index ir + nr*iz = izr), as it enters the initial channel.
  const Complex incoming =
      a.rotation_half[a.initial_channel + n_rot * ir] * (a.specular_to_scattering[ixy] * spec[izr]);

  for (int b = 0; b < a.n_mj_blocks; ++b) {
    const int first = a.block_first[b];
    const int n_j = a.block_size[b];
    Complex c[max_j_per_mj_block];
    for (int j = 0; j < n_j; ++j) {
      const int irot = first + j;
      c[j] = optical_z * (a.before_interaction[irot + n_rot * ir] * psi[a.psi(izr, irot, ixy)]);
      if (irot == a.initial_channel) c[j] = c[j] + incoming;
    }
    const double* m = a.matrix + a.block_matrix_offset[b];  // (theta, J) of this block
    for (int t = 0; t < a.n_theta; ++t) {
      Complex sum{0.0, 0.0};
      for (int j = 0; j < n_j; ++j) sum = sum + m[t * n_j + j] * c[j];
      dvr[dvr_index(a, col, t, a.block_phi_slot[b])] = sum;
    }
  }

  // The phi slots of |mJ| > MjMax carry no basis function.
  for (int t = 0; t < a.n_theta; ++t)
    for (int k = a.mj_max + 1; k <= a.n_phi - a.mj_max - 1; ++k) dvr[dvr_index(a, col, t, k)] = Complex{0.0, 0.0};
}

// One thread per DVR value (col, theta, phi); i runs col fastest, so the exp(i A) reads (izr
// fastest) are coalesced. Surface symmetry: the (x, y) point takes the value of its
// representative at the rotated / mirrored phi.
// CPU: PotentialOperator::apply_interaction (the loop between the transforms).
__global__ void multiply_by_exp_potential_kernel(PotentialArgs a, Complex* dvr) {
  const std::int64_t n_columns = a.n_slab * a.psi.n_xy;
  const std::int64_t i = thread_index();
  if (i >= n_columns * a.n_theta * a.n_phi) return;
  const std::int64_t col = i % n_columns;
  const int itheta = static_cast<int>((i / n_columns) % a.n_theta);
  const int iphi = static_cast<int>(i / (n_columns * a.n_theta));
  const std::int64_t izr = a.izr_first + col % a.n_slab;
  const std::int64_t ixy = col / a.n_slab;

  const std::int64_t n_tp = std::int64_t(a.n_theta) * a.n_phi;
  const int phi = a.phi_index[iphi + a.n_phi * a.operation_of[ixy]];
  const std::int64_t itp = phi + std::int64_t(a.n_phi) * itheta;
  const Complex b = a.exp_potential[izr + a.n_zr_potential * (itp + n_tp * a.representative_of[ixy])];
  dvr[i] = b * dvr[i];
}

// One thread per DVR column, as absorb_and_transform_to_dvr. For every FBR state J of a mJ
// block: c_J = sum_t M(t, J) psi(theta_t, phi slot), summed over theta in the CPU order; then
// rot(dt/2), and the result overwrites psi.
// CPU: AngularTransform::dvr_to_fbr and PotentialOperator::apply_rotation_half_step.
__global__ void transform_to_fbr_and_rotate_kernel(PotentialArgs a, const Complex* dvr, Complex* psi) {
  const std::int64_t col = thread_index();
  if (col >= a.n_slab * a.psi.n_xy) return;
  const std::int64_t izr = a.izr_first + col % a.n_slab;
  const std::int64_t ixy = col / a.n_slab;
  const std::int64_t ir = izr % a.nr;
  const std::int64_t n_rot = a.psi.n_rot;

  for (int b = 0; b < a.n_mj_blocks; ++b) {
    const int first = a.block_first[b];
    const int n_j = a.block_size[b];
    const double* m = a.matrix + a.block_matrix_offset[b];
    Complex c[max_j_per_mj_block];
    for (int j = 0; j < n_j; ++j) c[j] = Complex{0.0, 0.0};
    for (int t = 0; t < a.n_theta; ++t) {
      const Complex value = dvr[dvr_index(a, col, t, a.block_phi_slot[b])];
      for (int j = 0; j < n_j; ++j) c[j] = c[j] + m[t * n_j + j] * value;
    }
    for (int j = 0; j < n_j; ++j) {
      const int irot = first + j;
      psi[a.psi(izr, irot, ixy)] = a.rotation_half[irot + n_rot * ir] * c[j];
    }
  }
}

// One thread per (Z, r) point of the slab: the sum over (x, y) runs in the CPU order.
// CPU: PotentialOperator::apply_block, "Project the initial channel back onto the specular grid".
__global__ void project_onto_specular_kernel(PotentialArgs a, Complex* psi, Complex* spec) {
  const std::int64_t point = thread_index();
  if (point >= a.n_slab) return;
  const std::int64_t izr = a.izr_first + point;
  Complex projection{0.0, 0.0};
  for (std::int64_t ixy = 0; ixy < a.psi.n_xy; ++ixy)
    projection = projection + a.scattering_to_specular[ixy] * psi[a.psi(izr, a.initial_channel, ixy)];
  spec[izr] = projection;
  for (std::int64_t ixy = 0; ixy < a.psi.n_xy; ++ixy) {
    Complex& v = psi[a.psi(izr, a.initial_channel, ixy)];
    v = v - a.specular_to_scattering[ixy] * projection;
  }
}

}  // namespace

// Fortran: Absorption.f90 PerformEXPOpticalZ/R, GasPhase.f90 PerformEXPGasPhaseOnScatGrid
void apply_gas_phase_step(const PotentialArgs& args, Complex* psi) {
  const std::int64_t n = (args.psi.n_zr - args.n_zr_potential) * args.psi.n_rot * args.psi.n_xy;
  if (n == 0) return;
  gas_phase_step_kernel<<<blocks_for(n), threads_per_block>>>(args, psi);
  check_kernel_launch("apply_gas_phase_step");
}

// Fortran: Absorption.f90 PerformEXPOpticalZ/R, TwoGrids.f90 AddSpecToScat,
// transform.f90 FBRToDVRGenComplexTransform (Legendre part)
void absorb_and_transform_to_dvr(const PotentialArgs& args, const Complex* psi, const Complex* specular,
                                 Complex* dvr) {
  const std::int64_t n = args.n_slab * args.psi.n_xy;
  absorb_and_transform_to_dvr_kernel<<<blocks_for(n), threads_per_block>>>(args, psi, specular, dvr);
  check_kernel_launch("absorb_and_transform_to_dvr");
}

// Fortran: Potential.f90 PerformEXPInteractionOp
void multiply_by_exp_potential(const PotentialArgs& args, Complex* dvr) {
  const std::int64_t n = args.n_slab * args.psi.n_xy * args.n_theta * args.n_phi;
  multiply_by_exp_potential_kernel<<<blocks_for(n), threads_per_block>>>(args, dvr);
  check_kernel_launch("multiply_by_exp_potential");
}

// Fortran: transform.f90 DVRToFBRGenComplexTransform (Legendre part), EXPKinEn%halfStep%RotChan
void transform_to_fbr_and_rotate(const PotentialArgs& args, const Complex* dvr, Complex* psi) {
  const std::int64_t n = args.n_slab * args.psi.n_xy;
  transform_to_fbr_and_rotate_kernel<<<blocks_for(n), threads_per_block>>>(args, dvr, psi);
  check_kernel_launch("transform_to_fbr_and_rotate");
}

// Fortran: TwoGrids.f90 SubSpecFromScat
void project_onto_specular(const PotentialArgs& args, Complex* psi, Complex* specular) {
  project_onto_specular_kernel<<<blocks_for(args.n_slab), threads_per_block>>>(args, psi, specular);
  check_kernel_launch("project_onto_specular");
}

}  // namespace spodvr::cuda
