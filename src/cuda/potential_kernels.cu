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

// --- index helpers ----------------------------------------------------------------------------

// The (Z, r, x, y) point of DVR column col of the current slab: col = (izr - izr_first) + n_slab * ixy.
struct DvrColumn {
  std::int64_t izr;
  std::int64_t ixy;
};
__device__ inline DvrColumn dvr_column(const PotentialArgs& args, std::int64_t col) {
  return {args.izr_first + col % args.n_slab, col / args.n_slab};
}

// DVR columns of a slab (one per (Z, r, x, y) point).
__device__ inline std::int64_t n_dvr_columns(const PotentialArgs& args) { return args.n_slab * args.psi.n_xy; }

// Index of DVR value (itheta, iphi) of column col in the scratch array (col, theta, phi).
__device__ inline std::int64_t dvr_index(const PotentialArgs& args, std::int64_t col, int itheta, int iphi) {
  return col + n_dvr_columns(args) * (itheta + std::int64_t(args.n_theta) * iphi);
}

// Index into the tables per rotational state and r point, (irot, ir): before_interaction,
// gas_phase_step, rotation_half.
__device__ inline std::int64_t irot_ir_index(const PotentialArgs& args, std::int64_t irot, std::int64_t ir) {
  return irot + args.psi.n_rot * ir;
}

// --- kernels ----------------------------------------------------------------------------------

// The gas-phase region izr >= n_zr_potential: x over izr (coalesced), y over
// rot_xy = irot + n_rot * ixy (grid_2d). CPU: PotentialOperator::apply_block, branch iz >= nz_potential.
__global__ void gas_phase_step_kernel(PotentialArgs args, Complex* psi) {
  const std::int64_t izr = args.n_zr_potential + thread_index();
  if (izr >= args.psi.n_zr) return;
  const int ir = static_cast<int>(izr % args.nr);
  const double optical_z = args.z_optical[izr / args.nr];
  const int n_rot = static_cast<int>(args.psi.n_rot);
  for (std::int64_t rot_xy = blockIdx.y; rot_xy < args.psi.n_rot * args.psi.n_xy; rot_xy += gridDim.y) {
    const int irot = static_cast<int>(rot_xy % n_rot);
    Complex& v = psi[izr + args.psi.n_zr * rot_xy];
    v = optical_z * (args.gas_phase_step[irot_ir_index(args, irot, ir)] * v);
  }
}

// One thread per DVR column (izr, ixy) of the slab; neighbouring threads take neighbouring izr,
// so the psi reads (n_zr apart between rotational states) and the DVR writes are coalesced.
// For every mJ block:
//   fbr[j]           = the block's FBR coefficients with absorption and rot(dt/2) applied (and
//                      the specular wave added to the initial channel),
//   dvr(itheta, mJ)  = sum_j legendre(itheta, j) fbr[j], at the phi slot of mJ.
// The Legendre matrix is read at the same address by all threads of a warp (one cached load,
// broadcast). CPU: PotentialOperator::apply_block (up to apply_interaction) and
// AngularTransform::fbr_to_dvr.
__global__ void absorb_and_transform_to_dvr_kernel(PotentialArgs args, const Complex* psi, const Complex* spec,
                                                   Complex* dvr) {
  const std::int64_t col = thread_index();
  if (col >= n_dvr_columns(args)) return;
  const auto [izr, ixy] = dvr_column(args, col);
  const std::int64_t ir = izr % args.nr;
  const double optical_z = args.z_optical[izr / args.nr];

  // The specular wave of this (Z, r) point (index ir + nr*iz = izr), as it enters the initial channel.
  const Complex incoming = args.rotation_half[irot_ir_index(args, args.initial_channel, ir)] *
                           (args.specular_to_scattering[ixy] * spec[izr]);

  for (int block = 0; block < args.n_mj_blocks; ++block) {
    const int first_irot = args.block_first[block];
    const int n_j = args.block_size[block];
    Complex fbr[max_j_per_mj_block];
    for (int j = 0; j < n_j; ++j) {
      const int irot = first_irot + j;
      fbr[j] = optical_z * (args.before_interaction[irot_ir_index(args, irot, ir)] * psi[args.psi(izr, irot, ixy)]);
      if (irot == args.initial_channel) fbr[j] = fbr[j] + incoming;
    }
    const double* legendre = args.matrix + args.block_matrix_offset[block];  // (theta, j) of this block
    for (int itheta = 0; itheta < args.n_theta; ++itheta) {
      Complex sum{0.0, 0.0};
      for (int j = 0; j < n_j; ++j) sum = sum + legendre[itheta * n_j + j] * fbr[j];
      dvr[dvr_index(args, col, itheta, args.block_phi_slot[block])] = sum;
    }
  }

  // The phi slots of |mJ| > MjMax carry no basis function.
  for (int itheta = 0; itheta < args.n_theta; ++itheta)
    for (int iphi = args.mj_max + 1; iphi <= args.n_phi - args.mj_max - 1; ++iphi)
      dvr[dvr_index(args, col, itheta, iphi)] = Complex{0.0, 0.0};
}

// x over the DVR columns of the slab, y over the (theta, phi) points tp = itheta + n_theta * iphi
// (grid_2d); neighbouring threads take neighbouring izr, so the exp(i A) reads (izr fastest) and
// the DVR accesses are coalesced. Surface symmetry: see ExpPotentialTable.
// CPU: PotentialOperator::apply_interaction (the loop between the transforms).
__global__ void multiply_by_exp_potential_kernel(PotentialArgs args, Complex* dvr) {
  const std::int64_t col = thread_index();
  if (col >= n_dvr_columns(args)) return;
  const auto [izr, ixy] = dvr_column(args, col);
  for (int tp = blockIdx.y; tp < args.n_theta * args.n_phi; tp += gridDim.y) {
    const int itheta = tp % args.n_theta;
    const int iphi = tp / args.n_theta;
    Complex& v = dvr[dvr_index(args, col, itheta, iphi)];
    v = args.exp_potential.at(izr, ixy, itheta, iphi) * v;
  }
}

// One thread per DVR column, as absorb_and_transform_to_dvr_kernel. For every mJ block:
//   fbr[j] = sum_itheta legendre(itheta, j) dvr(itheta, mJ), summed over theta in the CPU order,
// then rot(dt/2), and the result overwrites psi.
// CPU: AngularTransform::dvr_to_fbr and PotentialOperator::apply_rotation_half_step.
__global__ void transform_to_fbr_and_rotate_kernel(PotentialArgs args, const Complex* dvr, Complex* psi) {
  const std::int64_t col = thread_index();
  if (col >= n_dvr_columns(args)) return;
  const auto [izr, ixy] = dvr_column(args, col);
  const std::int64_t ir = izr % args.nr;

  for (int block = 0; block < args.n_mj_blocks; ++block) {
    const int first_irot = args.block_first[block];
    const int n_j = args.block_size[block];
    const double* legendre = args.matrix + args.block_matrix_offset[block];  // (theta, j) of this block
    Complex fbr[max_j_per_mj_block];
    for (int j = 0; j < n_j; ++j) fbr[j] = Complex{0.0, 0.0};
    for (int itheta = 0; itheta < args.n_theta; ++itheta) {
      const Complex value = dvr[dvr_index(args, col, itheta, args.block_phi_slot[block])];
      for (int j = 0; j < n_j; ++j) fbr[j] = fbr[j] + legendre[itheta * n_j + j] * value;
    }
    for (int j = 0; j < n_j; ++j) {
      const int irot = first_irot + j;
      psi[args.psi(izr, irot, ixy)] = args.rotation_half[irot_ir_index(args, irot, ir)] * fbr[j];
    }
  }
}

// One thread per (Z, r) point of the slab: the sum over (x, y) runs in the CPU order.
// CPU: PotentialOperator::apply_block, "Project the initial channel back onto the specular grid".
__global__ void project_onto_specular_kernel(PotentialArgs args, Complex* psi, Complex* spec) {
  const std::int64_t point = thread_index();
  if (point >= args.n_slab) return;
  const std::int64_t izr = args.izr_first + point;
  Complex projection{0.0, 0.0};
  for (std::int64_t ixy = 0; ixy < args.psi.n_xy; ++ixy)
    projection = projection + args.scattering_to_specular[ixy] * psi[args.psi(izr, args.initial_channel, ixy)];
  spec[izr] = projection;
  for (std::int64_t ixy = 0; ixy < args.psi.n_xy; ++ixy) {
    Complex& v = psi[args.psi(izr, args.initial_channel, ixy)];
    v = v - args.specular_to_scattering[ixy] * projection;
  }
}

}  // namespace

// Fortran: Absorption.f90 PerformEXPOpticalZ/R, GasPhase.f90 PerformEXPGasPhaseOnScatGrid
void apply_gas_phase_step(const PotentialArgs& args, Complex* psi) {
  const std::int64_t n_points = args.psi.n_zr - args.n_zr_potential;
  if (n_points == 0) return;
  gas_phase_step_kernel<<<grid_2d(n_points, args.psi.n_rot * args.psi.n_xy), threads_per_block>>>(args, psi);
  check_kernel_launch("apply_gas_phase_step");
}

// Fortran: Absorption.f90 PerformEXPOpticalZ/R, TwoGrids.f90 AddSpecToScat,
// transform.f90 FBRToDVRGenComplexTransform (Legendre part)
void absorb_and_transform_to_dvr(const PotentialArgs& args, const Complex* psi, const Complex* specular,
                                 Complex* dvr) {
  const std::int64_t n_columns = args.n_slab * args.psi.n_xy;
  absorb_and_transform_to_dvr_kernel<<<blocks_for(n_columns), threads_per_block>>>(args, psi, specular, dvr);
  check_kernel_launch("absorb_and_transform_to_dvr");
}

// Fortran: Potential.f90 PerformEXPInteractionOp
void multiply_by_exp_potential(const PotentialArgs& args, Complex* dvr) {
  const std::int64_t n_columns = args.n_slab * args.psi.n_xy;
  multiply_by_exp_potential_kernel<<<grid_2d(n_columns, std::int64_t(args.n_theta) * args.n_phi),
                                     threads_per_block>>>(args, dvr);
  check_kernel_launch("multiply_by_exp_potential");
}

// Fortran: transform.f90 DVRToFBRGenComplexTransform (Legendre part), EXPKinEn%halfStep%RotChan
void transform_to_fbr_and_rotate(const PotentialArgs& args, const Complex* dvr, Complex* psi) {
  const std::int64_t n_columns = args.n_slab * args.psi.n_xy;
  transform_to_fbr_and_rotate_kernel<<<blocks_for(n_columns), threads_per_block>>>(args, dvr, psi);
  check_kernel_launch("transform_to_fbr_and_rotate");
}

// Fortran: TwoGrids.f90 SubSpecFromScat
void project_onto_specular(const PotentialArgs& args, Complex* psi, Complex* specular) {
  project_onto_specular_kernel<<<blocks_for(args.n_slab), threads_per_block>>>(args, psi, specular);
  check_kernel_launch("project_onto_specular");
}

}  // namespace spodvr::cuda
