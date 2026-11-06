#pragma once

// The CUDA kernels of the GPU backend, as plain launch functions (the .cu files define them).
// Every launcher starts its kernel(s) on the default stream and checks the launch; they run one
// after the other in launch order, like the CPU code. Thread -> data mapping and the CPU twin are
// documented at each definition.
//
// Pointers are device pointers (src/cuda/ exception to the pointer rule).

#include <cstdint>

#include "cuda/device_complex.hpp"

namespace spodvr::cuda {

// --- kinetic stages (kinetic_kernels.cu) ---------------------------------------------------
// The FFT scale factors of a transform pair (cuFFT does not scale). They are applied as two
// separate multiplications, like the CPU plans do: folding their product into the factor table
// would round it once more, a bias that repeats in every time step (it made the norm drift by
// 2.6e-13 relative from the CPU's in 1000 steps).
struct FftScales {
  double forward = 1.0;
  double backward = 1.0;
};
// data[i] = factors[i % period] * (backward * (forward * data[i])): the specular kinetic step
// (period = n).
void multiply_by_periodic_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t period,
                                  std::int64_t n);
// data[i] = factors[i / block] * (backward * (forward * data[i])): the momentum shift on (x, y)
// (block = n_zr * n_rot, the values of one ixy; scales 1 for the shift).
void multiply_by_blocked_factors(Complex* data, const Complex* factors, FftScales scales, std::int64_t block,
                                 std::int64_t n);
// The merged kinetic step on psi in (k_Z, k_r, k_x, k_y) space (layout (izr, irot, ixy)): first the
// (Z, r) factor, then the (x, y) factor, each with its own scales, in the CPU's stage order:
//   v = zr[izr] * (zr_scales.backward * (zr_scales.forward * v))
//   v = xy[ixy] * (xy_scales.backward * (xy_scales.forward * v))
// (docs/cuda.md, methodology exception "merged kinetic stages").
void multiply_by_kinetic_factors(Complex* psi, const DeviceLayout& layout, const Complex* zr_factors,
                                 FftScales zr_scales, const Complex* xy_factors, FftScales xy_scales);

// --- specular grid (specular_kernels.cu) ---------------------------------------------------
struct SpecularStepArgs {
  int nr;
  int nz;               // specular Z points
  int nz_potential;     // Z < zasymp
  const Complex* gas_phase_step;  // per ir: rot(dt/2) exp(-i vasy dt) rot(dt/2)
  int z_absorber_first, z_absorber_last;
  const double* z_optical;        // per iz
  int r_absorber_first, r_absorber_last;
  const double* r_optical;        // per ir (not applied for nr = 1)
};
// The gas-phase propagator beyond the interaction region, then the Z and r absorbers.
void apply_specular_gas_phase_and_absorbers(const SpecularStepArgs& args, Complex* specular);

// --- potential stage (potential_kernels.cu) ------------------------------------------------
// exp(i A), A = -V dt, at the (x, y) representatives of the surface symmetry, as the potential
// kernels and the cuFFT callback (fft_callbacks.cu) read it. Values (izr, itp, irep), izr fastest,
// with itp = iphi + n_phi * itheta. Point ixy takes the values of its representative
// representative_of[ixy] at the mapped phi phi_index[iphi + n_phi * operation_of[ixy]] (rotated
// or mirrored; the identity without symmetry).
struct ExpPotentialTable {
  const Complex* values;
  const int* representative_of;  // per ixy
  const int* operation_of;       // per ixy
  const int* phi_index;          // per symmetry operation: n_phi values
  std::int64_t n_zr_potential;
  int n_theta;
  int n_phi;

  // exp(i A) at DVR point (itheta, iphi) of the (Z, r) point izr and the (x, y) point ixy.
  SPODVR_HOST_DEVICE Complex at(std::int64_t izr, std::int64_t ixy, int itheta, int iphi) const {
    const std::int64_t mapped_phi = phi_index[n_phi * operation_of[ixy] + iphi];
    const std::int64_t itp = mapped_phi + std::int64_t(n_phi) * itheta;
    const std::int64_t n_tp = std::int64_t(n_theta) * n_phi;
    return values[izr + n_zr_potential * (n_tp * representative_of[ixy] + itp)];
  }
};

// Everything the potential kernels read. One (Z, r) slab of the interaction region at a time:
// points izr_first .. izr_first + n_slab - 1, i.e. n_slab * n_xy DVR columns, column
// col = (izr - izr_first) + n_slab * ixy. DVR scratch layout (col, theta, phi), phi slowest.
struct PotentialArgs {
  DeviceLayout psi;            // layout of the scattering wave function
  int nr;
  int n_theta;
  int n_phi;
  int initial_channel;         // Fortran: ichn - 1
  std::int64_t n_zr_potential; // (Z, r) points with Z < zasymp: izr < n_zr_potential

  std::int64_t izr_first;      // the slab
  std::int64_t n_slab;

  const double* z_optical;                 // exp(-V_opt,Z dt) per iz
  const Complex* before_interaction;       // (irot, ir), PotentialStepFactors
  const Complex* gas_phase_step;           // (irot, ir), PotentialStepFactors
  const Complex* rotation_half;            // (irot, ir): exp(-i E_rot dt/2)
  const Complex* specular_to_scattering;   // per ixy. Fortran: initXY
  const Complex* scattering_to_specular;   // per ixy. Fortran: XYProj

  // FBR <-> DVR (AngularTransform tables): mJ block b has the FBR states first[b] ..
  // first[b] + size[b] - 1, phi slot phi_slot[b] and its (theta, J) matrix at matrix_offset[b].
  int n_mj_blocks;
  int mj_max;
  const int* block_first;
  const int* block_size;
  const int* block_phi_slot;
  const int* block_matrix_offset;
  const double* matrix;

  ExpPotentialTable exp_potential;
};
// Largest number of J values in one mJ block that the DVR kernels support.
inline constexpr int max_j_per_mj_block = 64;

// Beyond the interaction region (izr >= n_zr_potential): absorption and the free molecule.
void apply_gas_phase_step(const PotentialArgs& args, Complex* psi);
// Interaction region, one slab: absorption, rot(dt/2), add the specular wave to the initial
// channel, FBR -> DVR in theta (the phi FFT follows separately).
void absorb_and_transform_to_dvr(const PotentialArgs& args, const Complex* psi, const Complex* specular, Complex* dvr);
// dvr *= exp(i A) at every (theta, phi) point of the slab. The propagation applies exp(i A) in the
// load callback of the inverse phi FFT instead; this kernel serves for nphi = 1 (no phi FFT) and
// as the reference of the callback (test).
void multiply_by_exp_potential(const PotentialArgs& args, Complex* dvr);
// DVR -> FBR in theta (after the phi FFT), then rot(dt/2); writes the slab of psi.
void transform_to_fbr_and_rotate(const PotentialArgs& args, const Complex* dvr, Complex* psi);
// Projects the initial channel onto the specular grid and subtracts it from the scattering wave.
void project_onto_specular(const PotentialArgs& args, Complex* psi, Complex* specular);

// --- whole wave function (wave_function_kernels.cu) ----------------------------------------
// densities[izr] = sum over (irot, ixy) of |psi|^2, summed in the host memory order.
void block_densities(const Complex* psi, const DeviceLayout& layout, double* densities);
// The nr (Z, r) blocks of Z index iz in host layout (irot, ixy, ir) (BK input).
void pack_z_plane(const Complex* psi, const DeviceLayout& layout, int nr, int iz, Complex* plane);
// Host layout <-> device layout for the points izr_first .. izr_first + n_points - 1: `chunk`
// holds them in host layout (irot, ixy, izr - izr_first).
void host_chunk_to_device(const Complex* chunk, const DeviceLayout& layout, std::int64_t izr_first,
                          std::int64_t n_points, Complex* psi);
void device_to_host_chunk(const Complex* psi, const DeviceLayout& layout, std::int64_t izr_first, std::int64_t n_points,
                          Complex* chunk);

}  // namespace spodvr::cuda
