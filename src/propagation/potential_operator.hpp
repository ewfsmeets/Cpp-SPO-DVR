#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "core/aligned_vector.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"
#include "symmetry/grid_symmetry.hpp"

namespace spodvr {

// The per-(irot, ir) factors of the potential step, each merged into one multiplication:
//   before_interaction: exp(-V_opt,r dt) rot(dt/2)                          (Z < zasymp)
//   gas_phase_step:     exp(-V_opt,r dt) rot(dt/2) exp(-i vasy dt) rot(dt/2)  (Z >= zasymp)
// Layout irot + n_rot * ir. Fortran: Absorption.f90 PerformEXPOpticalR, GasPhase.f90
// PerformEXPGasPhaseOnScatGrid, EXPKinEn%halfStep%RotChan
struct PotentialStepFactors {
  std::vector<std::complex<double>> before_interaction;
  std::vector<std::complex<double>> gas_phase_step;
};
[[nodiscard]] PotentialStepFactors potential_step_factors(const ScatteringShape& shape,
                                                          const PropagationTables& tables);

// exp(i A) as the Fortran code forms it: B = COS(A) + ima*SIN(A).
[[nodiscard]] inline std::complex<double> exp_i(double a) { return {std::cos(a), std::sin(a)}; }

// exp(-i V dt) on the scattering grid, per (Z, r) point:
//   optical potentials in Z and r;
//   for Z < zasymp (nzpot points): add the specular wave to the initial channel, then
//     rot(dt/2) . [FBR -> DVR . exp(-i V dt) . DVR -> FBR] . rot(dt/2),
//     then project the initial-channel part back onto the specular grid and subtract it;
//   beyond: rot(dt/2) . exp(-i vasy(r) dt) . rot(dt/2).
// Every (Z, r) point is independent (the specular wave is only touched at its own point).
// Fortran: Potential.f90 PerformEXPPotentialOp, PerformEXPInteractionOp; Absorption.f90
// PerformEXPOpticalZ/R; TwoGrids.f90 AddSpecToScat, SubSpecFromScat
class PotentialOperator {
 public:
  // Default limit for the precomputed exp(i A) table (see below).
  static constexpr std::size_t default_table_limit = std::size_t(4) << 30;  // 4 GiB

  // premultiplied_potential: A = -V dt on (itp, irep, izr) for the nzpot*nr points of the
  // interaction region (PotentialGrid after premultiply_with_time_step); irep runs over the
  // (x, y) representatives of symmetry, whose maps give A at every (x, y, phi) point.
  // exp(i A) = (cos A, sin A) is precomputed once (for the representatives) if that table needs
  // at most table_limit bytes; otherwise it is evaluated in every step (as the Fortran code does).
  // Both give bitwise identical results. Only in the second case is premultiplied_potential used
  // after construction (it must then stay alive).
  PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                    const PropagationTables& tables, const AngularTransform& transform,
                    std::span<const double> premultiplied_potential, GridSymmetry symmetry,
                    std::size_t table_limit = default_table_limit);

  // Without surface symmetry (A on the full (itp, ixy, izr) grid).
  PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                    const PropagationTables& tables, const AngularTransform& transform,
                    std::span<const double> premultiplied_potential, std::size_t table_limit = default_table_limit);

  [[nodiscard]] bool uses_precomputed_exponentials() const { return !exp_potential_.empty(); }
  // Size of the work array that apply_block needs (one per thread).
  [[nodiscard]] std::size_t work_size() const;

  // All (Z, r) points.
  void apply(WavePacket& psi) const;
  // One (Z, r) point izr: block = its (irot, ixy) values; work: work_size() elements.
  void apply_block(int izr, std::span<std::complex<double>> block, SpecularWave& specular,
                   std::span<std::complex<double>> work) const;

 private:
  void apply_interaction(std::span<std::complex<double>> block, int izr, std::span<std::complex<double>> work) const;
  // exp(i A) of the representatives of izr into scratch (on-the-fly path with symmetry).
  void exponentials_of_block(int izr, std::span<std::complex<double>> scratch) const;
  [[nodiscard]] std::size_t dvr_work_size() const;
  [[nodiscard]] bool reduces() const { return symmetry_.n_representatives() < symmetry_.n_xy(); }
  void apply_rotation_half_step(std::span<std::complex<double>> block, int ir) const;

  ScatteringShape shape_;
  int nz_potential_;
  int initial_channel_;  // Fortran: ichn - 1
  const PropagationTables& tables_;
  const AngularTransform& transform_;
  PotentialStepFactors factors_;  // (irot, ir), see potential_step_factors
  std::span<const double> potential_;              // A = -V dt (when not precomputed)
  ComplexVector exp_potential_;                    // exp(i A) (when precomputed)
  GridSymmetry symmetry_;                          // (x, y) representatives and phi maps
};

}  // namespace spodvr
