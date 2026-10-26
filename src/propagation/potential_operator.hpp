#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

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

  // premultiplied_potential: A = -V dt on (itp, ixy, izr) for the nzpot*nr points of the
  // interaction region (PotentialGrid after premultiply_with_time_step).
  // exp(i A) = (cos A, sin A) is precomputed once if that table needs at most table_limit
  // bytes; otherwise it is evaluated in every step (as the Fortran code does). Both give
  // bitwise identical results. Only in the second case is premultiplied_potential used
  // after construction (it must then stay alive).
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
  void apply_interaction(std::span<std::complex<double>> block, int izr, std::span<std::complex<double>> dvr) const;
  void apply_rotation_half_step(std::span<std::complex<double>> block, int ir) const;

  ScatteringShape shape_;
  int nz_potential_;
  int initial_channel_;  // Fortran: ichn - 1
  const PropagationTables& tables_;
  const AngularTransform& transform_;
  std::span<const double> potential_;              // A = -V dt (when not precomputed)
  std::vector<std::complex<double>> exp_potential_;  // exp(i A) (when precomputed)
};

}  // namespace spodvr
