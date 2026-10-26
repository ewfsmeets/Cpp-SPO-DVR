#pragma once

#include <span>

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
  // premultiplied_potential: -V dt on (itp, ixy, izr) for the nzpot*nr points of the
  // interaction region (PotentialGrid after premultiply_with_time_step).
  PotentialOperator(const ScatteringShape& shape, int nz_potential, int initial_channel,
                    const PropagationTables& tables, const AngularTransform& transform,
                    std::span<const double> premultiplied_potential);

  void apply(WavePacket& psi) const;

 private:
  void apply_interaction(std::span<std::complex<double>> block, int izr, std::span<std::complex<double>> dvr) const;
  void apply_rotation_half_step(std::span<std::complex<double>> block, int ir) const;

  ScatteringShape shape_;
  int nz_potential_;
  int initial_channel_;  // Fortran: ichn - 1
  const PropagationTables& tables_;
  const AngularTransform& transform_;
  std::span<const double> potential_;
};

}  // namespace spodvr
