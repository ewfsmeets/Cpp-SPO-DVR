#pragma once

#include <memory>

#include "core/aligned_vector.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"
#include "propagation/propagation_backend.hpp"
#include "propagation/specular_operator.hpp"
#include "propagation/split_operator.hpp"

namespace spodvr {

// The wave packet in host memory, propagated with the OpenMP operator stages
// (SplitOperatorPropagator). Fortran: SplitOperator.f90
class CpuBackend final : public PropagationBackend {
 public:
  // potential: premultiplied by -dt here (Fortran: vdvr = -V dt); initial: the specular wave of
  // the initial wave packet (the scattering grid starts empty).
  CpuBackend(const BackendInput& input, PotentialGrid potential, SpecularWave initial);

  void read_wave_packet(const std::function<void(WavePacket&)>& read) override { read(psi_); }
  [[nodiscard]] const WavePacket& wave_packet_on_host() override { return psi_; }
  void propagate_analysis_step() override { propagator_.propagate_analysis_step(psi_); }
  [[nodiscard]] const SpecularWave& specular() const override { return psi_.specular; }
  [[nodiscard]] std::span<const std::complex<double>> z_plane(int iz) const override;
  [[nodiscard]] std::vector<double> block_densities() const override;

 private:
  // Declaration order = construction order: the operators reference the members above them.
  AngularTransform transform_;
  AlignedVector<double> premultiplied_potential_;  // -V dt. Fortran: vdvr
  ScatteringKineticOperator kinetic_;
  PotentialOperator potential_;
  SpecularOperator specular_;
  SplitOperatorPropagator propagator_;
  WavePacket psi_;
};

// The CPU backend as a BackendFactory.
[[nodiscard]] std::unique_ptr<PropagationBackend> make_cpu_backend(const BackendInput& input, PotentialGrid potential,
                                                                   SpecularWave initial);

}  // namespace spodvr
