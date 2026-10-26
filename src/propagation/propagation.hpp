#pragma once

#include <functional>
#include <vector>

#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/averages.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/potential_operator.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/specular_operator.hpp"
#include "propagation/split_operator.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// A new wave-packet propagation (ProjOp = 2): sets up the initial wave packet and all
// operators, then runs the time loop of Main.f90 (nread = 0).
// Fortran: Main.f90 (SetupRDBasis ... SetupSplitOperator, DO propStep loop), Averages.f90
class WavePacketPropagation {
 public:
  // The potential grid is taken over and premultiplied by -dt.
  WavePacketPropagation(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                        PotentialGrid potential, const GasPhase& gas, const Fft& fft);

  [[nodiscard]] const RotationalBasis& rotational_basis() const { return basis_; }
  [[nodiscard]] const InitialWavePacketInfo& initial_wave_packet() const { return initial_.info; }
  [[nodiscard]] const AnalysisSurface& analysis_surface() const { return surface_; }
  [[nodiscard]] const PropagationTables& tables() const { return tables_; }
  [[nodiscard]] const WavePacket& wave_packet() const { return psi_; }

  // Propagates from tStart to tEnd. on_averages is called at tStart and every nave analysis
  // steps; on_surface_fixed(Z) when the dynamic analysis surface becomes fixed.
  void run(const std::function<void(const Averages&)>& on_averages,
           const std::function<void(double)>& on_surface_fixed);

 private:
  const Config& config_;
  const SpatialGrids& grids_;
  RotationalBasis basis_;
  AngularTransform transform_;
  PropagationTables tables_;
  std::vector<double> premultiplied_potential_;  // -V dt. Fortran: vdvr
  ScatteringKineticOperator kinetic_;
  PotentialOperator potential_;
  SpecularOperator specular_;
  SplitOperatorPropagator propagator_;
  InitialWavePacket initial_;
  WavePacket psi_;
  AnalysisSurface surface_;
};

}  // namespace spodvr
