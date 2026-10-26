#pragma once

#include <optional>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// The complete state of the analysis surface, for checkpoints. 0-based indices.
struct AnalysisSurfaceState {
  int index;            // Fortran: izanlySpecChan - 1
  int previous_index;   // Fortran: oldValZanly - 1
  bool fixed;           // Fortran: FixedZAnalysisSpec
  double norm_tolerance;  // Fortran: anlyTol
};

// Position of the analysis surface on the specular grid (ProjOp = 2). With DynamZSpAna it
// starts just before the specular absorber and moves inwards (towards the surface) as long as
// the wave packet still has significant norm beyond it; once it can no longer move it is fixed.
// The rest norm of the averages is summed up to this point. Indices are 0-based.
// Fortran: Analyzor.f90 SetupAnalysis (izanlySpecChan, anlyTol), FindNewAnalysisPoint
class AnalysisSurface {
 public:
  // initial: the initial specular wave function (used for the dynamic tolerance).
  AnalysisSurface(const Config& config, const SpatialGrids& grids, const SpecularWave& initial);

  [[nodiscard]] int index() const { return index_; }   // Fortran: izanlySpecChan - 1
  [[nodiscard]] bool fixed() const { return fixed_; }  // Fortran: FixedZAnalysisSpec
  [[nodiscard]] double tolerance() const { return tolerance_; }  // Fortran: anlyTol * delz * delr

  // Moves the surface for the current wave function (no-op once fixed). Returns the Z value
  // if the surface became fixed during this call. Fortran: FindNewAnalysisPoint
  std::optional<double> update(const SpecularWave& specular);

  [[nodiscard]] AnalysisSurfaceState state() const { return {index_, previous_index_, fixed_, norm_tolerance_}; }
  // Continues from a saved state (restart).
  void restore(const AnalysisSurfaceState& state);

 private:
  const UniformGrid& z_;
  int iz_analysis_;     // Fortran: izanly - 1
  int last_before_absorber_;  // Fortran: nschmi - 1
  double norm_tolerance_;     // tolerance on sum |psi|^2 (without dZ dR). Fortran: anlyTol
  double tolerance_;
  bool fixed_;
  int index_;
  int previous_index_;  // Fortran: oldValZanly - 1
  double volume_;       // dZ dR
};

}  // namespace spodvr
