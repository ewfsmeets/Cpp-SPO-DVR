#pragma once

#include "potential/potential.hpp"

namespace spodvr {

// Hard wall test potential: H2 gas-phase Rydberg potential in r plus a 1/(Z^4 + 0.001) wall in
// Z that is switched off (cosine switch) between Z = 0.1778 and Z = zasymp.
// Fortran: HardWallPES.f90 (the __EXP_HARD_WALL__ variant is not used by the active build)
class HardWallPotential final : public Potential {
 public:
  HardWallPotential();

  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override;
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "HardWallPES"; }

  // The gas-phase (Rydberg) potential V(r) in eV. Fortran: getVHardwallPESAsymptotic
  [[nodiscard]] static double gas_phase(double r);

 private:
  SystemParameters parameters_;
};

}  // namespace spodvr
