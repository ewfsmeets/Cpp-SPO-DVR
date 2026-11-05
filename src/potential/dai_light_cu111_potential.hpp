#pragma once

#include "potential/potential.hpp"

namespace spodvr {

// H2 on Cu(111): the analytic LEPS potential of Dai and Light, J. Chem. Phys. 102, 6280 (1995).
// Atom-surface Morse-type terms whose parameters are Fourier-interpolated between the top, bridge
// and centre sites of the hexagonal cell (60 deg skewed coordinates), plus the H2 molecular term.
// For Z >= zasymp the potential is the gas-phase LEPS molecular potential.
// Fortran: PESForCu111DaiLight.f90 (SystemID "H2Cu111DaiLight")
class DaiLightCu111Potential final : public Potential {
 public:
  DaiLightCu111Potential();

  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override;
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "H2Cu111DaiLight"; }

 private:
  SystemParameters parameters_;
};

}  // namespace spodvr
