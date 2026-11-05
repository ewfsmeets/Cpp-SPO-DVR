#pragma once

#include "potential/potential.hpp"

namespace spodvr {

// H2 on Cu(111): the corrugated LEPS potential of Persson et al., J. Chem. Phys. 110, 2240 (1999).
// Works in Angstrom and eV internally; the atom-surface parameters are Fourier expansions over
// the hexagonal cell. Note (Fortran quirk, kept): coordinates are converted to Angstrom with
// Persson's factor 0.529167 while the lattice constant uses cboan = 0.529177, so the potential
// is not exactly periodic on the grid cell (period mismatch 1.9e-5 relative).
// Fortran: PESForCu111Persson.f90 (SystemID "H2Cu111Persson")
class PerssonCu111Potential final : public Potential {
 public:
  PerssonCu111Potential();

  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override;
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "H2Cu111Persson"; }

 private:
  SystemParameters parameters_;
};

}  // namespace spodvr
