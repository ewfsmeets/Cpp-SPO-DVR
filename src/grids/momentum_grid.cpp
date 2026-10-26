// Momentum grids and kinetic energies of the FFT coordinates. Fortran: Setup.f90
#include "grids/momentum_grid.hpp"

#include <cmath>

#include "core/constants.hpp"

namespace spodvr {

// Fortran: Setup.f90 MiscellaneousSetup, "kinetic energy in Z" / "kinetic energies in r" loops
MomentumGrid::MomentumGrid(int n, double coordinate_spacing, double mass)
    : dk_(2.0 * Constants::pi / (n * coordinate_spacing)), k_(n, 0.0), ekin_(n, 0.0) {
  // k = 0 at index 0, +k at index i and -k at index n-i, and the Nyquist point +n/2 dk.
  const int half = n / 2;
  ekin_[half] = 0.5 / mass * std::pow(half * dk_, 2);
  k_[half] = half * dk_;
  for (int i = 1; i < half; ++i) {  // Fortran i = 2 .. n/2
    ekin_[i] = 0.5 / mass * std::pow(i * dk_, 2);
    k_[i] = i * dk_;
    ekin_[n - i] = ekin_[i];
    k_[n - i] = -k_[i];
  }
}

}  // namespace spodvr
