#include "grids/momentum_grid.hpp"

#include <cmath>

#include "core/constants.hpp"

namespace spodvr {

// Fortran: Setup.f90 MiscellaneousSetup, "kinetic energy in Z" / "kinetic energies in r" loops
MomentumGrid::MomentumGrid(int n, double coordinate_spacing, double mass)
    : dk_(2.0 * Constants::pi / (static_cast<double>(n) * coordinate_spacing)),
      k_(static_cast<std::size_t>(n), 0.0),
      ekin_(static_cast<std::size_t>(n), 0.0) {
  const auto half = static_cast<std::size_t>(n / 2);
  const auto nn = static_cast<std::size_t>(n);
  ekin_[0] = 0.0;
  k_[0] = 0.0;
  ekin_[half] = 0.5 / mass * std::pow(static_cast<double>(half) * dk_, 2);
  k_[half] = static_cast<double>(half) * dk_;
  for (std::size_t i = 1; i < half; ++i) {  // Fortran i = 2 .. n/2
    ekin_[i] = 0.5 / mass * std::pow(static_cast<double>(i) * dk_, 2);
    k_[i] = static_cast<double>(i) * dk_;
    ekin_[nn - i] = ekin_[i];
    k_[nn - i] = -k_[i];
  }
}

}  // namespace spodvr
