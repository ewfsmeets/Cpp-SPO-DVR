#include "potential/hard_wall_potential.hpp"

#include <cmath>

#include "core/constants.hpp"

namespace spodvr {

namespace {
constexpr double z_asymptotic = 7.0;  // Fortran: ZAsymptotic
}

// Fortran: HardWallPES.f90 module parameters and getParamsForHardwallPES
HardWallPotential::HardWallPotential()
    : parameters_{6.239660, 5.239660, z_asymptotic, Constants::pi / 2.0, "C2v"} {}

// Fortran: HardWallPES.f90 getVHardwallPESAsymptotic. The expression (including the real
// exponents r**2.0, r**3.0) is kept literally so that the result is bitwise identical.
double HardWallPotential::gas_phase(double r_in) {
  const double d = 4.8286;  // H2 Rydberg test parameters
  const double a1 = 2.282;
  const double a2 = 1.555;
  const double a3 = 0.7533;
  const double a4 = 2.23;
  const double re = 1.40;

  const double r = r_in - re;
  return -d * (1.0 + a1 * r + a2 * std::pow(r, 2.0) + a3 * std::pow(r, 3.0)) * std::exp(-a4 * r) + d;
}

// Fortran: HardWallPES.f90 getVHardwallPES
double HardWallPotential::value(double /*x*/, double /*y*/, double z, double r,
                                double /*theta*/, double /*phi*/, int /*iz*/, int /*ir*/) const {
  const double vr = gas_phase(r);
  const double vz = 1.0 / (z * z * z * z + 0.001);

  const double z_max = z_asymptotic;
  const double z_min = 0.1778;

  double v = vr;
  if (z < z_max) {
    if (z < z_min) {
      v = vr + vz;
    } else {
      v = vr + (1.0 - (0.5 - 0.5 * std::cos(Constants::pi * (z - z_min) / (z_max - z_min)))) * vz;
    }
  }
  return v;
}

}  // namespace spodvr
