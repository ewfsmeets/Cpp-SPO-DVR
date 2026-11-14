// Optical-potential strength and Gaussian wave-packet parameters. Fortran: MyFunctions.f90
#include "input/wave_packet_params.hpp"

#include <array>
#include <cmath>

namespace spodvr {

// Fortran: MyFunctions.f90 getACoeffOptPotential
double optical_potential_strength(double ekin, double length, double mass) {
  static constexpr std::array<double, 29> xtab = {1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5,
                                                  6.0,  6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5,  10.0, 10.5,
                                                  11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0};
  static constexpr std::array<double, 29> ftab = {2.750, 2.640, 2.510, 2.380, 2.230, 2.090, 1.940, 1.800, 1.670, 1.550,
                                                  1.440, 1.340, 1.250, 1.170, 1.100, 1.040, 0.986, 0.949, 0.903, 0.870,
                                                  0.840, 0.813, 0.788, 0.764, 0.741, 0.718, 0.697, 0.678, 0.663};

  const double momentum = std::sqrt(2.0 * mass * ekin);
  const double wavelength = 2.0 * std::acos(-1.0) / momentum;
  const double reduced_length = length / wavelength;

  double asel = 0.0;
  if (reduced_length < xtab.front()) {
    asel = ftab.front() * ekin;  // Fortran warns: extrapolating; reduced length too small
  } else if (reduced_length > xtab.back()) {
    asel = ftab.back() * ekin;  // Fortran warns: extrapolating; reduced length too big
  } else {
    for (std::size_t i = 0; i + 1 < xtab.size(); ++i) {
      if (reduced_length > xtab[i] && reduced_length <= xtab[i + 1]) {
        asel = (ftab[i] + (reduced_length - xtab[i]) * (ftab[i + 1] - ftab[i]) / (xtab[i + 1] - xtab[i])) * ekin;
        break;
      }
    }
    // Note: reduced_length == xtab[0] exactly matches no interval; the Fortran result is then
    // its uninitialised local (here 0).
  }
  return asel;
}

// Fortran: MyFunctions.f90 getAlphaCoeffWavePacket
double wave_packet_alpha(double low_energy, double high_energy, double mass) {
  const double kmin = std::sqrt(2.0 * mass * low_energy);
  const double kmax = std::sqrt(2.0 * mass * high_energy);
  const double k0 = (kmin + kmax) / 2.0;
  return std::pow(k0 - kmin, 2) / 4.0;
}

// Fortran: MyFunctions.f90 getAverageEnergyWavePacket
double wave_packet_average_energy(double low_energy, double high_energy, double mass) {
  const double kmin = -std::sqrt(2.0 * mass * low_energy);
  const double kmax = -std::sqrt(2.0 * mass * high_energy);
  const double k0 = (kmin + kmax) / 2.0;
  return std::pow(k0, 2) / (2.0 * mass);
}

}  // namespace spodvr
