#pragma once

#include <complex>

namespace spodvr {

// Strength A of the quadratic optical potential for a kinetic energy ekin (hartree), an
// absorber length (bohr) and a mass (m_e), interpolated from the Vibok & Balint-Kurti table.
// Fortran: MyFunctions.f90 getACoeffOptPotential
[[nodiscard]] double optical_potential_strength(double ekin, double length, double mass);

// Width parameter alpha of the Gaussian wave packet in k covering [low, high] (hartree).
// Fortran: MyFunctions.f90 getAlphaCoeffWavePacket
[[nodiscard]] double wave_packet_alpha(double low_energy, double high_energy, double mass);

// Energy (hartree) at the centre k0 of the wave packet. Fortran: getAverageEnergyWavePacket
[[nodiscard]] double wave_packet_average_energy(double low_energy, double high_energy, double mass);

// Momentum amplitude b(k) of the Gaussian wave packet covering [low, high] (hartree), centred at
// z0 (bohr): (2 xi2 / pi)^(1/4) exp(-xi2 (k0 - k)^2 + i (k0 - k) z0), xi2 = 1 / (4 alpha), so that
// the integral of |b|^2 over k is 1. The phase has the opposite sign of the packet's (flux Q9).
// Fortran: MyFunctions.f90 getBkCoeffWavePacket
[[nodiscard]] std::complex<double> wave_packet_momentum_amplitude(double momentum, double low_energy,
                                                                  double high_energy, double mass, double z0);

}  // namespace spodvr
