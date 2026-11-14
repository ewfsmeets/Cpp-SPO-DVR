#pragma once

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

}  // namespace spodvr
