#pragma once

namespace spodvr {

// Where an absorber's reduced length L/λ lies relative to the Vibok & Balint-Kurti table
// (L/λ = 1 .. 15). Below or above it, the strength is the table's value at the nearer end; at
// exactly 1 no interval of the table matches, and the strength is 0 (absorber D2).
enum class AbsorberTableRange { inside, below, above, at_first_point };

// An absorber's reduced length L/λ (λ = 2π/k, the de Broglie wavelength at the energy the
// absorber is tuned for) and where it lies relative to the table.
struct AbsorberReducedLength {
  double value;              // L/λ. Fortran: reduced_length
  AbsorberTableRange range;  // below: L/λ < 1, above: L/λ > 15 (the Fortran warns for both)
};

// L/λ for a kinetic energy ekin (hartree), an absorber length (bohr) and a mass (m_e).
// Fortran: MyFunctions.f90 getACoeffOptPotential (reduced_length and the branch it selects)
[[nodiscard]] AbsorberReducedLength absorber_reduced_length(double ekin, double length, double mass);

// Strength A of the quadratic optical potential for a kinetic energy ekin (hartree) and the
// absorber's reduced length (from absorber_reduced_length), interpolated from the table.
// Fortran: MyFunctions.f90 getACoeffOptPotential
[[nodiscard]] double optical_potential_strength(double ekin, const AbsorberReducedLength& reduced);

// The same for an absorber length (bohr) and a mass (m_e).
[[nodiscard]] double optical_potential_strength(double ekin, double length, double mass);

// Width parameter alpha of the Gaussian wave packet in k covering [low, high] (hartree).
// Fortran: MyFunctions.f90 getAlphaCoeffWavePacket
[[nodiscard]] double wave_packet_alpha(double low_energy, double high_energy, double mass);

// Energy (hartree) at the centre k0 of the wave packet. Fortran: getAverageEnergyWavePacket
[[nodiscard]] double wave_packet_average_energy(double low_energy, double high_energy, double mass);

}  // namespace spodvr
