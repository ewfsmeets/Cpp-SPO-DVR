#pragma once

#include <complex>
#include <span>
#include <string>
#include <vector>

#include "input/config.hpp"

namespace spodvr {

// Number of flux analysis energies, INT((EMaxFluxAna - EMinFluxAna) / DeltaEFluxAna) + 1 with the
// energies in hartree: the end point is dropped for most inputs (flux Q3). Shared by
// make_flux_energies and the memory estimate, so the two cannot disagree. 0 if doFluxAnalysis is off.
// Throws InputError (DeltaEFluxAna) above max_flux_energies.
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (numScatEnergies)
[[nodiscard]] int flux_energy_count(const InputParameters& input);

// Upper limit of the number of flux energies: far above any count whose coefficients fit in memory
// (each energy holds two complex values per point of the flux surface), and safe for int.
inline constexpr int max_flux_energies = 1'000'000;

// The energies of the flux analysis and the per-energy factors of its time-to-energy transform,
// in atomic units. The Fortran recomputes k_n and b(k_n) at every analysis step from the same
// operands; here they are computed once (the same values).
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor, storeTimeIndependentWF, computeFluxProbability
struct FluxEnergies {
  std::vector<double> z_energy;  // E_n = E_min + n dE: the initial state's energy in Z. Fortran: ZScatEnergy
  std::vector<double> total;     // (E_n + parallel energy) + rovibrational energy. Fortran: totalScatEnergy
  std::vector<double> momentum;  // k_n = sqrt(2 M E_n). Fortran: kZ_initial
  std::vector<std::complex<double>> amplitude;  // b(k_n) of the initial wave packet. Fortran: bkz
  std::vector<std::string> warnings;            // EMaxFluxAna not reached (flux Q3)

  [[nodiscard]] int size() const { return static_cast<int>(z_energy.size()); }
};

// The flux energies of the input: EMinFluxAna, EMaxFluxAna, DeltaEFluxAna and the wave packet's
// EMinWP, EMaxWP converted with units::precise (flux Q5), M = total mass, b(k) centred at Z0.
// parallel_energy: PropagationTables::parallel_energy (Fortran edifin); initial_state_energy:
// GasPhase::initial.energy (Fortran eviroi); both hartree.
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor, storeTimeIndependentWF
[[nodiscard]] FluxEnergies make_flux_energies(const Config& config, double parallel_energy,
                                              double initial_state_energy);

// Weights D(ir), ir = 0 .. nr-1, of the Fourier-grid derivative d psi / dr at the flux surface
// r_c = r(ir_flux), with ir_flux 0-based (Fortran irFluxAnalysis - 1) and dr in bohr:
// D(ir) = ((2c) x) / dr with c = pi / nr, d = ir_flux - ir and x = (-1)^d / (2 tan(c d)); D(ir_flux) = 0.
// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr (ddr(:, irFluxAnalysis))
[[nodiscard]] std::vector<double> flux_derivative_weights(int nr, double dr, int ir_flux);

// d psi / dr at the flux surface from one r column (column[ir] = psi(r_ir)) and the weights above:
// the sum in the Fortran order, from 0 over every ir ascending, the zero term D(ir_flux) included.
// The reference of the backend tests (CPU and CUDA), which compare with it by EXPECT_EQ: the backends
// loop ir-outer over the strided wave function, with the same order per element.
// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr
[[nodiscard]] std::complex<double> flux_r_derivative(std::span<const double> weights,
                                                     std::span<const std::complex<double>> column);

}  // namespace spodvr
