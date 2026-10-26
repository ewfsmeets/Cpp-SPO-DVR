#pragma once

#include <span>
#include <string>
#include <vector>

#include "fft/fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "linalg/symmetric_eigensolver.hpp"

namespace spodvr {

// Rovibrational eigenstates of the free molecule for j = 0 .. JMax and all nr vibrational
// levels. Wave functions are stored in the Fortran order (ir, j, iv). Atomic units.
// Fortran: GasPhase.f90 eivaa(0:jmaxw, nr), eivea(nr, 0:jmaxw, nr)
class RovibrationalStates {
 public:
  RovibrationalStates(int n_r, int j_max)
      : n_r_(n_r), j_max_(j_max),
        energies_(static_cast<std::size_t>(n_r) * static_cast<std::size_t>(j_max + 1)),
        wave_functions_(energies_.size() * static_cast<std::size_t>(n_r)) {}

  [[nodiscard]] int n_r() const { return n_r_; }
  [[nodiscard]] int j_max() const { return j_max_; }
  [[nodiscard]] int n_v() const { return n_r_; }

  [[nodiscard]] double energy(int v, int j) const { return energies_[energy_index(v, j)]; }
  double& energy(int v, int j) { return energies_[energy_index(v, j)]; }
  [[nodiscard]] std::span<const double> wave_function(int v, int j) const {
    return std::span<const double>(wave_functions_).subspan(wave_function_offset(v, j), static_cast<std::size_t>(n_r_));
  }
  std::span<double> wave_function(int v, int j) {
    return std::span<double>(wave_functions_).subspan(wave_function_offset(v, j), static_cast<std::size_t>(n_r_));
  }

 private:
  [[nodiscard]] std::size_t energy_index(int v, int j) const {
    return static_cast<std::size_t>(j) + static_cast<std::size_t>(j_max_ + 1) * static_cast<std::size_t>(v);
  }
  [[nodiscard]] std::size_t wave_function_offset(int v, int j) const {
    return static_cast<std::size_t>(n_r_) * energy_index(v, j);
  }

  int n_r_;
  int j_max_;
  std::vector<double> energies_;
  std::vector<double> wave_functions_;
};

// One row of the Fortran "asymptotic molecular vibrational wave functions" table.
struct RovibrationalStateSummary {
  int v;
  int j;
  double energy;            // hartree
  double relative_energy;   // E(v,j) - E(vin,jin), hartree
  double psi_first;         // psi(r_1)
  double psi_last;          // psi(r_nr)
  double norm;              // sum psi^2 * dr
  double r_mean;            // sum r psi^2 * dr
};

struct GasPhase {
  RovibrationalStates states;
  int n_v_analysis;         // number of vibrational levels analysed. Fortran: nviba (ivmaxAnalysis)
  int j_max_analysis;       // highest j analysed. Fortran: jana (jmaxAnalysis)
  int v_max_open;           // Fortran: ivmaxAnalysis (count of open vibrational levels)
  int j_max_open;           // Fortran: jmaxAnalysis
  std::vector<RovibrationalStateSummary> summary;  // v < n_v_analysis, j = 0 .. JMax

  // Initial state (vin, jin). Fortran: eviroi, vpotin, erotin, evibin, ekvroi, rinvj
  double initial_energy;
  double initial_potential_energy;
  double initial_rotational_energy;
  double initial_vibrational_energy;
  double initial_rovibrational_kinetic_energy;
  double initial_r_mean;

  std::vector<double> vibrational_energies;  // E(v,0) - E(0,0), v < n_v_analysis. Fortran: aevib
  std::vector<double> rotational_energies;   // E(0,j) - E(0,0), j = 0..JMax. Fortran: aerot
  std::vector<std::string> warnings;
};

// Fortran: GasPhase.f90 SetupAsympWF. asymptotic_potential is vasy (hartree) on the r grid.
[[nodiscard]] GasPhase compute_gas_phase(const Config& config, const SpatialGrids& grids,
                                         const AngularGrid& angular,
                                         std::span<const double> asymptotic_potential, const Fft& fft,
                                         const SymmetricEigensolver& solver);

}  // namespace spodvr
