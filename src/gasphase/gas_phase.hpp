// Rovibrational states of the free molecule (the asymptotic channels). Fortran: GasPhase.f90
#pragma once

#include <span>
#include <string>
#include <utility>
#include <vector>

#include "fft/fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "linalg/symmetric_eigensolver.hpp"

namespace spodvr {

// Rovibrational eigenstates of the free molecule for j = 0 .. JMax and all nr vibrational
// levels. Wave functions are stored in the Fortran order (ir, j, v). Atomic units.
// Fortran: GasPhase.f90 eivaa(0:jmaxw, nr), eivea(nr, 0:jmaxw, nr)
class RovibrationalStates {
 public:
  RovibrationalStates(int n_r, int j_max)
      : n_r_(n_r), j_max_(j_max), energies_(n_r * (j_max + 1)), wave_functions_(n_r * n_r * (j_max + 1)) {}

  [[nodiscard]] int n_r() const { return n_r_; }
  [[nodiscard]] int j_max() const { return j_max_; }
  [[nodiscard]] int n_v() const { return n_r_; }  // one vibrational level per r grid point

  [[nodiscard]] double energy(int v, int j) const { return energies_[state_index(v, j)]; }
  double& energy(int v, int j) { return energies_[state_index(v, j)]; }

  [[nodiscard]] std::span<const double> wave_function(int v, int j) const {
    return std::span<const double>(wave_functions_).subspan(n_r_ * state_index(v, j), n_r_);
  }
  std::span<double> wave_function(int v, int j) {
    return std::span<double>(wave_functions_).subspan(n_r_ * state_index(v, j), n_r_);
  }

 private:
  [[nodiscard]] std::size_t state_index(int v, int j) const { return j + (j_max_ + 1) * v; }

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

// Expectation values of the initial rovibrational state (vin, jin), hartree / bohr.
// Fortran: eviroi, vpotin, erotin, evibin, ekvroi, rinvj
struct InitialStateProperties {
  double energy = 0.0;
  double potential_energy = 0.0;
  double rotational_energy = 0.0;
  double vibrational_energy = 0.0;               // energy - potential - rotational
  double rovibrational_kinetic_energy = 0.0;     // vibrational + rotational
  double r_mean = 0.0;
};

struct GasPhase {
  explicit GasPhase(RovibrationalStates solved_states) : states(std::move(solved_states)) {}

  RovibrationalStates states;
  int v_max_open = 0;       // number of open vibrational levels. Fortran: ivmaxAnalysis
  int j_max_open = 0;       // highest open j. Fortran: jmaxAnalysis
  int n_v_analysis = 0;     // number of vibrational levels analysed. Fortran: nviba
  int j_max_analysis = 0;   // highest j analysed. Fortran: jana
  std::vector<RovibrationalStateSummary> summary;  // v < n_v_analysis, j = 0 .. JMax
  InitialStateProperties initial;

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
