// Rovibrational states of the free molecule and the open-channel limits used by the analysis.
// Fortran: GasPhase.f90 SetupAsympWF
#include "gasphase/gas_phase.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/machine_epsilon.hpp"
#include "core/unit_conversion.hpp"
#include "gasphase/fourier_grid_hamiltonian.hpp"

namespace spodvr {

namespace {

// Solves the FGH eigenproblem for every j with the centrifugal term j(j+1)/(2 mu r^2) added to
// the gas-phase potential. Fortran: SetupAsympWF loop over Jlevel, CopyAsympPotentialToFGHType
RovibrationalStates solve_rovibrational_states(const FourierGridHamiltonian& fgh, int j_max,
                                               std::span<const double> asymptotic_potential,
                                               double reduced_mass, const SymmetricEigensolver& solver) {
  const int nr = fgh.size();
  RovibrationalStates states(nr, j_max);
  std::vector<double> potential(nr);

  for (int j = 0; j <= j_max; ++j) {
    for (int i = 0; i < nr; ++i) {
      if (j == 0) {
        potential[i] = asymptotic_potential[i];
        continue;
      }
      double r = fgh.grid()[i];
      if (r == 0.0) r = r + fortran_machine_epsilon();  // prevent division by zero
      potential[i] = asymptotic_potential[i] + (j * (j + 1)) / 2.0 / reduced_mass / (r * r);
    }

    const EigenDecomposition eigen = fgh.eigenstates(potential, solver);
    for (int v = 0; v < nr; ++v) {
      states.energy(v, j) = eigen.values[v];
      std::ranges::copy(eigen.vector(v), states.wave_function(v, j).begin());
    }
  }
  return states;
}

// Number of open vibrational levels: count v (for j = 0) until E(v,0) - E(0,0) exceeds the
// available energy. Fortran: SetupAsympWF (E. Pijper part), ivmaxAnalysis
int count_open_vibrational_levels(const RovibrationalStates& states, double available_energy,
                                  std::vector<std::string>& warnings) {
  const double e00 = states.energy(0, 0);
  int iv = 1;  // 1-based, as in the Fortran loop
  while (states.energy(iv - 1, 0) - e00 <= available_energy) {
    iv = iv + 1;
    if (iv == states.n_v()) {
      warnings.emplace_back("Vibrational basis too small to describe all open channels!");
      break;
    }
  }
  return iv - 1;
}

// Highest open j (for v = 0), stepping by delta_j (2 for a homonuclear molecule).
// Fortran: SetupAsympWF (E. Pijper part), jmaxAnalysis
int highest_open_rotational_level(const RovibrationalStates& states, int delta_j, double available_energy,
                                  std::vector<std::string>& warnings) {
  const double e00 = states.energy(0, 0);
  int j = 0;
  while (states.energy(0, j) - e00 <= available_energy) {
    j = j + delta_j;
    if (j > states.j_max()) {
      warnings.emplace_back("Rotational basis too small to describe all open channels!");
      break;
    }
  }
  return j - delta_j;
}

// Norm and <r> of one state with the scattering r grid spacing (as in the Fortran table).
RovibrationalStateSummary summarise_state(const RovibrationalStates& states, int v, int j,
                                          const UniformGrid& r_grid, double e_initial) {
  const auto psi = states.wave_function(v, j);
  double norm = 0.0;
  double r_mean = 0.0;
  for (int ir = 0; ir < r_grid.size(); ++ir) {
    const double y = psi[ir];
    norm = norm + y * y;
    r_mean = r_mean + r_grid[ir] * y * y;
  }
  norm = norm * r_grid.spacing();
  r_mean = r_mean * r_grid.spacing();
  return {v, j, states.energy(v, j), states.energy(v, j) - e_initial, psi.front(), psi.back(), norm, r_mean};
}

// <V>, <E_rot> and <r> of the initial state, normalised with its norm.
InitialStateProperties initial_state_properties(const RovibrationalStates& states, int v, int j,
                                                const UniformGrid& r_grid, double norm,
                                                std::span<const double> asymptotic_potential,
                                                double reduced_mass) {
  const auto psi = states.wave_function(v, j);
  double erot = 0.0;
  double r_in = 0.0;
  double vpot = 0.0;
  for (int ir = 0; ir < r_grid.size(); ++ir) {
    const double y = psi[ir];
    const double r = r_grid[ir];
    r_in = r_in + r * (y * y);
    erot = erot + 1.0 / (2.0 * reduced_mass * (r * r)) * (j * (j + 1)) * (y * y);
    vpot = vpot + asymptotic_potential[ir] * (y * y);
  }
  const double dr = r_grid.spacing();
  InitialStateProperties initial;
  initial.energy = states.energy(v, j);
  initial.rotational_energy = erot * dr / norm;
  initial.potential_energy = vpot * dr / norm;
  initial.r_mean = r_in * dr / norm;
  initial.vibrational_energy = initial.energy - initial.potential_energy - initial.rotational_energy;
  initial.rovibrational_kinetic_energy = initial.vibrational_energy + initial.rotational_energy;
  return initial;
}

}  // namespace

GasPhase compute_gas_phase(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                           std::span<const double> asymptotic_potential, const Fft& fft,
                           const SymmetricEigensolver& solver) {
  const InputParameters& p = config.input;
  const auto& units = units::legacy;
  const int nr = p.nr;
  const double reduced_mass = config.derived.reduced_mass;

  // The FGH uses its own grid on [rmin, rmin + nr*dr). Fortran: InitFourierGridHamiltonian
  const FourierGridHamiltonian fgh(nr, p.r_start, p.r_start + nr * p.dr, reduced_mass, fft);
  if (fgh.size() != nr || asymptotic_potential.size() != static_cast<std::size_t>(nr))
    throw std::invalid_argument("compute_gas_phase: nR must be a valid FFT size");

  GasPhase gas(solve_rovibrational_states(fgh, p.j_max, asymptotic_potential, reduced_mass, solver));
  const RovibrationalStates& states = gas.states;
  const double e00 = states.energy(0, 0);
  const double e_initial = states.energy(p.v_in, p.j_in);

  // Open channels: maximum wave-packet energy plus the internal and parallel energy of the
  // initial state, relative to the ground state.
  const double available_energy =
      units.ev_to_hartree(p.e_max_wave_packet) + (e_initial - e00) + units.ev_to_hartree(p.e_parallel);
  gas.v_max_open = count_open_vibrational_levels(states, available_energy, gas.warnings);
  gas.j_max_open = highest_open_rotational_level(states, angular.delta_j(), available_energy, gas.warnings);

  // VAna < 1 (or absent) and absent JAna mean: determined by the program.
  gas.n_v_analysis = (p.v_analysis && *p.v_analysis >= 1) ? *p.v_analysis : gas.v_max_open;
  gas.j_max_analysis = (p.j_analysis && *p.j_analysis >= 0) ? *p.j_analysis : gas.j_max_open;

  for (int v = 0; v < gas.n_v_analysis; ++v) {
    for (int j = 0; j <= p.j_max; ++j) {
      gas.summary.push_back(summarise_state(states, v, j, grids.r, e_initial));
      if (v == p.v_in && j == p.j_in)
        gas.initial = initial_state_properties(states, v, j, grids.r, gas.summary.back().norm,
                                               asymptotic_potential, reduced_mass);
    }
  }

  gas.initial.energy = e_initial;

  // Warnings if not all channels within the 95% confidence energy are analysed (they do not
  // change the analysis).
  const double e_check = e_initial + units.ev_to_hartree(config.derived.confidence_e_max);
  const int next_j = gas.j_max_analysis + angular.delta_j();
  const int j_limit = p.homonuclear ? p.j_max - 1 : p.j_max;
  if (gas.j_max_analysis < j_limit && e_check > states.energy(0, next_j))
    gas.warnings.emplace_back("Not all open j channels analysed!");
  if (gas.n_v_analysis < nr && e_check > states.energy(gas.n_v_analysis, 0))
    gas.warnings.emplace_back("Not all open v channels analysed!");

  for (int j = 0; j <= p.j_max; ++j) gas.rotational_energies.push_back(states.energy(0, j) - e00);
  for (int v = 0; v < gas.n_v_analysis; ++v) gas.vibrational_energies.push_back(states.energy(v, 0) - e00);
  return gas;
}

}  // namespace spodvr
