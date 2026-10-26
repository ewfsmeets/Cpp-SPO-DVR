#include "gasphase/gas_phase.hpp"

#include <stdexcept>

#include "core/machine_epsilon.hpp"
#include "core/unit_conversion.hpp"
#include "gasphase/fourier_grid_hamiltonian.hpp"

namespace spodvr {

GasPhase compute_gas_phase(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                           std::span<const double> asymptotic_potential, const Fft& fft,
                           const SymmetricEigensolver& solver) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const auto& units = units::legacy;
  const int nr = p.nr;
  const double reduced_mass = d.reduced_mass;

  const FourierGridHamiltonian fgh(nr, p.r_start, p.r_start + static_cast<double>(nr) * p.dr,
                                   reduced_mass, fft);
  if (fgh.size() != nr || asymptotic_potential.size() != static_cast<std::size_t>(nr))
    throw std::invalid_argument("compute_gas_phase: nR must be a valid FFT size");

  GasPhase gas{RovibrationalStates(nr, p.j_max), 0, 0, 0, 0, {}, 0, 0, 0, 0, 0, 0, {}, {}, {}};
  RovibrationalStates& states = gas.states;

  // Fortran: SetupAsympWF loop over Jlevel with CopyAsympPotentialToFGHType
  std::vector<double> potential(static_cast<std::size_t>(nr));
  for (int j = 0; j <= p.j_max; ++j) {
    for (int i = 0; i < nr; ++i) {
      const auto ui = static_cast<std::size_t>(i);
      if (j != 0) {
        double r = fgh.grid()[ui];
        if (r == 0.0) r = r + fortran_machine_epsilon();  // prevent division by zero
        potential[ui] = asymptotic_potential[ui] +
                        static_cast<double>(j * (j + 1)) / 2.0 / reduced_mass / (r * r);
      } else {
        potential[ui] = asymptotic_potential[ui];
      }
    }
    const EigenDecomposition eigen = fgh.eigenstates(potential, solver);
    for (int v = 0; v < nr; ++v) {
      states.energy(v, j) = eigen.values[static_cast<std::size_t>(v)];
      const auto source = eigen.vector(v);
      const auto target = states.wave_function(v, j);
      std::copy(source.begin(), source.end(), target.begin());
    }
  }

  // Highest open vibrational and rotational states. Fortran: SetupAsympWF (E. Pijper part)
  const double e00 = states.energy(0, 0);
  const double e_initial = states.energy(p.v_in, p.j_in);
  const double e_parallel = units.ev_to_hartree(p.e_parallel);
  const double total_energy = units.ev_to_hartree(p.e_max_wave_packet) + (e_initial - e00) + e_parallel;

  int iv = 1;  // 1-based, as in the Fortran loop
  while (true) {
    if (states.energy(iv - 1, 0) - e00 > total_energy) break;
    iv = iv + 1;
    if (iv == nr) {
      gas.warnings.emplace_back("Vibrational basis too small to describe all open channels!");
      break;
    }
  }
  gas.v_max_open = iv - 1;
  // VAna < 1 (or absent) means: determined by the program.
  gas.n_v_analysis = (p.v_analysis && *p.v_analysis >= 1) ? *p.v_analysis : gas.v_max_open;

  const int delta_j = angular.delta_j();
  int j = 0;
  while (true) {
    if (states.energy(0, j) - e00 > total_energy) break;
    j = j + delta_j;
    if (j > p.j_max) {
      gas.warnings.emplace_back("Rotational basis too small to describe all open channels!");
      break;
    }
  }
  gas.j_max_open = j - delta_j;
  gas.j_max_analysis = (p.j_analysis && *p.j_analysis >= 0) ? *p.j_analysis : gas.j_max_open;

  // Table of the analysed states and the initial-state expectation values.
  const double dr = p.dr;
  for (int v = 0; v < gas.n_v_analysis; ++v) {
    for (int jj = 0; jj <= p.j_max; ++jj) {
      const auto psi = states.wave_function(v, jj);
      double norm = 0.0;
      double r_mean = 0.0;
      for (int ir = 0; ir < nr; ++ir) {
        const double y = psi[static_cast<std::size_t>(ir)];
        norm = norm + y * y;
        r_mean = r_mean + grids.r[ir] * y * y;
      }
      norm = norm * dr;
      r_mean = r_mean * dr;
      gas.summary.push_back({v, jj, states.energy(v, jj), states.energy(v, jj) - e_initial, psi.front(),
                             psi.back(), norm, r_mean});

      if (jj == p.j_in && v == p.v_in) {
        double erot = 0.0;
        double r_in = 0.0;
        double vpot = 0.0;
        for (int ir = 0; ir < nr; ++ir) {
          const double y = psi[static_cast<std::size_t>(ir)];
          const double r = grids.r[ir];
          r_in = r_in + r * (y * y);
          erot = erot + 1.0 / (2.0 * reduced_mass * (r * r)) * static_cast<double>(jj * (jj + 1)) * (y * y);
          vpot = vpot + asymptotic_potential[static_cast<std::size_t>(ir)] * (y * y);
        }
        gas.initial_rotational_energy = erot * dr / norm;
        gas.initial_potential_energy = vpot * dr / norm;
        gas.initial_r_mean = r_in * dr / norm;
        gas.initial_vibrational_energy =
            states.energy(v, jj) - gas.initial_potential_energy - gas.initial_rotational_energy;
        gas.initial_rovibrational_kinetic_energy =
            gas.initial_vibrational_energy + gas.initial_rotational_energy;
      }
    }
  }

  // Warnings if not all open channels are analysed (they do not change the analysis).
  const double e_check = e_initial + units.ev_to_hartree(d.confidence_e_max);
  if (p.homonuclear) {
    if (gas.j_max_analysis < p.j_max - 1 && e_check > states.energy(0, gas.j_max_analysis + 2))
      gas.warnings.emplace_back("Not all open j channels analysed!");
  } else {
    if (gas.j_max_analysis < p.j_max && e_check > states.energy(0, gas.j_max_analysis + 1))
      gas.warnings.emplace_back("Not all open j channels analysed!");
  }
  if (gas.n_v_analysis < nr && e_check > states.energy(gas.n_v_analysis, 0))
    gas.warnings.emplace_back("Not all open v channels analysed!");

  gas.initial_energy = e_initial;
  for (int jj = 0; jj <= p.j_max; ++jj) gas.rotational_energies.push_back(states.energy(0, jj) - e00);
  for (int v = 0; v < gas.n_v_analysis; ++v) gas.vibrational_energies.push_back(states.energy(v, 0) - e00);
  return gas;
}

}  // namespace spodvr
