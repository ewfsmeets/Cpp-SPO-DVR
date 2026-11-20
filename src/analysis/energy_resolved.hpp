#pragma once

#include <complex>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "gasphase/gas_phase.hpp"
#include "input/config.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// The energy-resolved partial sums are taken at the first analysis step at or after this fraction of tEnd.
constexpr double partial_sums_fraction = 0.8;

// The energy-resolved sums of a propagation (EnergyResolvedSites, EnergyResolvedEnergies; C++ only, no Fortran
// counterpart): at each chosen (x, y) grid point (site) and collision energy E,
//   s(Z, r, irot; E) = sum over the analysis steps of psi_scat(Z, r, irot; t) exp(i E_tot t) dt_ana,
// and per energy the same sum of the specular wave chi(Z, r), sigma(Z, r; E). The stationary scattering state at
// energy E follows from them with the S matrix's prefactor (tools/plots, compute/energy_resolved.py;
// docs/output_files.md "energy_resolved/").
struct EnergyResolvedSums {
  std::vector<std::string> site_names;  // EnergyResolvedSites, in input order
  std::vector<int> site_ixy;            // their (x, y) grid points, ixy = iy + ny * ix
  std::vector<double> collision_ev;     // EnergyResolvedEnergies
  std::vector<double> total_hartree;    // E_tot = E_coll + E(vin, jin), as the S matrix's energies
  int n_z = 0;                          // scattering grid
  int n_r = 0;
  int n_rot = 0;
  int n_z_specular = 0;
  double analysis_time_step = 0.0;  // nanspo * dt (a.u.)
  double last_time = 0.0;           // the time of the last analysis step summed (a.u.); t_start before any
  int steps_summed = 0;             // analysis steps summed (all segments of a restarted run)
  std::vector<std::complex<double>> sums;           // [site, energy, Z, r, rot], rot fastest
  std::vector<std::complex<double>> specular_sums;  // [energy, Z_specular, r], r fastest

  [[nodiscard]] std::size_t site_block() const { return std::size_t(n_z) * n_r * n_rot; }
  [[nodiscard]] std::size_t specular_block() const { return std::size_t(n_z_specular) * n_r; }
};

// The complex values the sums of a run hold in host memory: every (site, energy) pair's nZ nR nrot sums, every
// energy's specular sums (nZ_specular nR) and the column buffer of one site (the memory estimate's item).
[[nodiscard]] std::size_t energy_resolved_memory_values(const InputParameters& p, int n_rot, int n_z_specular);

// The sums before the first analysis step of a new run (all zero; last_time = t_start).
[[nodiscard]] EnergyResolvedSums empty_energy_resolved_sums(const Config& config, const GasPhase& gas, int n_rot,
                                                            int n_z_specular, double t_start);

// Adds every analysis step to the sums. A new run starts from empty_energy_resolved_sums, a restart from the sums
// of the earlier segment (read from its energy_resolved/).
class EnergyResolvedAccumulator {
 public:
  // read_column(ixy, out): the scattering wave at grid point ixy, layout (izr, irot) with izr fastest
  // (PropagationBackend::site_column).
  using ColumnReader = std::function<void(int ixy, std::span<std::complex<double>> out)>;

  explicit EnergyResolvedAccumulator(EnergyResolvedSums start);

  // Adds the analysis step that ends at `time` (as SMatrixAccumulator::add: sum + value * exp(i E_tot t) * dt_ana).
  void add(double time, const ColumnReader& read_column, const SpecularWave& specular);

  [[nodiscard]] const EnergyResolvedSums& sums() const { return sums_; }

 private:
  EnergyResolvedSums sums_;
  std::vector<std::complex<double>> column_;  // one site's (izr, irot): the only buffer the sums need
};

}  // namespace spodvr
