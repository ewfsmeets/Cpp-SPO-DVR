#pragma once

#include <complex>
#include <cstdint>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/config.hpp"

namespace spodvr {

// Energies at which the S-matrix is computed on the fly. NrOfEs = -95 (nen95) chooses the
// spectral resolution pi / tEnd of the propagation instead of NrOfEs equidistant energies.
// Fortran: Analyzor.f90 SetupAnalysis (enrgi, denrgy, nenrgy, ensctr, enscev, enscau)
struct AnalysisEnergies {
  std::vector<double> collision;     // translational energy (eV). Fortran: ensctr
  std::vector<double> total;         // collision + E(vin, jin) (eV). Fortran: enscev
  std::vector<double> total_hartree; // Fortran: enscau

  [[nodiscard]] int size() const { return static_cast<int>(collision.size()); }
};

[[nodiscard]] AnalysisEnergies make_analysis_energies(const Config& config, const GasPhase& gas);

// S-matrix element, probability and phase of every analysed state at one energy, plus the
// sums and averages of Fortran smatrix_combine. Layouts: (is, v) = is + nstata * v,
// (irotana, v) = irotana + nana * v, (order, irotana, v) = order + (mdifa+1) * (irotana + nana * v).
// Energies in eV.
struct EnergyResult {
  double total_energy;      // Fortran: enscev (senrg)
  double collision_energy;  // Fortran: ensctr (etrans)
  double z_energy;          // collision energy minus the initial parallel energy. Fortran: etransZ
  double total_momentum;    // sqrt(2 M E_coll), a.u. Fortran: rkt

  std::vector<std::complex<double>> s_matrix;  // (is, v); 0 for closed channels
  std::vector<double> probability;              // |S|^2 (is, v)
  std::vector<double> phase;                    // arg S (0 if Re S = 0, as in Fortran)
  std::vector<std::uint8_t> open;               // channel energetically open (is, v)

  std::vector<double> p_vibration;   // summed over rotation and diffraction (v). Fortran: sumvb
  std::vector<double> p_rotation;    // summed over diffraction (irotana, v). Fortran: sumrot
  std::vector<double> p_order;       // summed per diffraction order (order, irotana, v). Fortran: sumrd
  std::vector<double> e_rotation_v;  // <E_rot> of the molecules scattered into v. Fortran: averov
  std::vector<double> e_diffraction_v;  // <E_dif> per v. Fortran: avediv
  double p_total;          // sum of all probabilities. Fortran: sump
  double p_reaction;       // 1 - p_total
  double e_vibration;      // <E_vib> over all scattered molecules. Fortran: avevib
  double e_rotation;       // Fortran: averot
  double e_diffraction;    // Fortran: avedif
};

struct ScatteringResults {
  std::vector<EnergyResult> energies;
};

// Accumulates the energy-resolved S-matrix from the BK coefficients during the propagation,
// S(is, v; E) = sum_t c(is, v; t) exp(i E t) nanspo dt, and turns it into S-matrix elements and
// probabilities at the end. Fortran: Analyzor.f90 BalintKurtiCoeffs (accumulation),
// StateToStateProbs, smatrix_calc and smatrix_combine (nread = 0, 1)
class SMatrixAccumulator {
 public:
  SMatrixAccumulator(const Config& config, const AnalysisBasis& basis, const GasPhase& gas);

  [[nodiscard]] const AnalysisEnergies& energies() const { return energies_; }

  // Adds the coefficients of the analysis step that ends at `time`.
  void add(const BkCoefficients& bk, double time);

  [[nodiscard]] ScatteringResults results() const;

 private:
  [[nodiscard]] EnergyResult result_at(int ien) const;

  const Config& config_;
  const AnalysisBasis& basis_;
  const GasPhase& gas_;
  AnalysisEnergies energies_;
  double analysis_time_step_;  // nanspo * dt. Fortran: dtmul
  std::vector<std::complex<double>> sums_;  // (is, v, ien). Fortran: smtrx
};

}  // namespace spodvr
