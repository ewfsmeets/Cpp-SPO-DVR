#pragma once

#include <ostream>
#include <string>
#include <string_view>

#include "analysis/analysis_basis.hpp"
#include "analysis/non_finite.hpp"
#include "analysis/s_matrix.hpp"
#include "input/config.hpp"
#include "report/format.hpp"

namespace spodvr {

// The analysis basis, energy grid and analysis surface. Fortran: the stdout of RotDifBasis.f90
// and Analyzor.f90 SetupAnalysis
void print_analysis_setup(std::ostream& out, const Style& style, const Config& config, const AnalysisBasis& basis,
                          const AnalysisEnergies& energies);

// Compact result tables (all details are in the data files): one summary row per energy, the
// rovibrationally resolved probabilities of the open (v, j) states, and the diffraction orders
// of the initial (v, j). Replaces the per-energy blocks of Fortran smatrix_combine.
void print_scattering_results(std::ostream& out, const Style& style, const Config& config, const AnalysisBasis& basis,
                              const ScatteringResults& results);

// The note above the result tables when they show only every k-th of more than 50 energies,
// e.g. "One energy in 3 is shown; all are in the data files."; empty when every energy is shown.
[[nodiscard]] std::string energy_thinning_note(int n_energies);

// One line naming the non-finite number `value` of `results` (find_non_finite): the result set (`set`, e.g. "the
// on-the-fly analysis"), the energy, the quantity with its state labels and the number itself, e.g. "the on-the-fly
// analysis, energy 3 of 41 (E_coll = 0.220000 eV): S of (v, j, mJ, n, m) = (0, 0, 0, 0, 0) is (nan, nan)".
[[nodiscard]] std::string describe_non_finite(std::string_view set, const NonFiniteValue& value,
                                              const AnalysisBasis& basis, const ScatteringResults& results);

}  // namespace spodvr
