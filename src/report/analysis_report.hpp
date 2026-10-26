#pragma once

#include <ostream>

#include "analysis/analysis_basis.hpp"
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

}  // namespace spodvr
