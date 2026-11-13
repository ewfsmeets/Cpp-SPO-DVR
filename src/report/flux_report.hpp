#pragma once

#include <ostream>

#include "analysis/flux_analysis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/memory_estimate.hpp"
#include "report/format.hpp"

namespace spodvr {

// The FLUX ANALYSIS set-up block: the flux surface r_c and its index, the energies used (with the
// warnings of `energies`, flux Q3), initial_state_ratio = flux_initial_state_ratio (flux Q7; a warning
// above flux_initial_state_ratio_limit), and from the memory rule the slot count, the memory per
// energy and the most energies that fit. Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (its printout)
void print_flux_analysis_setup(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                               const FluxEnergies& energies, double initial_state_ratio, const FluxMemoryCheck& memory);

// The reaction probability per flux energy; above max_printed_energies energies every k-th
// (energy_stride), all of them in flux_reaction.tsv.
void print_flux_results(std::ostream& out, const Style& style, const FluxResults& results);

// The memory rule of the flux analysis, below the memory estimate: the host memory of the run without
// it (set-up, time loop; a GPU run's own), its coefficients and slots, the physical memory and the most
// energies that fit: the numbers of require_flux_memory's error.
void print_flux_memory(std::ostream& out, const Style& style, const FluxMemoryCheck& memory);

}  // namespace spodvr
