#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <ostream>
#include <span>
#include <string_view>

#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/memory_estimate.hpp"
#include "report/format.hpp"

namespace spodvr {

// Presentation of the setup results. Numerical modules return data; only these functions print.
void print_input_summary(std::ostream& out, const Style& style, const Config& config);
void print_grids(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                 const AngularGrid& angular);
// The surface symmetry used for the potential (UseSymmetry): point group, usable and dropped
// operations, reduction, and the result of the PES self-check.
void print_surface_symmetry(std::ostream& out, const Style& style, const GridSymmetry& symmetry,
                            const SymmetryCheck& check);
// The expected memory, with a warning if the peak exceeds the physical memory (if known).
void print_memory_estimate(std::ostream& out, const Style& style, const MemoryEstimate& estimate,
                           std::optional<std::size_t> physical_bytes, std::string_view title = "MEMORY ESTIMATE");
void print_asymptotic_potential(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                                std::span<const double> asymptotic);
void print_gas_phase(std::ostream& out, const Style& style, const Config& config, const GasPhase& gas);

// Full-precision data files (regression gates). Fortran equivalents: the vasy table in stdout and
// debug_asyp.dat.
void write_asymptotic_potential_file(const std::filesystem::path& path, const SpatialGrids& grids,
                                     std::span<const double> asymptotic);
void write_gas_phase_file(const std::filesystem::path& path, const SpatialGrids& grids, const GasPhase& gas);

}  // namespace spodvr
