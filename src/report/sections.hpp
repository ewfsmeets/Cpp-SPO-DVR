#pragma once

#include <filesystem>
#include <ostream>

#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "report/format.hpp"

namespace spodvr {

// Presentation of the setup results. Numerical modules return data; only these functions print.
void print_input_summary(std::ostream& out, const Style& style, const Config& config);
void print_grids(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                 const AngularGrid& angular);
void print_asymptotic_potential(std::ostream& out, const Style& style, const Config& config,
                                const SpatialGrids& grids, const PotentialGrid& potential);
void print_gas_phase(std::ostream& out, const Style& style, const Config& config, const GasPhase& gas);

// Full-precision data files (regression gates). Fortran equivalents: the vasy table in stdout and
// debug_asyp.dat.
void write_asymptotic_potential_file(const std::filesystem::path& path, const SpatialGrids& grids,
                                     const PotentialGrid& potential);
void write_gas_phase_file(const std::filesystem::path& path, const SpatialGrids& grids, const GasPhase& gas);

}  // namespace spodvr
