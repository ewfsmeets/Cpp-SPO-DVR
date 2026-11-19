#pragma once

#include <filesystem>
#include <functional>
#include <optional>

#include "core/program_info.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"

namespace spodvr {

// Writes run_info.json: what a plot or a reader needs to interpret the other files of a run, in
// one place: program, input, times, molecule, initial state and wave packet, grids, rotational
// basis, lattice, analysis surface and energies, absorbers, potential settings and the unit
// factors. Atomic units (bohr, hartree, a.u. of time, electron masses) unless a field name ends
// in a unit (_ev, _amu, _deg, _rad). Without a gas phase (`gas` empty; the write2DPES stop) its
// fields are null.
// Format: docs/output_files.md, section "Run metadata".
void write_run_info(const std::filesystem::path& path, const Config& config, const SpatialGrids& grids,
                    const AngularGrid& angular, const ProgramInfo& program, const std::filesystem::path& input_file,
                    std::optional<std::reference_wrapper<const GasPhase>> gas);

}  // namespace spodvr
