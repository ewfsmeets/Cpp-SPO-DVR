#pragma once
// energy_resolved/ (EnergyResolvedSites, EnergyResolvedEnergies; C++ only): the energy-resolved sums of the run, an
// array directory (io/series_writer.hpp ArraysWriter) written in two steps, the partial sums at the first analysis
// step at or after 0.8 tEnd and the final sums at the end of the run; their reader for a restart; the projected size
// and cost (docs/output_files.md, "energy_resolved/").

#include <complex>
#include <filesystem>
#include <optional>
#include <ostream>
#include <span>

#include "analysis/energy_resolved.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

class EnergyResolvedFiles {
 public:
  // Creates energy_resolved/ in output_dir with its manifest (the sites and their s2s factors, the energies, the
  // grids, the packet). `layout`: the sums the run starts from (their shape and keys); s2s: the specular-to-
  // scattering factor of every (x, y) grid point (PropagationTables::specular_to_scattering).
  EnergyResolvedFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
                      const EnergyResolvedSums& layout, std::span<const std::complex<double>> s2s);

  // The partial sums and the time of their last analysis step (once).
  void write_partial(const EnergyResolvedSums& sums);
  // The final sums, the time of their last analysis step and the steps summed (and the partial sums, the same,
  // when none were written: a segment without analysis steps); closes the directory.
  void write_final(const EnergyResolvedSums& sums);

 private:
  ArraysWriter arrays_;
  bool partial_written_ = false;
};

// Restart, before the set-up: RestartDir's sums must have this input's sites (names and grid points) and energies
// (InputError naming EnergyResolvedSites otherwise; CheckpointError when the manifest cannot be read).
void check_energy_resolved_restart(const std::filesystem::path& restart_dir, const InputParameters& p);

// Restart: the sums of the earlier segment in restart_dir/energy_resolved/, which must hold the same sites and
// energies (InputError) and grids as `expected` (the empty sums of this run), end at tStart and hold as many analysis
// steps as the earlier BK records read for the restart (CheckpointError otherwise: sums of another run).
[[nodiscard]] EnergyResolvedSums read_energy_resolved_sums(const std::filesystem::path& restart_dir,
                                                           const Config& config, const EnergyResolvedSums& expected,
                                                           int earlier_steps);

// The measured cost at the end of the run: the pairs summed and the wall-clock time of the sums (column copies and
// accumulation).
void print_energy_resolved_cost(std::ostream& out, std::size_t pairs, double seconds, double propagation_seconds);

// The size of energy_resolved/ and its projected cost per analysis step, or nothing without the keys. Above
// energy_resolved_warning_fraction of an analysis step the line warns.
[[nodiscard]] std::optional<PlotDataSize> energy_resolved_data_size(const Config& config, const SpatialGrids& grids,
                                                                    int n_rot);

// Above this share of an analysis step the projected cost of the sums warns (the PLOT DATA line).
inline constexpr double energy_resolved_warning_fraction = 0.10;

}  // namespace spodvr
