// Input checks of the plot-data keys (C++ only): writeDensities, writePotential, PlotSites, dtFrames, FrameExtras,
// NormBudget, EnergyResolvedSites, EnergyResolvedEnergies. Each check names its key (exit 2).
#include "input/validate_outputs.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <set>

#include "core/real_compare.hpp"
#include "input/validation_helpers.hpp"

namespace spodvr {

using validation::check;

namespace {

constexpr double frame_grid_tolerance = 1e-9;  // tStart / dtFrames and tEnd / dtFrames: whole to this

constexpr double on_grid_tolerance = 1e-9;  // |fraction * n - nearest integer|: a site lies on a grid point

// The grid index nearest to a fractional coordinate (0 <= fraction < 1) on n points.
int nearest_index(double fraction, int n) { return static_cast<int>(std::lround(fraction * n)) % n; }

// One coordinate of a site: in [0, 1) and on one of the n grid points (the error names the nearest one).
void check_site_coordinate(const PlotSite& site, const std::string& axis, double fraction, int n,
                           const std::string& n_key) {
  check(
      !(fraction >= 0.0 && fraction < 1.0), "PlotSites",
      std::format("site '{}': {} = {} is not in [0, 1) (a fraction of the lattice vector)", site.name, axis, fraction));
  const double scaled = fraction * n;
  const double nearest = std::round(scaled);
  const int index = nearest_index(fraction, n);
  check(!(std::abs(scaled - nearest) <= on_grid_tolerance), "PlotSites",
        std::format("site '{}': {} = {} is not on the grid ({} = {}); the nearest point is {} = {} (index {})",
                    site.name, axis, fraction, n_key, n, axis, static_cast<double>(index) / n, index));
  check(nearest == n, "PlotSites",
        std::format("site '{}': {} = {} rounds to 1, the next cell's {} = 0: give 0", site.name, axis, fraction, axis));
}

// Whether t / dt is a whole number to `tolerance` (in units of dt, so that it holds for any size of t and dt).
bool whole_multiple(double t, double dt, double tolerance) {
  const double ratio = t / dt;
  return std::abs(ratio - validation::nint(ratio)) <= tolerance;
}

// dtFrames (C++ only checks): frames line up with the analysis steps (to an absolute 1e-12, as the Fortran's dtAna
// checks) and, across restarts, with tStart and tEnd (relative to the frame interval).
void check_frames(const InputParameters& p, bool propagates) {
  using fortran_compare::rlt;
  using fortran_compare::rne;
  using validation::nint;
  if (!p.dt_frames) return;
  const double dt = *p.dt_frames;
  check(p.write_2d_pes, "dtFrames", "write2DPES stops the run after the cut of the PES, before the propagation");
  check(!propagates, "dtFrames", "only a propagating run (a new run or a restart) writes frames");
  check(!(dt > 0.0), "dtFrames", "dtFrames <= 0.");
  check(rlt(dt, p.dt_analysis), "dtFrames", "dtFrames < dtAna.");
  // The program holds the analyses per frame in an int (analysis D4); the rounding of the next check also needs a
  // ratio within the range of a long.
  check(dt / p.dt_analysis > std::numeric_limits<int>::max(), "dtFrames",
        std::format("dtFrames / dtAna = {} is above {}, the largest number of analyses per frame the program can hold.",
                    dt / p.dt_analysis, std::numeric_limits<int>::max()));
  check(rne(nint(dt / p.dt_analysis) * p.dt_analysis, dt), "dtFrames", "dtFrames must be a whole multiple of dtAna.");
  check(!whole_multiple(p.t_start, dt, frame_grid_tolerance), "dtFrames",
        "tStart must be a whole multiple of dtFrames (the frames of a restart continue those of the run before)");
  check(!whole_multiple(p.t_end, dt, frame_grid_tolerance), "dtFrames", "tEnd must be a whole multiple of dtFrames.");
}

void check_sites(const InputParameters& p) {
  if (p.plot_sites.empty()) return;
  check(!p.write_potential && !p.dt_frames && p.energy_resolved_sites.empty(), "PlotSites",
        "no plot key uses the sites: set writePotential, dtFrames or EnergyResolvedSites");
  check(static_cast<int>(p.plot_sites.size()) > max_plot_sites, "PlotSites",
        std::format("{} sites: at most {}", p.plot_sites.size(), max_plot_sites));
  std::set<std::string> names;
  for (const PlotSite& site : p.plot_sites) {
    check(site.name.empty(), "PlotSites", "a site without a name");
    check(!names.insert(site.name).second, "PlotSites", std::format("two sites are named '{}'", site.name));
    check_site_coordinate(site, "x", site.x, p.nx, "nX");
    check_site_coordinate(site, "y", site.y, p.ny, "nY");
  }
}

// EnergyResolvedSites and EnergyResolvedEnergies: both or neither, in a propagating run; sites named in PlotSites
// (once each), collision energies inside the packet's range EMinWP..EMaxWP (+ EParallel), at most
// max_energy_resolved_pairs (site, energy) pairs. Their cost and memory are projected at set-up (the PLOT DATA line
// warns above 10 % of an analysis step; the memory estimate holds the sums).
void check_energy_resolved(const InputParameters& p, bool propagates) {
  const bool sites = !p.energy_resolved_sites.empty();
  const bool energies = !p.energy_resolved_energies.empty();
  if (!sites && !energies) return;
  check(!energies, "EnergyResolvedSites", "set EnergyResolvedEnergies too (the sums need sites and energies)");
  check(!sites, "EnergyResolvedEnergies", "set EnergyResolvedSites too (the sums need sites and energies)");
  check(p.write_2d_pes, "EnergyResolvedSites",
        "write2DPES stops the run after the cut of the PES, before the propagation");
  check(!propagates, "EnergyResolvedSites", "only a propagating run (a new run or a restart) sums the wave packet");
  std::set<std::string> seen;
  for (const std::string& name : p.energy_resolved_sites) {
    const bool known = std::ranges::any_of(p.plot_sites, [&](const PlotSite& site) { return site.name == name; });
    check(!known, "EnergyResolvedSites", std::format("'{}' is not a site of PlotSites", name));
    check(!seen.insert(name).second, "EnergyResolvedSites", std::format("'{}' is listed twice", name));
  }
  const double low = p.e_min_wave_packet + p.e_parallel;
  const double high = p.e_max_wave_packet + p.e_parallel;
  for (const double energy : p.energy_resolved_energies)
    check(fortran_compare::rlt(energy, low) || fortran_compare::rgt(energy, high), "EnergyResolvedEnergies",
          std::format("{:g} eV lies outside the packet's collision energies {:g} .. {:g} eV (EMinWP, EMaxWP + "
                      "EParallel)",
                      energy, low, high));
  const std::size_t pairs = p.energy_resolved_sites.size() * p.energy_resolved_energies.size();
  check(pairs > max_energy_resolved_pairs, "EnergyResolvedSites",
        std::format("{} (site, energy) pairs: at most {}", pairs, max_energy_resolved_pairs));
}

}  // namespace

const std::vector<std::string>& frame_extra_names() {
  static const std::vector<std::string> names = {"angular", "diffraction", "vibration", "rotation_moments", "momentum"};
  return names;
}

const std::vector<std::string>& plot_site_fields() {
  static const std::vector<std::string> fields = {"name", "x", "y"};
  return fields;
}

bool is_plot_site_field(const std::string& field) {
  return std::ranges::find(plot_site_fields(), field) != plot_site_fields().end();
}

int plot_site_ixy(const PlotSite& site, int nx, int ny) {
  return nearest_index(site.y, ny) + ny * nearest_index(site.x, nx);
}

void validate_output_keys(const Config& config) {
  const InputParameters& p = config.input;
  const int mode = config.derived.run_mode;
  const bool propagates = mode == 0 || mode == 1;
  // write2DPES stops the run right after the PES is set up (as the Fortran does): no plot data follow.
  check(p.write_2d_pes && p.write_densities, "writeDensities",
        "write2DPES stops the run after the cut of the PES, before the propagation");
  check(p.write_2d_pes && p.write_potential, "writePotential",
        "write2DPES stops the run after the cut of the PES, before the 6D potential");
  check(p.write_densities && !propagates, "writeDensities",
        "only a propagating run (a new run or a restart) writes densities");
  check(p.write_potential && mode == 2, "writePotential",
        "an analysis-only run (onlyAnalysis) computes no 6D potential to write");
  check_frames(p, propagates);
  if (!p.frame_extras.empty()) {  // FrameExtras: optional datasets of the frames
    check(!p.dt_frames, "FrameExtras", "the extras are datasets of frames/: set dtFrames");
    std::set<std::string> seen;
    for (const std::string& name : p.frame_extras) {
      check(std::ranges::find(frame_extra_names(), name) == frame_extra_names().end(), "FrameExtras",
            std::format("unknown extra '{}': known are angular, diffraction, vibration, rotation_moments, momentum",
                        name));
      check(!seen.insert(name).second, "FrameExtras", std::format("'{}' is listed twice", name));
    }
  }
  check_sites(p);
  check_energy_resolved(p, propagates);
  if (p.norm_budget) {  // NormBudget: one record per analysis window of the propagation
    check(p.write_2d_pes, "NormBudget", "write2DPES stops the run after the cut of the PES, before the propagation");
    check(!propagates, "NormBudget", "only a propagating run (a new run or a restart) has a norm budget");
  }
}

}  // namespace spodvr
