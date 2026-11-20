// Input checks of the plot-data keys (C++ only): writeDensities, writePotential, PlotSites. Each check names its
// key (exit 2).
#include "input/validate_outputs.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <set>

#include "input/validation_helpers.hpp"

namespace spodvr {

using validation::check;

namespace {

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

void check_sites(const InputParameters& p) {
  if (p.plot_sites.empty()) return;
  check(!p.write_potential, "PlotSites", "no plot key uses the sites: set writePotential");
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

}  // namespace

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
  check_sites(p);
}

}  // namespace spodvr
