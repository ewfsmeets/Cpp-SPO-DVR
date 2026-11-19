// potential/ (writePotential): V of the run on its grid at set-up, an array directory (io/series_writer.hpp
// ArraysWriter), read through PotentialGrid::operator() like the propagation's gather. C++ only.
#include "report/potential_files.hpp"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "core/unit_conversion.hpp"
#include "input/validate_outputs.hpp"
#include "io/series_writer.hpp"

namespace spodvr {

namespace {

// The datasets of potential/, in the order write_potential_files writes them (v_sites only with sites).
std::vector<DatasetSpec> potential_datasets(const Config& config, const SpatialGrids& grids,
                                            const AngularGrid& angular) {
  const InputParameters& p = config.input;
  const int nzpot = grids.nz_potential;
  std::vector<DatasetSpec> datasets = {
      {"v_min.bin",
       "<f8",
       {nzpot, p.nr, p.nx, p.ny},
       {"Z", "r", "x", "y"},
       "hartree",
       "the minimum of V over the DVR orientations (theta, phi nodes) at each grid point of Z < zasymp, as the "
       "propagation uses V (switching and the VMin/VMax limits included)"},
      {"v_min_node.bin",
       "<i4",
       {nzpot, p.nr, p.nx, p.ny, 2},
       {"Z", "r", "x", "y", "node"},
       "",
       "the orientation of that minimum: (itheta, iphi), indices into theta_nodes_rad and phi_nodes_rad"},
  };
  const int n_sites = static_cast<int>(p.plot_sites.size());
  if (n_sites > 0)
    datasets.push_back({"v_sites.bin",
                        "<f8",
                        {n_sites, nzpot, p.nr, angular.n_theta(), angular.n_phi()},
                        {"site", "Z", "r", "theta", "phi"},
                        "hartree",
                        "V at every DVR orientation at the (x, y) grid point of each site of PlotSites"});
  return datasets;
}

nlohmann::ordered_json potential_manifest_fields(const Config& config, const SpatialGrids& grids,
                                                 const AngularGrid& angular) {
  const InputParameters& p = config.input;
  nlohmann::ordered_json sites = nlohmann::ordered_json::array();
  for (const PlotSite& site : p.plot_sites) {
    const int ixy = plot_site_ixy(site, p.nx, p.ny);
    sites.push_back(
        {{"name", site.name}, {"x", site.x}, {"y", site.y}, {"ix", ixy / p.ny}, {"iy", ixy % p.ny}, {"ixy", ixy}});
  }
  const std::span<const double> thetas = angular.thetas();
  const std::span<const double> weights = angular.theta_weights();
  const std::span<const double> phis = angular.phis();
  return {{"grids",
           {{"z_start", grids.z_scattering.start()},
            {"dz", grids.z_scattering.spacing()},
            {"nz_potential", grids.nz_potential},
            {"r_start", grids.r.start()},
            {"dr", grids.r.spacing()},
            {"nr", p.nr},
            {"nx", p.nx},
            {"ny", p.ny}}},
          {"theta_nodes_rad", std::vector<double>(thetas.begin(), thetas.end())},
          {"theta_weights", std::vector<double>(weights.begin(), weights.end())},
          {"phi_nodes_rad", std::vector<double>(phis.begin(), phis.end())},
          {"homonuclear", angular.homonuclear()},
          {"sites", sites},
          {"limits",
           {{"use_switching", p.use_switching},
            {"z_min_switch", p.z_min_switch},
            {"z_max_switch", p.z_max_switch},
            {"impose_v_min", p.impose_v_min},
            {"v_min_hartree", config.derived.v_min},
            {"impose_v_max", p.impose_v_max},
            {"v_max_hartree", config.derived.v_max}}},
          {"hartree_to_ev", units::legacy.hartree_to_ev(1.0)}};
}

// v_min.bin and v_min_node.bin: one (Z, r) point at a time, every (x, y) point.
void write_minimum(ArraysWriter& files, const PotentialGrid& potential, int nzr, int n_xy, int n_phi) {
  std::vector<double> minimum(static_cast<std::size_t>(n_xy));
  std::vector<int> orientation(static_cast<std::size_t>(n_xy));
  std::vector<std::int32_t> nodes(2 * static_cast<std::size_t>(n_xy));
  for (int izr = 0; izr < nzr; ++izr) {
    minimum_over_orientations(potential, izr, minimum, orientation);
    for (int ixy = 0; ixy < n_xy; ++ixy) {
      nodes[2 * ixy] = orientation[ixy] / n_phi;      // itheta
      nodes[2 * ixy + 1] = orientation[ixy] % n_phi;  // iphi
    }
    files.file(0).write_array(std::span<const double>(minimum));
    files.file(1).write_array(std::span<const std::int32_t>(nodes));
  }
}

// v_sites.bin: per site, every (Z, r) point and orientation (itheta outer, iphi inner).
void write_sites(BinaryWriter& file, const PotentialGrid& potential, const InputParameters& p, int nzr,
                 int n_orientations) {
  std::vector<double> row(static_cast<std::size_t>(n_orientations));
  for (const PlotSite& site : p.plot_sites) {
    const int ixy = plot_site_ixy(site, p.nx, p.ny);
    for (int izr = 0; izr < nzr; ++izr) {
      for (int itp = 0; itp < n_orientations; ++itp) row[itp] = potential(itp, ixy, izr);
      file.write_array(std::span<const double>(row));
    }
  }
}

}  // namespace

void write_potential_files(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
                           const AngularGrid& angular, const PotentialGrid& potential) {
  const InputParameters& p = config.input;
  ArraysWriter files(output_dir / "potential", "potential", potential_datasets(config, grids, angular),
                     potential_manifest_fields(config, grids, angular));
  const int nzr = grids.nz_potential * p.nr;  // izr = ir + nr * iz, i.e. (Z, r) in C order
  write_minimum(files, potential, nzr, p.nx * p.ny, angular.n_phi());
  if (!p.plot_sites.empty()) write_sites(files.file(2), potential, p, nzr, angular.n_theta() * angular.n_phi());
  files.finish();
}

std::optional<PlotDataSize> potential_data_size(const Config& config, const SpatialGrids& grids,
                                                const AngularGrid& angular) {
  if (!config.input.write_potential) return std::nullopt;
  std::uintmax_t bytes = 0;
  for (const DatasetSpec& d : potential_datasets(config, grids, angular)) bytes += dataset_bytes(d);
  return PlotDataSize{"potential/", "writePotential", bytes, "written once at set-up (no time-loop work)"};
}

}  // namespace spodvr
