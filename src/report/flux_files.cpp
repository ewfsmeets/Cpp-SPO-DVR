// Data files of the flux analysis: the reaction probability and its (x, y) density as tidy tables,
// and JSON metadata. Replaces the Fortran flux_reaction.dat and flux_density.dat (FluxAnalyzor.f90).
#include "report/flux_files.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "core/unit_conversion.hpp"
#include "report/format.hpp"
#include "report/tsv_file.hpp"
#include "symmetry/surface_lattice.hpp"

namespace spodvr {

namespace {

const std::string energy_comment =
    "E_z: energy of the motion normal to the surface (eV; Fortran ZScatEnergy), "
    "EMinFluxAna + n DeltaEFluxAna.";

// Fortran: computeFluxProbability, fileFluxReaction
void write_reaction(const std::filesystem::path& path, const FluxResults& results) {
  TsvFile file(path,
               {"Reaction probability from the flux through the surface r = r_c (rFluxAnalysis), one row per "
                "flux energy.",
                energy_comment, "P_flux: the reaction probability at E_z."},
               {"E_z", "P_flux"});
  for (std::size_t n = 0; n < results.energy.size(); ++n)
    file.write_row({full_precision(results.energy[n]), full_precision(results.probability[n])});
}

// Fortran: computeFluxProbability, fileFluxDensity (its header lines become the comments)
void write_density(const std::filesystem::path& path, const Config& config, const SpatialGrids& grids,
                   const FluxResults& results) {
  const SystemParameters& system = config.derived.system;
  const SurfaceLattice lattice(system);
  TsvFile file(
      path,
      {"Reaction probability density from the flux through r = r_c, one row per flux energy and (x, y) point.",
       energy_comment,
       "ix, iy: grid indices (0-based); x, y: skewed grid coordinates along the lattice vectors (bohr);",
       "X, Y: Cartesian coordinates (bohr), X = x + y cos(skewing angle), Y = y sin(skewing angle).",
       "density: P(x, y; E_z) (1/bohr^2); P_flux(E_z) in flux_reaction.tsv is its sum times dx dy sin(skewing angle).",
       "Surface grid: nx = " + std::to_string(grids.x.size()) + ", dx = " + full_precision(grids.x.spacing()) +
           " bohr, ny = " + std::to_string(grids.y.size()) + ", dy = " + full_precision(grids.y.spacing()) +
           " bohr, skewing angle " + full_precision(rad_to_deg(system.skewing_angle)) + " deg, point group " +
           system.point_group + "."},
      {"E_z", "ix", "iy", "x", "y", "X", "Y", "density"});
  for (int n = 0; n < static_cast<int>(results.energy.size()); ++n)
    for (int ix = 0; ix < grids.x.size(); ++ix)
      for (int iy = 0; iy < grids.y.size(); ++iy) {
        const CartesianVector cartesian = lattice.to_cartesian(grids.x[ix], grids.y[iy]);
        file.write_row({full_precision(results.energy[n]), std::to_string(ix), std::to_string(iy),
                        full_precision(grids.x[ix]), full_precision(grids.y[iy]), full_precision(cartesian.x),
                        full_precision(cartesian.y),
                        full_precision(results.density[results.density_index(ix, iy, n)])});
      }
}

void write_metadata(const std::filesystem::path& path, const Config& config, const SpatialGrids& grids,
                    const FluxResults& results, const std::vector<std::filesystem::path>& files) {
  const InputParameters& p = config.input;
  const SystemParameters& system = config.derived.system;
  const int ir_flux = *grids.ir_flux_analysis;
  nlohmann::json j;
  j["program"] = {{"name", "spodvr"}, {"version", SPODVR_VERSION}};
  j["title"] = p.title;
  j["system_id"] = p.system_id;
  j["units"] = {{"energy", "eV"}, {"length", "bohr"}, {"density", "1/bohr^2"}};
  j["initial_state"] = {{"v", p.v_in}, {"j", p.j_in}, {"mJ", p.mj_in}};
  j["flux_surface"] = {{"r", grids.r[ir_flux]}, {"ir", ir_flux}, {"z_points", grids.nz_potential}};
  j["energies"] = {{"count", results.energy.size()},
                   {"E_z_min", results.energy.empty() ? 0.0 : results.energy.front()},
                   {"E_z_max", results.energy.empty() ? 0.0 : results.energy.back()},
                   {"EMinFluxAna", p.e_min_flux_analysis},
                   {"EMaxFluxAna", p.e_max_flux_analysis},
                   {"DeltaEFluxAna", p.delta_e_flux_analysis}};
  j["surface_grid"] = {{"nx", grids.x.size()},
                       {"ny", grids.y.size()},
                       {"dx", grids.x.spacing()},
                       {"dy", grids.y.spacing()},
                       {"skewing_angle_deg", rad_to_deg(system.skewing_angle)},
                       {"point_group", system.point_group}};
  nlohmann::json list = nlohmann::json::array();
  for (const auto& f : files) list.push_back(f.filename().string());
  j["files"] = list;
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << j.dump(2) << '\n';
}

}  // namespace

std::vector<std::filesystem::path> write_flux_files(const std::filesystem::path& dir, const Config& config,
                                                    const SpatialGrids& grids, const FluxResults& results) {
  if (!grids.ir_flux_analysis) throw std::invalid_argument("write_flux_files: no flux surface (doFluxAnalysis)");
  const std::size_t n_energies = results.energy.size();
  if (results.probability.size() != n_energies || results.density_index.n0 != std::size_t(grids.x.size()) ||
      results.density_index.n1 != std::size_t(grids.y.size()) || results.density_index.n2 != n_energies ||
      results.density.size() != results.density_index.size())
    throw std::invalid_argument("write_flux_files: the results do not match the surface grid or their energies");
  std::vector<std::filesystem::path> files = {dir / "flux_reaction.tsv", dir / "flux_density.tsv"};
  write_reaction(files[0], results);
  write_density(files[1], config, grids, results);
  files.push_back(dir / "flux_analysis.json");
  write_metadata(files.back(), config, grids, results, files);
  return files;
}

}  // namespace spodvr
