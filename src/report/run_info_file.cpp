// run_info.json: the metadata of a run (program, input, grids, basis, lattice, absorbers, wave packet,
// analysis), for plots and readers of the other output files. C++ only (no Fortran counterpart).
#include "report/run_info_file.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "analysis/s_matrix.hpp"
#include "core/index_range.hpp"
#include "core/unit_conversion.hpp"
#include "grids/rotational_basis.hpp"
#include "input/validate_outputs.hpp"
#include "report/pes_cut_file.hpp"

namespace spodvr {

namespace {

// Keeps the fields in the order they are written below (nlohmann::json would sort them by name).
using json = nlohmann::ordered_json;

// An optional value as JSON: the value, or null.
template <typename T>
json or_null(const std::optional<T>& value) {
  return value ? json(*value) : json(nullptr);
}

std::string run_mode_name(int run_mode) {
  switch (run_mode) {
    case 0: return "new";
    case 1: return "restart";
    case 2: return "analysis only";
    case 6: return "potential only";
    default: return "unknown";
  }
}

json times(const InputParameters& p, const DerivedParameters& d) {
  return {{"t_start", p.t_start},         {"t_end", p.t_end},
          {"dt_step", p.time_step},       {"n_steps", p.n_steps},
          {"dt_analysis", p.dt_analysis}, {"steps_per_analysis", d.steps_per_analysis},
          {"dt_averages", p.dt_averages}, {"analyses_per_average", d.analyses_per_average}};
}

json molecule(const InputParameters& p, const DerivedParameters& d) {
  return {{"mass1_amu", p.mass1},
          {"mass2_amu", p.mass2},
          {"homonuclear", p.homonuclear},
          {"total_mass", d.total_mass},
          {"reduced_mass", d.reduced_mass}};
}

// The initial rovibrational state; its energy E(vin, jin) needs the gas phase.
json initial_state(const InputParameters& p, std::optional<std::reference_wrapper<const GasPhase>> gas) {
  json state = {{"v", p.v_in}, {"j", p.j_in}, {"mJ", p.mj_in}, {"energy", nullptr}, {"energy_ev", nullptr}};
  if (gas) {
    state["energy"] = gas->get().initial.energy;
    state["energy_ev"] = units::legacy.hartree_to_ev(gas->get().initial.energy);
  }
  return state;
}

// The Gaussian packet (2α/π)^{1/4} exp(−α(Z−Z0)² − i kz Z): sigma_z is the width of |ψ|²; kx, ky are
// components along the reciprocal lattice axes (process.cpp), k_parallel the Cartesian magnitude.
json wave_packet(const InputParameters& p, const DerivedParameters& d) {
  return {{"e_min_ev", p.e_min_wave_packet},
          {"e_max_ev", p.e_max_wave_packet},
          {"translational_energy_ev", d.translational_energy},
          {"z0", p.z0},
          {"alpha", d.wave_packet_alpha},
          {"sigma_z", 1.0 / (2.0 * std::sqrt(d.wave_packet_alpha))},
          {"k_total", d.k_total},
          {"kz", d.kz},
          {"kx", d.kx},
          {"ky", d.ky},
          {"k_parallel", d.k_total * std::sin(d.theta_incidence)},
          {"theta_incidence_rad", d.theta_incidence},
          {"phi_incidence_rad", d.phi_incidence},
          {"e_parallel_ev", p.e_parallel},
          {"phi_parallel_deg", p.phi_parallel}};
}

// The energies of the on-the-fly S matrix (they need E(vin, jin) from the gas phase).
json analysis_energies(const Config& config, std::optional<std::reference_wrapper<const GasPhase>> gas) {
  if (!gas) return {{"collision_ev", nullptr}, {"total_ev", nullptr}};
  const AnalysisEnergies energies = make_analysis_energies(config, gas->get());
  return {{"collision_ev", energies.collision}, {"total_ev", energies.total}};
}

json grid_sizes(const InputParameters& p, const DerivedParameters& d, const SpatialGrids& grids,
                const AngularGrid& angular) {
  return {{"z_start", grids.z_scattering.start()},
          {"dz", grids.z_scattering.spacing()},
          {"nz", grids.z_scattering.size()},
          {"nz_specular", grids.z_specular.size()},
          {"r_start", grids.r.start()},
          {"dr", grids.r.spacing()},
          {"nr", grids.r.size()},
          {"nx", p.nx},
          {"ny", p.ny},
          {"dx", d.dx},
          {"dy", d.dy},
          {"j_max", p.j_max},
          {"mj_max", d.mj_max_basis},
          {"n_theta", angular.n_theta()},
          {"n_phi", angular.n_phi()},
          {"nz_potential", grids.nz_potential},
          {"xy_volume_element", d.xy_volume_element},
          {"unit_cell_area", d.unit_cell_area},
          {"unit_cell_norm", d.unit_cell_norm}};
}

// [J, mJ] of every rotational state irot, in the order of the wave function's first index.
json rotational_basis(const InputParameters& p, const AngularGrid& angular) {
  const RotationalBasis basis(angular, p.j_in, p.mj_in);
  json states = json::array();
  for (int irot = 0; irot < basis.size(); ++irot) states.push_back({basis.j(irot), basis.mj(irot)});
  return {{"states", states}, {"initial_index", basis.initial_index()}};
}

json lattice(const DerivedParameters& d) {
  return {{"a", d.system.lattice_constant_x},
          {"b", d.system.lattice_constant_y},
          {"skew_angle_rad", d.system.skewing_angle},
          {"point_group", d.system.point_group},
          {"delta_kx", d.delta_kx},
          {"delta_ky", d.delta_ky}};
}

json analysis_surface(const InputParameters& p, const DerivedParameters& d, const SpatialGrids& grids) {
  return {{"z_ana", grids.z_analysis},
          {"iz_ana", grids.iz_analysis},
          {"z_ana_specular", or_null(p.z_analysis_specular)},
          {"dynamic_specular_analysis", p.dynamical_zsp_analysis},
          {"fixed_z_analysis_specular", d.fixed_z_analysis_specular}};
}

// One absorber: its range, effective strength (hartree, after the energy table) and grid points.
json absorber(const std::string& name, const std::string& grid, const std::string& coordinate, double start, double end,
              double strength, const std::optional<IndexRange>& points) {
  json range = points ? json{{"first", points->first}, {"last", points->last}} : json(nullptr);
  return {{"name", name}, {"grid", grid},         {"coordinate", coordinate}, {"start", start},
          {"end", end},   {"strength", strength}, {"points", range}};
}

json absorbers(const InputParameters& p, const DerivedParameters& d, const SpatialGrids& grids) {
  json list = json::array();
  list.push_back(absorber("ZOpt", "scattering", "Z", p.z_opt_start, p.z_opt_end, d.z_opt_factor, grids.z_optical));
  list.push_back(absorber("ROpt", "both", "r", p.r_opt_start, p.r_opt_end, d.r_opt_factor, grids.r_optical));
  list.push_back(
      absorber("ZSpOpt", "specular", "Z", p.zsp_opt_start, p.zsp_opt_end, d.zsp_opt_factor, grids.z_specular_optical));
  if (p.z_opt_low.enabled())
    list.push_back(absorber("ZOptLow", "scattering", "Z", *p.z_opt_low.start, *p.z_opt_low.end, d.z_opt_low_factor,
                            grids.z_optical_low));
  if (p.zsp_opt_low.enabled())
    list.push_back(absorber("ZSpOptLow", "specular", "Z", *p.zsp_opt_low.start, *p.zsp_opt_low.end,
                            d.zsp_opt_low_factor, grids.z_specular_optical_low));
  return {{"power", d.optical_potential_power}, {"list", list}};
}

// The plot-data keys of the run and its sites (with their grid points); written only when a plot key is set, so
// runs without them keep their run_info.json.
json plot_data(const InputParameters& p) {
  json sites = json::array();
  for (const PlotSite& site : p.plot_sites) {
    const int ixy = plot_site_ixy(site, p.nx, p.ny);
    sites.push_back({{"name", site.name}, {"x", site.x}, {"y", site.y}, {"ix", ixy / p.ny}, {"iy", ixy % p.ny}});
  }
  return {{"write_densities", p.write_densities}, {"write_potential", p.write_potential}, {"sites", sites}};
}

// The cut of a write2DPES run (written only by such a run, which stops after it).
json pes_cut(const InputParameters& p) {
  return {{"file", std::string(pes_cut_file_name)},
          {"x", p.x_2d_pes},
          {"y", p.y_2d_pes},
          {"theta_deg", p.theta_2d_pes},
          {"phi_deg", p.phi_2d_pes}};
}

json potential_settings(const InputParameters& p, const DerivedParameters& d) {
  return {{"use_switching", p.use_switching},
          {"z_min_switch", p.z_min_switch},
          {"z_max_switch", p.z_max_switch},
          {"impose_v_min", p.impose_v_min},
          {"v_min", d.v_min},
          {"v_min_ev", p.v_min},
          {"impose_v_max", p.impose_v_max},
          {"v_max", d.v_max},
          {"v_max_ev", p.v_max},
          {"use_symmetry", p.use_symmetry},
          {"symmetry_tolerance_ev", p.symmetry_tolerance}};
}

// The factors the program uses for its conversions (the Fortran's legacy set, units::legacy).
json unit_factors() {
  return {{"hartree_to_ev", units::legacy.hartree_to_ev(1.0)},
          {"amu_to_me", units::legacy.amu_to_me(1.0)},
          {"bohr_to_angstrom", units::legacy.bohr_to_angstrom(1.0)},
          {"hartree_to_wavenumber", units::legacy.hartree_to_wavenumber(1.0)}};
}

}  // namespace

void write_run_info(const std::filesystem::path& path, const Config& config, const SpatialGrids& grids,
                    const AngularGrid& angular, const ProgramInfo& program, const std::filesystem::path& input_file,
                    std::optional<std::reference_wrapper<const GasPhase>> gas) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  json j;
  j["format_version"] = 1;
  j["program"] = {{"name", "spodvr"}, {"version", program.version}, {"git", program.git}};
  j["units"] = {{"length", "bohr"},
                {"time", "atomic units"},
                {"energy", "hartree, unless the field name ends in _ev"},
                {"mass", "electron masses, unless the field name ends in _amu"},
                {"angle", "rad, unless the field name ends in _deg"}};
  j["title"] = p.title;
  j["system_id"] = p.system_id;
  j["run_mode"] = {{"code", d.run_mode}, {"name", run_mode_name(d.run_mode)}};
  j["input_file"] = input_file.string();
  j["new_calculation"] = p.new_calculation;
  j["restart_dir"] = p.restart_dir.empty() ? json(nullptr) : json(p.restart_dir);
  j["times"] = times(p, d);
  j["molecule"] = molecule(p, d);
  j["initial_state"] = initial_state(p, gas);
  j["wave_packet"] = wave_packet(p, d);
  j["analysis_energies"] = analysis_energies(config, gas);
  j["grids"] = grid_sizes(p, d, grids, angular);
  j["rotational_basis"] = rotational_basis(p, angular);
  j["lattice"] = lattice(d);
  j["z_asymptotic"] = d.z_asymptotic;
  j["z_asymptotic_unswitched"] = d.z_asymptotic_unswitched;
  j["analysis"] = analysis_surface(p, d, grids);
  j["absorbers"] = absorbers(p, d, grids);
  j["potential"] = potential_settings(p, d);
  j["unit_factors"] = unit_factors();
  if (p.write_densities || p.write_potential || !p.plot_sites.empty()) j["plot_data"] = plot_data(p);
  if (p.write_2d_pes) j["pes_cut"] = pes_cut(p);
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << j.dump(2) << '\n';
}

}  // namespace spodvr
