// energy_resolved/ (EnergyResolvedSites, EnergyResolvedEnergies): the energy-resolved sums of the run as an array
// directory, their reader for a restart, and their projected size and cost. C++ only.
#include "report/energy_resolved_files.hpp"

#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "input/validate_outputs.hpp"
#include "report/format.hpp"

namespace spodvr {

namespace {

constexpr std::string_view directory_name = "energy_resolved";

// The cost per complex value of one site's sums (nZ nR nrot) per (site, energy) pair and analysis step, the column
// copy and the accumulation together: measured on the thesis DZ grid with one site and two energies
// (docs/performance.md "Plot data on the thesis DZ grid"; the CPU's serial strided gather dominates there).
constexpr double pair_seconds_per_value_cpu = 8.3e-9;
constexpr double pair_seconds_per_value_gpu = 2.7e-9;

// The datasets, in the order the files are written: the final sums, the partial sums and the scalars.
std::vector<DatasetSpec> energy_resolved_datasets(const EnergyResolvedSums& s) {
  const int n_sites = static_cast<int>(s.site_names.size());
  const int n_e = static_cast<int>(s.collision_ev.size());
  const std::vector<int> site_shape = {n_sites, n_e, s.n_z, s.n_r, s.n_rot};
  const std::vector<std::string> site_axes = {"site", "energy", "Z", "r", "rot"};
  const std::vector<int> specular_shape = {n_e, s.n_z_specular, s.n_r};
  const std::vector<std::string> specular_axes = {"energy", "Z_specular", "r"};
  const std::string sum_units =
      "bohr^-2 a.u. (psi in the program's normalisation, sum |psi|^2 dZ dr dA = 1, times time)";
  return {
      {"sums.bin", "<c16", site_shape, site_axes, sum_units,
       "sum over the analysis steps of psi_scat at each site times exp(i E_tot t) dt_ana (the whole run)"},
      {"specular_sums.bin", "<c16", specular_shape, specular_axes, sum_units,
       "the same sum of the specular wave chi (the initial channel), per energy"},
      {"sums_partial.bin", "<c16", site_shape, site_axes, sum_units,
       "sums.bin up to the first analysis step at or after 0.8 tEnd of this segment (partial_time)"},
      {"specular_sums_partial.bin", "<c16", specular_shape, specular_axes, sum_units,
       "specular_sums.bin up to partial_time"},
      {"partial_time.bin", "<f8", {}, {}, "a.u.", "the time of the last analysis step in the partial sums"},
      {"last_time.bin", "<f8", {}, {}, "a.u.", "the time of the last analysis step in the final sums (tEnd)"},
      {"steps_summed.bin", "<i4", {}, {}, "", "the analysis steps in the final sums (all segments of the run)"},
  };
}

nlohmann::ordered_json sites_json(const EnergyResolvedSums& s, std::span<const std::complex<double>> s2s, int ny) {
  nlohmann::ordered_json sites = nlohmann::ordered_json::array();
  for (std::size_t k = 0; k < s.site_names.size(); ++k) {
    const int ixy = s.site_ixy[k];
    const std::complex<double> factor = s2s[std::size_t(ixy)];
    sites.push_back({{"name", s.site_names[k]},
                     {"ix", ixy / ny},
                     {"iy", ixy % ny},
                     {"ixy", ixy},
                     {"specular_to_scattering", {factor.real(), factor.imag()}}});
  }
  return sites;
}

nlohmann::ordered_json manifest_fields(const Config& config, const SpatialGrids& grids, const EnergyResolvedSums& s,
                                       std::span<const std::complex<double>> s2s) {
  const DerivedParameters& d = config.derived;
  return {{"sites", sites_json(s, s2s, config.input.ny)},
          {"collision_ev", s.collision_ev},
          {"total_hartree", s.total_hartree},
          {"analysis_time_step", s.analysis_time_step},
          {"segment_t_start", config.input.t_start},
          {"partial_fraction", partial_sums_fraction},
          {"grids",
           {{"z_start", grids.z_scattering.start()},
            {"dz", grids.z_scattering.spacing()},
            {"nz", s.n_z},
            {"nz_specular", s.n_z_specular},
            {"r_start", grids.r.start()},
            {"dr", grids.r.spacing()},
            {"nr", s.n_r},
            {"n_rot", s.n_rot}}},
          {"wave_packet", {{"alpha", d.wave_packet_alpha}, {"kz", d.kz}, {"z0", config.input.z0}}}};
}

}  // namespace

EnergyResolvedFiles::EnergyResolvedFiles(const std::filesystem::path& output_dir, const Config& config,
                                         const SpatialGrids& grids, const EnergyResolvedSums& layout,
                                         std::span<const std::complex<double>> s2s)
    : arrays_(output_dir / directory_name, "energy_resolved", energy_resolved_datasets(layout),
              manifest_fields(config, grids, layout, s2s)) {}

void EnergyResolvedFiles::write_partial(const EnergyResolvedSums& sums) {
  arrays_.file(2).write_array(std::span<const std::complex<double>>(sums.sums));
  arrays_.file(3).write_array(std::span<const std::complex<double>>(sums.specular_sums));
  arrays_.file(4).write(sums.last_time);
  partial_written_ = true;
}

void EnergyResolvedFiles::write_final(const EnergyResolvedSums& sums) {
  if (!partial_written_) write_partial(sums);  // a segment without analysis steps (nStep = 0): its sums are both
  arrays_.file(0).write_array(std::span<const std::complex<double>>(sums.sums));
  arrays_.file(1).write_array(std::span<const std::complex<double>>(sums.specular_sums));
  arrays_.file(5).write(sums.last_time);
  arrays_.file(6).write(static_cast<std::int32_t>(sums.steps_summed));
  arrays_.finish();
}

namespace {

nlohmann::json restart_manifest(const std::filesystem::path& restart_dir) {
  const std::filesystem::path path = restart_dir / directory_name / "manifest.json";
  std::ifstream file(path);
  if (!file) throw CheckpointError("cannot open '" + path.string() + "'");
  const nlohmann::json manifest = nlohmann::json::parse(file, nullptr, false);
  if (manifest.is_discarded()) throw CheckpointError("'" + path.string() + "' is not JSON");
  return manifest;
}

// The sites (name and grid point) and the energies of RestartDir's sums must be this run's: InputError otherwise.
void check_same_keys(const nlohmann::json& manifest, const std::filesystem::path& dir, const InputParameters& p) {
  std::vector<std::pair<std::string, int>> sites;
  for (const auto& site : manifest.value("sites", nlohmann::json::array()))
    sites.emplace_back(site.value("name", ""), site.value("ixy", -1));
  std::vector<std::pair<std::string, int>> expected;
  for (const std::string& name : p.energy_resolved_sites)
    for (const PlotSite& site : p.plot_sites)
      if (site.name == name) expected.emplace_back(name, plot_site_ixy(site, p.nx, p.ny));
  if (sites != expected || manifest.value("collision_ev", std::vector<double>{}) != p.energy_resolved_energies)
    throw InputError("EnergyResolvedSites", "the sums in '" + dir.string() +
                                                "' have other sites (names or grid points) or energies than this "
                                                "restart; EnergyResolvedSites, their PlotSites and "
                                                "EnergyResolvedEnergies must stay the same");
}

}  // namespace

void check_energy_resolved_restart(const std::filesystem::path& restart_dir, const InputParameters& p) {
  check_same_keys(restart_manifest(restart_dir), restart_dir / directory_name, p);
}

EnergyResolvedSums read_energy_resolved_sums(const std::filesystem::path& restart_dir, const Config& config,
                                             const EnergyResolvedSums& expected, int earlier_steps) {
  const std::filesystem::path dir = restart_dir / directory_name;
  const nlohmann::json manifest = restart_manifest(restart_dir);
  check_same_keys(manifest, dir, config.input);
  const auto& g = manifest.value("grids", nlohmann::json::object());
  if (g.value("nz", -1) != expected.n_z || g.value("nr", -1) != expected.n_r ||
      g.value("n_rot", -1) != expected.n_rot || g.value("nz_specular", -1) != expected.n_z_specular)
    throw CheckpointError("the energy-resolved sums in '" + dir.string() + "' are on other grids than this restart");
  EnergyResolvedSums s = expected;
  BinaryReader(dir / "sums.bin").read_array(std::span<std::complex<double>>(s.sums));
  BinaryReader(dir / "specular_sums.bin").read_array(std::span<std::complex<double>>(s.specular_sums));
  s.last_time = BinaryReader(dir / "last_time.bin").read<double>();
  s.steps_summed = BinaryReader(dir / "steps_summed.bin").read<std::int32_t>();
  const double t_start = config.input.t_start;  // as the checkpoint's check (run.cpp)
  if (fortran_compare::rne(s.last_time, t_start))
    throw CheckpointError("tStart = " + full_precision(t_start) + " but the energy-resolved sums in '" + dir.string() +
                          "' end at t = " + full_precision(s.last_time));
  // Sums left by an earlier run in the same directory: they do not belong to the BK records read with them.
  if (s.steps_summed != earlier_steps)
    throw CheckpointError("the energy-resolved sums in '" + dir.string() + "' hold " + std::to_string(s.steps_summed) +
                          " analysis steps, the BK file " + std::to_string(earlier_steps) +
                          ": they belong to another run");
  return s;
}

void print_energy_resolved_cost(std::ostream& out, std::size_t pairs, double seconds, double propagation_seconds) {
  out << "  Energy-resolved sums (EnergyResolvedSites): " << pairs << " (site, energy) pairs, " << fixed(seconds, 2)
      << " s, " << percent(seconds / propagation_seconds) << " of the propagation\n";
}

std::optional<PlotDataSize> energy_resolved_data_size(const Config& config, const SpatialGrids& grids, int n_rot) {
  const InputParameters& p = config.input;
  if (p.energy_resolved_sites.empty()) return std::nullopt;
  const std::size_t pairs = p.energy_resolved_sites.size() * p.energy_resolved_energies.size();
  const double column = double(p.nz) * p.nr * n_rot;
  const double step_values = column * p.nx * p.ny * config.derived.steps_per_analysis;
  const double cpu_share =
      double(pairs) * pair_seconds_per_value_cpu * column / (step_seconds_per_value_cpu * step_values);
  const double gpu_share =
      double(pairs) * pair_seconds_per_value_gpu * column / (step_seconds_per_value_gpu * step_values);
  std::string cost = std::to_string(pairs) + " (site, energy) pairs, about " + percent(cpu_share) + " (CPU) or " +
                     percent(gpu_share) + " (GPU) of each analysis step";
  if (std::max(cpu_share, gpu_share) > energy_resolved_warning_fraction)
    cost += "; above 10 %: fewer sites or energies are cheaper";
  // Two copies of the sums (final and partial; the in-memory column buffer left out) and the scalars.
  const std::size_t sums_values = energy_resolved_memory_values(p, n_rot, static_cast<int>(grids.z_specular.size())) -
                                  std::size_t(p.nz) * p.nr * n_rot;
  const std::uintmax_t bytes =
      2 * sums_values * sizeof(std::complex<double>) + sizeof(double) + sizeof(double) + sizeof(std::int32_t);
  return PlotDataSize{"energy_resolved/", "EnergyResolvedSites", bytes, cost};
}

}  // namespace spodvr
