// Data files of the asymptotic analysis: tidy tab-separated tables and JSON metadata.
// Replaces the Fortran fort.25 / fil10 / unit 3 output (Analyzor.f90 smatrix_combine, probaf).
#include "report/analysis_files.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "report/format.hpp"

namespace spodvr {

namespace {

// A tab-separated table: '#' comment lines, one header line, then rows.
class TsvFile {
 public:
  TsvFile(const std::filesystem::path& path, const std::vector<std::string>& comments,
          const std::vector<std::string>& columns)
      : file_(path) {
    if (!file_) throw std::runtime_error("cannot open output file '" + path.string() + "'");
    for (const std::string& line : comments) file_ << "# " << line << '\n';
    write_row(columns);
  }

  void write_row(const std::vector<std::string>& cells) {
    for (std::size_t i = 0; i < cells.size(); ++i) file_ << (i == 0 ? "" : "\t") << cells[i];
    file_ << '\n';
  }

 private:
  std::ofstream file_;
};

const std::string energy_comment =
    "Energies in eV. E_coll: collision (translational) energy; E_total: E_coll + E(vin, jin).";

// Labels v, j, mJ of a state, as table cells.
std::vector<std::string> rovib_cells(const AnalysisBasis& basis, int is, int v) {
  return {std::to_string(v), std::to_string(basis.j(is)), std::to_string(basis.mj(is))};
}

void write_summary(const std::filesystem::path& path, const AnalysisBasis& basis, const ScatteringResults& results) {
  std::vector<std::string> columns = {"E_coll", "E_total", "E_z", "P_scattered", "P_reaction"};
  for (int v = 0; v < basis.n_v(); ++v) columns.push_back("P_v" + std::to_string(v));
  for (const char* name : {"E_vib_mean", "E_rot_mean", "E_dif_mean"}) columns.push_back(name);
  TsvFile file(path,
               {"Scattering summary, one row per analysis energy.", energy_comment,
                "E_z: E_coll minus the initial parallel kinetic energy.",
                "P_scattered: sum of all state-to-state probabilities; P_reaction = 1 - P_scattered;",
                "P_v<v>: summed over rotation and diffraction; E_*_mean: average final vibrational, rotational",
                "and diffraction (parallel kinetic) energy of the scattered molecules (eV)."},
               columns);
  for (const EnergyResult& e : results.energies) {
    std::vector<std::string> row = {full_precision(e.collision_energy), full_precision(e.total_energy),
                                    full_precision(e.z_energy), full_precision(e.p_total),
                                    full_precision(e.p_reaction)};
    for (const double p : e.p_vibration) row.push_back(full_precision(p));
    for (const double x : {e.e_vibration, e.e_rotation, e.e_diffraction}) row.push_back(full_precision(x));
    file.write_row(row);
  }
}

void write_s_matrix(const std::filesystem::path& path, const AnalysisBasis& basis, const ScatteringResults& results) {
  TsvFile file(path,
               {"S-matrix elements of the open channels (closed channels are omitted).", energy_comment,
                "Final state: vibration v, rotation j, mJ, diffraction channel (n, m) of the given order.",
                "P = |S|^2; phase = arg S in rad (0 if Re S = 0, as in the Fortran code)."},
               {"E_coll", "E_total", "v", "j", "mJ", "n", "m", "order", "Re_S", "Im_S", "P", "phase"});
  const std::vector<int> order = basis.regular_state_order();
  for (const EnergyResult& e : results.energies)
    for (int v = 0; v < basis.n_v(); ++v)
      for (const int is : order) {
        const int k = is + basis.n_states() * v;
        if (!e.open[k]) continue;
        const DiffractionChannel& c = basis.channel(is);
        std::vector<std::string> row = {full_precision(e.collision_energy), full_precision(e.total_energy)};
        for (std::string& cell : rovib_cells(basis, is, v)) row.push_back(std::move(cell));
        for (const int q : {c.n, c.m, c.order}) row.push_back(std::to_string(q));
        for (const double x : {e.s_matrix[k].real(), e.s_matrix[k].imag(), e.probability[k], e.phase[k]})
          row.push_back(full_precision(x));
        file.write_row(row);
      }
}

void write_rovib(const std::filesystem::path& path, const AnalysisBasis& basis, const ScatteringResults& results) {
  TsvFile file(path, {"Rovibrational state-resolved probabilities, summed over diffraction.", energy_comment},
               {"E_coll", "v", "j", "mJ", "P"});
  for (const EnergyResult& e : results.energies)
    for (int v = 0; v < basis.n_v(); ++v)
      for (const int irotana : basis.regular_rotation_order()) {
        std::vector<std::string> row = {full_precision(e.collision_energy)};
        for (std::string& cell : rovib_cells(basis, irotana * basis.n_diffraction(), v)) row.push_back(std::move(cell));
        row.push_back(full_precision(e.p_rotation[irotana + basis.n_rotation() * v]));
        file.write_row(row);
      }
}

void write_diffraction(const std::filesystem::path& path, const AnalysisBasis& basis,
                       const ScatteringResults& results) {
  TsvFile file(path,
               {"Probabilities per diffraction order, for every rovibrational state.", energy_comment,
                "P is summed over all diffraction channels (n, m) of the given order."},
               {"E_coll", "v", "j", "mJ", "order", "P"});
  const int n_orders = basis.max_order() + 1;
  for (const EnergyResult& e : results.energies)
    for (int v = 0; v < basis.n_v(); ++v)
      for (const int irotana : basis.regular_rotation_order())
        for (int order = 0; order < n_orders; ++order) {
          std::vector<std::string> row = {full_precision(e.collision_energy)};
          for (std::string& cell : rovib_cells(basis, irotana * basis.n_diffraction(), v))
            row.push_back(std::move(cell));
          row.push_back(std::to_string(order));
          row.push_back(full_precision(e.p_order[order + n_orders * (irotana + basis.n_rotation() * v)]));
          file.write_row(row);
        }
}

void write_metadata(const std::filesystem::path& path, const Config& config, const AnalysisBasis& basis,
                    const ScatteringResults& results, const std::vector<std::filesystem::path>& files) {
  const InputParameters& p = config.input;
  const int is0 = basis.initial_state();
  nlohmann::json j;
  j["program"] = {{"name", "spodvr"}, {"version", SPODVR_VERSION}};
  j["title"] = p.title;
  j["system_id"] = p.system_id;
  j["units"] = {{"energy", "eV"}, {"length", "bohr"}, {"time", "atomic units"}, {"phase", "rad"}};
  j["initial_state"] = {{"v", p.v_in}, {"j", p.j_in}, {"mJ", p.mj_in}, {"n", 0}, {"m", 0}};
  j["analysis"] = {{"z_analysis", p.z_analysis},
                   {"z0", p.z0},
                   {"n_vibrational_levels", basis.n_v()},
                   {"j_max", basis.j_max()},
                   {"max_diffraction_order", basis.max_order()},
                   {"n_rotational_states", basis.n_rotation()},
                   {"n_diffraction_channels", basis.n_diffraction()},
                   {"n_states", basis.n_states()},
                   {"initial_state_index", is0}};
  j["energies"] = {{"count", results.energies.size()},
                   {"collision_min", results.energies.empty() ? 0.0 : results.energies.front().collision_energy},
                   {"collision_max", results.energies.empty() ? 0.0 : results.energies.back().collision_energy}};
  nlohmann::json list = nlohmann::json::array();
  for (const auto& f : files) list.push_back(f.filename().string());
  j["files"] = list;
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << j.dump(2) << '\n';
}

}  // namespace

std::vector<std::filesystem::path> write_scattering_files(const std::filesystem::path& dir, const Config& config,
                                                          const AnalysisBasis& basis,
                                                          const ScatteringResults& results) {
  std::vector<std::filesystem::path> files = {dir / "scattering_summary.tsv", dir / "s_matrix.tsv"};
  write_summary(files[0], basis, results);
  write_s_matrix(files[1], basis, results);
  if (config.input.write_probabilities) {
    files.push_back(dir / "probabilities_rovib.tsv");
    write_rovib(files.back(), basis, results);
  }
  if (config.input.write_rd_probabilities) {
    files.push_back(dir / "probabilities_diffraction.tsv");
    write_diffraction(files.back(), basis, results);
  }
  files.push_back(dir / "analysis.json");
  write_metadata(files.back(), config, basis, results, files);
  return files;
}

}  // namespace spodvr
