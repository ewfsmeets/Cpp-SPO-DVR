// Screen output and data files of the fine analysis. Fortran: the stdout of AsympInfo.f90
// (nread = 2), Analyzor.f90 probaf (fil10, unit 3) and StateToStateProbs (nread = 2)
#include "report/fine_analysis_report.hpp"

#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include "report/analysis_files.hpp"
#include "report/analysis_report.hpp"

namespace spodvr {

void print_fine_analysis_grid(std::ostream& out, const Style& style, const Config& config, const FineAnalysisGrid& g) {
  const DerivedParameters& d = config.derived;
  print_section(out, style, "FINE ANALYSIS");
  KeyValueList list;
  list.add("BK records in the file", std::to_string(g.n_records_in_file));
  list.add("Records up to tEnd (ndatol)", std::to_string(g.n_records_available));
  list.add("Records used in the transform (ndused)", std::to_string(g.n_records_used));
  list.add("Effective propagation time", fixed(g.effective_time, 4), "a.u.");
  list.add("Spectral resolution", fixed(g.spectral_resolution * 1000.0, 6), "meV");
  list.add("Meaningful energy range", fixed(g.meaningful_min, 6) + " .. " + fixed(g.meaningful_max, 6), "eV");
  list.add("95% confidence interval", fixed(d.confidence_e_min, 6) + " .. " + fixed(d.confidence_e_max, 6), "eV");
  list.add("Grid points in the interval (nen95)", std::to_string(g.size()));
  if (g.size() > 0)
    list.add("Grid energies (E_z)", fixed(g.z_energy(0), 6) + " .. " + fixed(g.z_energy(g.size() - 1), 6), "eV");
  list.print(out, style);
  if (g.n_records_used != g.n_records_available)
    out << "\n  The transform uses the first " << g.n_records_used << " records (the largest allowed FFT length).\n";
}

void print_fine_analysis_results(std::ostream& out, const Style& style, const Config& config,
                                 const AnalysisBasis& basis, const FineAnalysisResults& r) {
  const int n = r.grid.size();
  if (n == 0) {
    out << "\n  " << style.warning("Warning")
        << ": no grid point lies in the 95% confidence interval (nen95 = 0), so the fine analysis gives no "
           "results. The propagation is probably too short.\n";
    return;
  }

  // Spectrum, thinned to about 50 rows.
  out << '\n' << style.emphasis("  Spectrum on the fine-analysis grid") << "\n\n";
  const int initial = basis.initial_state() + basis.n_states() * (config.derived.v_in_index - 1);
  const int stride = n > 50 ? (n + 49) / 50 : 1;
  Table table({"k", "E_z (eV)", "P_scattered", "P_reaction", "phase initial (rad)"});
  for (int i = 0; i < n; i += stride) {
    const EnergyResult& e = r.spectrum.energies[i];
    table.add_row({std::to_string(r.grid.first_index + i), fixed(e.z_energy, 7), fixed(e.p_total, 7),
                   scientific(e.p_reaction, 6), fixed(e.phase[initial], 6)});
  }
  table.print(out, style);
  if (stride > 1) out << "  1 in " << stride << " grid points is shown; all are in fine_analysis/spectrum/.\n";
  if (r.large_phase_steps > 0)
    out << "\n  " << style.warning("Warning") << ": " << r.large_phase_steps
        << " populated channel(s) change their phase by more than pi/2 between grid points; the interpolation "
           "to the analysis energies is unreliable (propagate longer).\n";

  if (!r.skipped.empty()) {
    out << "\n  Analysis energies without a fine-analysis result:\n";
    for (const SkippedEnergy& s : r.skipped) out << "    E_z = " << fixed(s.z_energy, 6) << " eV: " << s.reason << '\n';
  }
  if (!r.interpolated.energies.empty()) {
    out << '\n' << style.emphasis("  Interpolated to the analysis energies") << '\n';
    print_scattering_results(out, style, config, basis, r.interpolated);
  }
}

std::vector<std::filesystem::path> write_fine_analysis_files(const std::filesystem::path& dir, const Config& config,
                                                             const AnalysisBasis& basis, const FineAnalysisResults& r) {
  std::filesystem::create_directories(dir / "spectrum");
  std::vector<std::filesystem::path> files =
      write_scattering_files(dir / "spectrum", config, basis, r.spectrum, config.input.write_s_matrix);
  if (!r.interpolated.energies.empty()) {
    const auto interpolated = write_scattering_files(dir, config, basis, r.interpolated);
    files.insert(files.end(), interpolated.begin(), interpolated.end());
  }

  const FineAnalysisGrid& g = r.grid;
  nlohmann::json j;
  j["records_in_file"] = g.n_records_in_file;
  j["records_available"] = g.n_records_available;
  j["records_used"] = g.n_records_used;
  j["transform_length"] = g.transform_length;
  j["time_step"] = g.time_step;
  j["effective_time"] = g.effective_time;
  j["spectral_resolution_eV"] = g.spectral_resolution;
  j["first_index"] = g.first_index;
  j["last_index"] = g.last_index;
  j["n_points"] = g.size();
  j["large_phase_steps"] = r.large_phase_steps;
  nlohmann::json skipped = nlohmann::json::array();
  for (const SkippedEnergy& s : r.skipped) skipped.push_back({{"E_z", s.z_energy}, {"reason", s.reason}});
  j["skipped_energies"] = skipped;
  files.push_back(dir / "grid.json");
  std::ofstream file(files.back());
  if (!file) throw std::runtime_error("cannot open output file '" + files.back().string() + "'");
  file << j.dump(2) << '\n';
  return files;
}

}  // namespace spodvr
