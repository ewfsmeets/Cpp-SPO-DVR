// Program driver: input -> grids -> potential -> gas-phase states -> propagation -> output.
// Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

#include "app/banner.hpp"
#include "core/error.hpp"
#include "core/numerical_libraries.hpp"
#include "core/real_compare.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/checkpoint.hpp"
#include "io/fortran_files.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/propagation.hpp"
#include "report/analysis_files.hpp"
#include "report/analysis_report.hpp"
#include "report/propagation_report.hpp"
#include "report/sections.hpp"

namespace spodvr {

namespace {

// Checkpoint file names in an output / restart directory: native, and as written by the
// Fortran code (fort.9, fort.21; a Fortran restart reads them as fort.8, fort.22).
constexpr std::string_view checkpoint_file = "checkpoint.bin";
constexpr std::string_view bk_file = "bk_coefficients.bin";
constexpr std::string_view fortran_wave_function_file = "fort.9";
constexpr std::string_view fortran_bk_file = "fort.21";

// BK coefficients of the earlier runs, in time order.
using BkRecords = std::vector<std::pair<double, BkCoefficients>>;

// The first of the files that exists in dir (the Fortran input names fort.8 / fort.22 or the
// output names fort.9 / fort.21).
std::filesystem::path first_existing(const std::filesystem::path& dir, std::string_view a, std::string_view b) {
  if (std::filesystem::exists(dir / a)) return dir / a;
  if (std::filesystem::exists(dir / b)) return dir / b;
  throw CheckpointError("neither " + std::string(a) + " nor " + std::string(b) + " found in '" + dir.string() + "'");
}

// Restart (nread = 1): reads the wave packet and analysis surface of the earlier run and its BK
// coefficients, which are added to the S-matrix and returned (to be copied into the new BK files).
// Fortran: WaveFunction.f90 ReadWaveFunc, AsympInfo.f90 ReadOrWriteAsympInfo (nread = 1),
// Analyzor.f90 GetSMatrixElemsOfOldRun
BkRecords restore_from_checkpoint(std::ostream& out, const Config& config, WavePacketPropagation& propagation) {
  const InputParameters& p = config.input;
  const std::filesystem::path dir = p.restart_dir;
  BkRecords records;
  auto add_record = [&](double time, const BkCoefficients& bk) {
    propagation.add_earlier_analysis_step(time, bk);
    records.emplace_back(time, bk);
  };

  double checkpoint_time = 0.0;
  if (p.restart_format == "native") {
    const CheckpointInfo info = read_checkpoint(dir / checkpoint_file, config, propagation.wave_packet_for_restart());
    propagation.restore_analysis_surface(info.surface);
    checkpoint_time = info.time;
    const auto last = read_bk_file(dir / bk_file, config, propagation.analysis_basis(), add_record);
    if (last && fortran_compare::rne(*last, checkpoint_time))
      throw CheckpointError("the BK coefficient file in '" + dir.string() + "' ends at t = " + full_precision(*last) +
                            ", the checkpoint at t = " + full_precision(checkpoint_time));
  } else {
    // The Fortran wave function file has no time: the BK records give it (n deltat).
    const AnalysisSurfaceState fresh = propagation.analysis_surface().state();
    const AnalysisSurfaceState state = read_fortran_wave_function(first_existing(dir, "fort.8", "fort.9"),
                                                                  propagation.wave_packet_for_restart(), fresh);
    propagation.restore_analysis_surface(state);
    read_fortran_bk_file(first_existing(dir, "fort.22", "fort.21"), config, propagation.rotational_basis(),
                         propagation.analysis_basis(), add_record);
    checkpoint_time = records.empty() ? 0.0 : records.back().first;
  }
  if (fortran_compare::rne(checkpoint_time, p.t_start))
    throw CheckpointError("tStart = " + full_precision(p.t_start) + " but the checkpoint in '" + dir.string() +
                          "' is at t = " + full_precision(checkpoint_time));
  out << "\n  Restarted from " << dir.string() << " (" << p.restart_format
      << " files) at t = " << fixed(checkpoint_time, 2) << " a.u.\n";
  return records;
}

// A restart must not overwrite the files it continues from.
void check_restart_directory(const Config& config, const std::filesystem::path& output_dir) {
  if (config.derived.run_mode != 1) return;
  const std::filesystem::path dir = config.input.restart_dir;
  if (std::filesystem::exists(dir) && std::filesystem::exists(output_dir) &&
      std::filesystem::equivalent(dir, output_dir))
    throw InputError("RestartDir",
                     "the output directory must differ from RestartDir (the checkpoint would be overwritten)");
}

// The BK coefficient files of this run: native always, Fortran fort.21 on request.
struct BkOutput {
  BkFileWriter native;
  std::optional<FortranBkFileWriter> fortran;

  void append(double time, const BkCoefficients& bk) {
    native.append(time, bk);
    if (fortran) fortran->append(bk);
  }
};

// A propagation from tStart to tEnd (new run or restart) followed by the asymptotic analysis.
// Fortran: Main.f90 from SetupRDBasis to StateToStateProbs
void run_propagation(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                     const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                     const std::filesystem::path& output_dir) {
  check_restart_directory(config, output_dir);
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft);
  print_propagation_setup(out, style, grids, propagation);
  print_analysis_setup(out, style, config, propagation.analysis_basis(), propagation.analysis_energies());

  const BkRecords earlier =
      config.derived.run_mode == 1 ? restore_from_checkpoint(out, config, propagation) : BkRecords{};
  const bool write_native = config.input.checkpoint_format != "fortran";
  const bool write_fortran = config.input.checkpoint_format != "native";
  BkOutput bk_out{BkFileWriter(output_dir / bk_file, config, propagation.analysis_basis()), std::nullopt};
  if (write_fortran) {
    const int n_records = static_cast<int>(earlier.size()) + config.input.n_steps / config.derived.steps_per_analysis;
    bk_out.fortran.emplace(output_dir / fortran_bk_file, config, propagation.rotational_basis(),
                           propagation.analysis_basis(), gas, n_records);
  }
  for (const auto& [time, bk] : earlier) bk_out.append(time, bk);

  const AveragesTable table(out, style);
  table.print_header();
  std::vector<Averages> averages;
  const auto start = std::chrono::steady_clock::now();
  PropagationObserver observer;
  observer.on_averages = [&](const Averages& a) {
    table.print_row(a);
    averages.push_back(a);
  };
  observer.on_surface_fixed = [&](double z) {
    table.print_note("analysis surface on the specular grid fixed at Z = " + fixed(z, 6) + " bohr");
  };
  observer.on_analysis_step = [&](double time, const BkCoefficients& bk) { bk_out.append(time, bk); };
  propagation.run(observer);
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  out << "\n  Propagation time: " << fixed(elapsed.count(), 1) << " s\n";

  // Checkpoint for a restart. Fortran: SaveWaveFunc
  const AnalysisSurfaceState surface = propagation.analysis_surface().state();
  if (write_native)
    write_checkpoint(output_dir / checkpoint_file, config, config.input.t_end, propagation.wave_packet(), surface);
  if (write_fortran)
    write_fortran_wave_function(output_dir / fortran_wave_function_file, propagation.wave_packet(), surface);

  // Asymptotic analysis: S-matrix and probabilities at the analysis energies.
  const ScatteringResults results = propagation.scattering_results();
  print_scattering_results(out, style, config, propagation.analysis_basis(), results);

  write_averages_file(output_dir / "averages.dat", averages);
  const auto files = write_scattering_files(output_dir, config, propagation.analysis_basis(), results);
  out << "\n  Data files written to " << output_dir.string() << ": averages.dat";
  for (const auto& f : files) out << ", " << f.filename().string();
  out << ", " << bk_file;
  if (write_native) out << ", " << checkpoint_file;
  if (write_fortran) out << ", " << fortran_bk_file << ", " << fortran_wave_function_file;
  out << '\n';
}

}  // namespace

// Fortran: Main.f90 (new runs, restarts and potential-only runs).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color) {
  const Style style{color};
  try {
    if (command_line.threads) omp_set_num_threads(*command_line.threads);

    out << banner << '\n';
    out << style.emphasis("SPO-DVR (C++) version " SPODVR_VERSION) << "  -  " << omp_get_max_threads()
        << " OpenMP thread(s)\n";
    out << "  " << numerical_libraries_summary() << '\n';

    const Config config = build_config(read_input_file(command_line.input_file));
    print_input_summary(out, style, config);

    const auto potential = make_potential(config.input.system_id);
    const AngularGrid angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear,
                              config.input.j_in, config.input.mj_in);
    const SpatialGrids grids = make_spatial_grids(config);
    print_grids(out, style, config, grids, angular);

    if (config.derived.run_mode == 2)
      throw NotSupportedError("analysis-only runs (onlyAnalysis) are not implemented yet");

    PotentialGrid potential_grid = compute_potential_grid(*potential, config, grids, angular);
    print_asymptotic_potential(out, style, config, grids, potential_grid);

    const FftwFft fft;
    const LapackEigensolver eigensolver;
    const GasPhase gas = compute_gas_phase(config, grids, angular, potential_grid.asymptotic, fft, eigensolver);
    print_gas_phase(out, style, config, gas);

    std::filesystem::create_directories(command_line.output_dir);
    write_asymptotic_potential_file(command_line.output_dir / "asymptotic_potential.dat", grids, potential_grid);
    write_gas_phase_file(command_line.output_dir / "gas_phase_states.dat", grids, gas);
    out << "\n  Data files written to " << command_line.output_dir.string() << "\n";

    if (config.derived.run_mode == 6) {
      out << '\n' << style.emphasis("Potential-only run finished.") << '\n';
      return 0;
    }

    run_propagation(out, style, config, grids, angular, std::move(potential_grid), gas, fft, command_line.output_dir);
    out << '\n' << style.emphasis("Run finished.") << '\n';
    return 0;
  } catch (const InputError& e) {
    err << style.warning("Input error") << ": " << e.what() << '\n';
    return 2;
  } catch (const CheckpointError& e) {
    err << style.warning("Checkpoint error") << ": " << e.what() << '\n';
    return 4;
  } catch (const NotSupportedError& e) {
    err << style.warning("Not supported") << ": " << e.what() << '\n';
    return 3;
  } catch (const std::exception& e) {
    err << style.warning("Error") << ": " << e.what() << '\n';
    return 1;
  }
}

}  // namespace spodvr
