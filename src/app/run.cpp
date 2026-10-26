// Program driver: input -> grids -> potential -> gas-phase states -> propagation -> output.
// Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <string_view>
#include <vector>

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

// Names of the native checkpoint files in an output / restart directory.
constexpr std::string_view checkpoint_file = "checkpoint.bin";
constexpr std::string_view bk_file = "bk_coefficients.bin";

// Restart (nread = 1): reads the wave packet and analysis surface of the earlier run and copies
// its BK coefficients into the new BK file (they also enter the S-matrix).
// Fortran: WaveFunction.f90 ReadWaveFunc, AsympInfo.f90 ReadOrWriteAsympInfo (nread = 1),
// Analyzor.f90 GetSMatrixElemsOfOldRun
void restore_from_checkpoint(std::ostream& out, const Config& config, WavePacketPropagation& propagation,
                             BkFileWriter& bk_out) {
  const InputParameters& p = config.input;
  if (p.restart_format != "native") throw NotSupportedError("RestartFormat \"fortran\" is not implemented yet");
  const std::filesystem::path dir = p.restart_dir;
  const CheckpointInfo info = read_checkpoint(dir / checkpoint_file, config, propagation.wave_packet_for_restart());
  if (fortran_compare::rne(info.time, p.t_start))
    throw CheckpointError("tStart = " + full_precision(p.t_start) + " but the checkpoint '" +
                          (dir / checkpoint_file).string() + "' is at t = " + full_precision(info.time));
  propagation.restore_analysis_surface(info.surface);

  const auto last = read_bk_file(dir / bk_file, config, propagation.analysis_basis(),
                                 [&](double time, const BkCoefficients& bk) {
                                   propagation.add_earlier_analysis_step(time, bk);
                                   bk_out.append(time, bk);
                                 });
  if (last && fortran_compare::rne(*last, info.time))
    throw CheckpointError("the BK coefficient file '" + (dir / bk_file).string() + "' ends at t = " +
                          full_precision(*last) + ", the checkpoint at t = " + full_precision(info.time));
  out << "\n  Restarted from " << dir.string() << " at t = " << fixed(info.time, 2) << " a.u.\n";
}

// The BK file of this run must not overwrite the one being restarted from.
void check_restart_files_differ(const Config& config, const std::filesystem::path& output_dir) {
  if (config.derived.run_mode != 1) return;
  const std::filesystem::path in = std::filesystem::path(config.input.restart_dir) / bk_file;
  const std::filesystem::path out = output_dir / bk_file;
  if (std::filesystem::exists(in) && std::filesystem::exists(out) && std::filesystem::equivalent(in, out))
    throw InputError("RestartDir", "the output directory must differ from RestartDir (the checkpoint would be overwritten)");
}

// A propagation from tStart to tEnd (new run or restart) followed by the asymptotic analysis.
// Fortran: Main.f90 from SetupRDBasis to StateToStateProbs
void run_propagation(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                     const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                     const std::filesystem::path& output_dir) {
  if (config.input.checkpoint_format != "native")
    throw NotSupportedError("CheckpointFormat \"fortran\" / \"both\" is not implemented yet");
  check_restart_files_differ(config, output_dir);
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft);
  print_propagation_setup(out, style, grids, propagation);
  print_analysis_setup(out, style, config, propagation.analysis_basis(), propagation.analysis_energies());

  BkFileWriter bk_out(output_dir / bk_file, config, propagation.analysis_basis());
  if (config.derived.run_mode == 1) restore_from_checkpoint(out, config, propagation, bk_out);

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
  write_checkpoint(output_dir / checkpoint_file, config, config.input.t_end, propagation.wave_packet(),
                   propagation.analysis_surface().state());

  // Asymptotic analysis: S-matrix and probabilities at the analysis energies.
  const ScatteringResults results = propagation.scattering_results();
  print_scattering_results(out, style, config, propagation.analysis_basis(), results);

  write_averages_file(output_dir / "averages.dat", averages);
  const auto files = write_scattering_files(output_dir, config, propagation.analysis_basis(), results);
  out << "\n  Data files written to " << output_dir.string() << ": averages.dat";
  for (const auto& f : files) out << ", " << f.filename().string();
  out << ", " << bk_file << ", " << checkpoint_file << '\n';
}

}  // namespace

// Fortran: Main.f90 (new runs, restarts and potential-only runs).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color) {
  const Style style{color};
  try {
    if (command_line.threads) omp_set_num_threads(*command_line.threads);

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
